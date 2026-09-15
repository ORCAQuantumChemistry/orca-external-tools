#!/usr/bin/env python3
"""
Provides the class for running a server loaded with a calculator

Provides
--------
class: OtoolServer
    Main class for running a server
function: main
    Main function for executing
"""

from __future__ import annotations

import argparse
import concurrent
import concurrent.futures
import contextlib
import gc
import importlib
import io
import logging
import multiprocessing as mp
import os
import signal
import sys
import threading
import time
import traceback
import typing
from argparse import Action, ArgumentParser
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from concurrent.futures import BrokenExecutor, ProcessPoolExecutor
from contextlib import redirect_stdout
from pathlib import Path
from types import FrameType
from typing import Any

import psutil
from flask import Flask, Response, jsonify, request
from waitress import serve

from oet import __version__ as version
from oet.core.base_calc import CALCULATOR_CLASSES, BaseCalc
from oet.core.misc import get_ncores_from_input

if typing.TYPE_CHECKING:
    from argparse import Namespace


# Cache of initialized calculators: each worker process populates its own copy of the cache
# This guarantees that for calculators which are inadvertently not thread-safe
# (because they temporarily hold calculation data in member variables)
# we never use the same instance by multiple processes in parallel.
# key: (module, class, frozenset(setup_items)) -> calc instance
# OrderedDict, where the most currently used entry is moved to the end
_WORKER_CALC_CACHE: OrderedDict[tuple[str, str, frozenset[tuple[str, Any]]], Any] = OrderedDict()

# Logger
logger = logging.getLogger(__name__)


def _pop_one_worker(protected_key: tuple[str, str, frozenset[tuple[str, Any]]] | None) -> bool:
    """
    Removes left most worker that isn't protected

    Parameters
    ----------
    protected_key: tuple[str, str, frozenset[tuple[str, Any]]] | None
        Protected keys that should not be deleted.

    Returns
    -------
    bool
        Whether a worker could be deleted or not.
    """
    for k in list(_WORKER_CALC_CACHE.keys()):
        if protected_key is not None and k == protected_key:
            # Move protected to MRU end so we don't keep looping on it
            _WORKER_CALC_CACHE.move_to_end(k)
            continue
        _WORKER_CALC_CACHE.pop(k, None)
        logger.debug(f"PID {os.getpid()}: Deleted calculator with id: {k}")
        return True
    return False


def _evict_until_within_limits(
    mem_limit_mib: int, protected_key: tuple[str, str, frozenset[tuple[str, Any]]] | None
) -> None:
    """
    Evict least-recently-used calculators until the memory limit is satisfied.
    `protect_key` (if given) will not be evicted (skip over it).

    Parameters
    ----------
    mem_limit_mib: int
        Maximum memory usable by server in MiB
    protected_key: tuple[str, str, frozenset[tuple[str, Any]]]
        The one that should not be deleted from the _WORKER_CALC_CACHE
    """
    # RSS-based eviction (only if we can read RSS and a limit is set)
    # Get memory usage in megabyte
    rss = psutil.Process(os.getpid()).memory_info().rss / (1024**2)
    # Return if memory can't be measured
    if rss is None:
        print("No memory use could be detected. Take care, memory usage might stack.")
        return
    # Evict until RSS is below budget (best effort; GC happens after pops)
    # Add a small grace (5%) to avoid thrashing around the threshold.
    grace = int(mem_limit_mib * 0.05)
    target = mem_limit_mib - grace
    # Avoid negative target
    target = max(target, int(mem_limit_mib * 0.5))
    logger.debug(f"PID {os.getpid()}: Memory use before: {rss} / {target}")
    # Try evicting a few items and rechecking RSS
    # (popping, then letting GC run once in a while)
    attempts = 0
    while rss > target and _WORKER_CALC_CACHE:
        attempts += 1
        if not _pop_one_worker(protected_key=protected_key):
            break
        # Run garbage collector
        gc.collect()
        # Refresh current memory usage
        rss = psutil.Process(os.getpid()).memory_info().rss / (1024**2)
    logger.debug(f"PID {os.getpid()}: Memory use after: {rss} / {target}")


class CalculatorRuntimeException(RuntimeError):
    """Custom exception class to pass STDOUT to the client along with any runtime error."""

    def __init__(self, stdout: str):
        self.stdout = stdout


def _run_calc_in_process(
    calc_module: str, calc_class: str, run_kwargs: dict[str, Any], max_memory_per_thread: int
) -> str:
    """
    Worker entrypoint. Runs in a separate process.
    We lazily create & cache a calculator instance per unique setup.

    Parameters
    ----------
    calc_module: str
        Module to load for calculator
    calc_class: str
        Calculator type
    run_kwargs: dict
        Infos from client about the run
    max_memory_per_thread: int
        Maximum memory in MiB

    Returns
    -------
    str
        The STDOUT of the calculator's `run` function
    """

    # Make run_kwargs with only method-specific settings
    # It is assumed that all not parsed arguments are also method specific
    args_parsed_frozen = tuple(sorted(vars(run_kwargs["args_parsed"]).items()))
    args_not_parsed_frozen = tuple(run_kwargs["args_not_parsed"])
    method_specific_args = {
        "args_parsed": args_parsed_frozen,
        "args_not_parsed": args_not_parsed_frozen,
    }
    # Use those to check if a calculator with these settings already exists.
    # If not, make another one and delete old calculators if it would exceed the current memory limit.
    key = (calc_module, calc_class, frozenset(method_specific_args.items()))
    calc = _WORKER_CALC_CACHE.get(key)
    if calc is None:
        mod = importlib.import_module(calc_module)
        Cls = getattr(mod, calc_class)
        calc = Cls()
        _WORKER_CALC_CACHE[key] = calc
        logger.debug(f"PID {os.getpid()}: Initialized new calculator with id: {key}")
    else:
        logger.debug(f"PID {os.getpid()}: Using existing calculator with id: {key}")
    # Mark as most recently used (move it to the end of ordered dict)
    _WORKER_CALC_CACHE.move_to_end(key)

    # Evict old entries if max memory is used, but never evict the one that is used
    _evict_until_within_limits(max_memory_per_thread, protected_key=key)
    # Run calc.run and return its STDOUT
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            calc.run(**run_kwargs)
    except Exception as e:
        # Attach STDOUT to the exception
        raise CalculatorRuntimeException(buf.getvalue()) from e
    return buf.getvalue()


class CoreLimiter:
    """
    Enforces a global core budget across concurrent jobs.
    Blocks until enough cores are free.
    """

    def __init__(self, total_cores: int) -> None:
        # Total cores of server
        self.total = int(total_cores)
        # Available cores of server
        self.available = int(total_cores)
        # Waiting
        self._cv = threading.Condition()

    def acquire(self, n: int) -> None:
        """
        Wait until the required cores for the job are available again

        Parameters
        ----------
        n: int
            Number of cores that are requested for the job
        """
        n = int(n)
        with self._cv:
            if n > self.total:
                # Fail fast: this job can never run on this server
                raise ValueError(f"Requested {n} cores but only {self.total} total available.")
            while n > self.available:
                self._cv.wait()
            self.available -= n

    def release(self, n: int) -> None:
        """
        Releases number of cores after the job is done

        Parameters
        ----------
        n: int
            Number of cores to release
        """
        n = int(n)
        with self._cv:
            self.available += n
            self.available = min(self.available, self.total)
            self._cv.notify_all()


class CalculatorClass:
    """
    Stores the calculator class and is used to build the full parser multiple times
    """

    def __init__(self, calc_type: str) -> None:
        """
        Set the calculator types and import it

        Parameters
        ----------
        calc_type: str
            Type of calculator
        """
        # Determine module and calculator name to be imported from dict
        import_module, calculator_name = CALCULATOR_CLASSES[calc_type]
        # Module to import from
        self.import_module = import_module
        # Name of calculator class
        self.calculator_name = calculator_name
        # Import module
        mod = importlib.import_module(import_module)
        # Set calculator type
        self._cls: type[BaseCalc] = getattr(mod, calculator_name)

    def build_full_parser(self, parser: ArgumentParser) -> None:
        """
        Extends the given argument parser by calling the class method `extend_parser` of the calculator class

        parser: ArgumentParser
            Parser to be extended
        """
        if self._cls is None:
            raise RuntimeError("Call set_calculator_type() first.")
        self._cls.extend_parser(parser)


class OtoolServer:
    def __init__(
        self,
        calc_class: CalculatorClass,
        total_cores: int,
        executor: ProcessPoolExecutor,
        max_memory_per_thread: int,
    ):
        # Calculator class used to build parser
        self.calc_class = calc_class
        # Core limiter to prevent using more cores than the server has
        self.core_limiter = CoreLimiter(total_cores)
        # Executor to run the actual calculation
        self.executor = executor
        # Mechanism to handle a broken executor
        self._executor_lock = threading.Lock()
        self._executor_broken = False
        # Maximum memory of the server in MiB
        self.max_memory_per_thread = max_memory_per_thread

    def handle_client(self, content: Mapping[str, Any]) -> dict[str, Any]:
        """
        Takes the input from client and runs the calculation

        Parameters
        ----------
        content: Mapping[str, Any]
            Message received from client

        Returns
        -------
        dict[str, Any]
            Message to be sent to client
        """
        # First, check if the executor was broken when handling another request and if so, bail out directly
        with self._executor_lock:
            if self._executor_broken:
                raise RuntimeError("Server unavailable")

        arguments: Sequence[str] = content["arguments"]
        working_dir = Path(content["directory"]).resolve()

        # Parse client args
        inputfile, args, args_not_parsed = self.parse_client_input(arguments)

        # Make inputfile absolute (per-request, thread-safe)
        inputfile_path = (working_dir / inputfile).resolve()

        # Get per-job core demand
        ncores_job = get_ncores_from_input(inputfile_path)

        # Gate by cores
        self.core_limiter.acquire(ncores_job)
        try:
            # Submit the job to a separate process
            run_kwargs = {
                "inputfile": str(inputfile_path),
                "args_parsed": args,
                "args_not_parsed": args_not_parsed,
                "directory": str(working_dir),
            }
            fut = self.executor.submit(
                _run_calc_in_process,
                self.calc_class.import_module,
                self.calc_class.calculator_name,
                run_kwargs,
                self.max_memory_per_thread,
            )
            # Will raise if the worker raised
            output = fut.result()
        except BrokenExecutor as e:
            # If a worker process gets killed, e.g. by the OS's OOM killer, the process pool gets broken
            # and any subsequent requests will hang. The executor must be stopped and possibly restarted.
            # Mark it as broken to prevent other requests from hanging.
            with self._executor_lock:
                if not self._executor_broken:
                    self._executor_broken = True

                    # Since new requests are also likely to get killed, just kill the server completely.
                    # No need to do it asynchronously, so as not to leave and client requests hanging.
                    def _shutdown() -> None:
                        # Short delay to flush responses
                        time.sleep(0.5)
                        os.kill(os.getpid(), signal.SIGTERM)

                    threading.Thread(target=_shutdown, daemon=True).start()
            # Pass the exception up and to the client
            raise RuntimeError(
                "A worker process was terminated unexpectedly: server shutting down"
            ) from e
        except:
            raise
        else:
            return {"status": "Success", "stdout": output}
        finally:
            self.core_limiter.release(ncores_job)

    def parse_client_input(self, arguments: Sequence[str]) -> tuple[str, Namespace, list[str]]:
        """
        Handles the input sent by client

        Parameters
        ----------
        arguments: list[str]
            Arguments from client

        Returns
        -------
        str
            inputfile name
        Namespace
            parsed settings
        list[str]
            not parsed settings
        """
        # First server-related settings
        parser = ArgumentParser(
            prog="oet_server",
            description="Client arguments parser.",
        )
        parser.add_argument("inputfile")

        # Let the calculator define its per-request flags
        self.calc_class.build_full_parser(parser)

        args, remaining_args = parser.parse_known_args(arguments)

        inputfile = args.inputfile
        delattr(args, "inputfile")

        return inputfile, args, remaining_args

def create_app(server: OtoolServer) -> Flask:
    """
    Takes the OtoolServer and returns a Flask application
    """
    app = Flask(__name__)

    @app.get("/healthz")
    def healthz() -> Response:
        return jsonify({"status": "OK"})

    @app.post("/calculate")
    def calculate() -> Response:
        try:
            data = request.get_json(force=True, silent=False)
            if not isinstance(data, dict):
                return jsonify(
                    {
                        "status": "Error",
                        "error_message": "Invalid JSON payload",
                        "error_type": "ValueError",
                    }
                )

            arguments = data.get("arguments")
            directory = data.get("directory")

            if not isinstance(arguments, list) or not isinstance(directory, str):
                return jsonify(
                    {
                        "status": "Error",
                        "error_message": "Payload must have list 'arguments' and str 'directory'",
                        "error_type": "ValueError",
                    }
                )

            # Validate directory exists and is a dir
            p = Path(directory)
            if not p.exists() or not p.is_dir():
                return jsonify(
                    {
                        "status": "Error",
                        "error_message": f"Invalid directory: {p}",
                        "error_type": "ValueError",
                    }
                )

            # Delegate to server
            result = server.handle_client({"arguments": arguments, "directory": directory})
            return jsonify(result)

        except Exception as e:
            output = {
                "status": "Error",
                "error_message": str(e),
                "error_type": type(e).__name__,
                "traceback": traceback.format_exc(),
            }
            # attach STDOUT if possible
            if isinstance(e, CalculatorRuntimeException):
                output["stdout"] = e.stdout
            return jsonify(output)

    return app


def get_available_methods() -> list[str]:
    """
    Attempt to initialize all calculators in `CALCULATOR_CLASSES`
    and return the keys for those that don't throw an error.
    Note that this function takes a few seconds due to heavy imports.
    """
    available = []
    # discard stdour and stderr
    with (
        open(os.devnull, "w") as devnull,
        contextlib.redirect_stdout(devnull),
        contextlib.redirect_stderr(devnull),
    ):
        for method in CALCULATOR_CLASSES:
            with contextlib.suppress(Exception):
                CalculatorClass(method)
                available.append(method)
    return available


class PrintAvailableMethods(Action):
    """Custom parser action to print all available methods and exit"""

    def __call__(
        self,
        parser: ArgumentParser,
        _: Namespace,
        __: str | Sequence[Any] | None,
        ___: str | None = None,
    ) -> None:
        parser.exit(0, "\n".join(get_available_methods()) + "\n")


def cleanup_and_exit(
    signum: int,
    _: FrameType | None,
    executor: concurrent.futures.ProcessPoolExecutor,
    parser: argparse.ArgumentParser,
) -> None:
    """
    Handle termination signals by shutting down the executor and exiting the process.

    Parameters
    ----------
    signum : int
        The signal number received by the handler (e.g., ``signal.SIGTERM``).
    frame : FrameType or None
        The current stack frame at the moment the signal was received, supplied
        automatically by the Python signal machinery. Typically unused.
    executor : concurrent.futures.ProcessPoolExecutor
        The process pool executor whose worker processes should be terminated
        cleanly before exiting.
    parser : argparse.ArgumentParser
        The argument parser whose ``exit()`` method is used to terminate the program.
    """
    print(f"PID {os.getpid()} received signal {signum}, shutting down...")
    executor.shutdown(wait=False, cancel_futures=True)
    parser.exit(0)


def worker_initializer() -> None:
    """Initialize a worker process by restoring default signal handling."""
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    if sys.platform != "win32":
        signal.signal(signal.SIGHUP, signal.SIG_DFL)


def main() -> None:
    """
    Main routine of oet_server
    """
    # First parse arguments
    parser = ArgumentParser(
        prog="oet_server",
        description="Starts a server with the selected method that performs the energy/gradient calculation.",
        epilog="Specific keywords of the calculators should also set here.",
    )
    parser.add_argument("--version", action="version", version=version)
    parser.add_argument(
        "method",
        choices=CALCULATOR_CLASSES.keys(),
        help="Type of calculation to execute.",
    )
    # Server related things already here to have it in help
    parser.add_argument(
        "-b",
        "--bind",
        metavar="hostname:port",
        default="127.0.0.1:8888",
        dest="host_port",
        help="Server bind address and port. Default: 127.0.0.1:8888.",
    )
    parser.add_argument(
        "-n",
        "--nthreads",
        metavar="nthreads",
        type=int,
        default=1,
        dest="nthreads",
        help="Number of threads to use. Default: 1",
    )
    parser.add_argument(
        "-l",
        "--list-methods",
        nargs=0,
        action=PrintAvailableMethods,
        help="List available calculators and exit. Note that this takes a few seconds.",
    )
    parser.add_argument(
        "-m",
        "--memory-per-thread",
        metavar="memory_per_thread",
        type=int,
        default=500,
        dest="memory_per_thread",
        help="Maximum memory per thread in mebibyte(MiB). Default 500.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="Print additional output.",
    )

    # First parse only the known args to get method/nthreads/etc.
    args, ignored_args = parser.parse_known_args()

    # Logging for printout of infos
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    # Suppress the warnings (e.g. jobs are queued)
    logging.getLogger("waitress.queue").setLevel(logging.ERROR)

    if ignored_args:
        logger.warning("The following arguments will be ignored: " + " ".join(ignored_args))

    signal.signal(signal.SIGTERM, lambda s, f: cleanup_and_exit(s, f, executor, parser))
    signal.signal(signal.SIGINT, lambda s, f: cleanup_and_exit(s, f, executor, parser))
    if sys.platform != "win32":
        signal.signal(signal.SIGHUP, lambda s, f: cleanup_and_exit(s, f, executor, parser))

    # Create workers
    workers = args.nthreads
    # Use spawn start method so CUDA can be initialized safely in worker processes
    mp_ctx = mp.get_context("spawn")
    # Initialize the ProcessPool
    executor = ProcessPoolExecutor(
        max_workers=workers, initializer=worker_initializer, mp_context=mp_ctx
    )

    # Make a CalculatorClass getting the hooks on calculators argument parsing
    # Info on calculator type is store in the object for client requests
    calcClass = CalculatorClass(args.method)

    # Then initialize a server instance that uses the calc_class
    server = OtoolServer(
        calc_class=calcClass,
        total_cores=args.nthreads,  # or another CLI flag like --total-cores
        executor=executor,
        max_memory_per_thread=args.memory_per_thread,
    )

    # Start the server
    host, port = args.host_port.split(":")
    app = create_app(server)
    # Waitress will create a thread per incoming request (bounded by 'threads').
    # Each request acquires/releases a calculator safely.
    serve(app, host=host, port=int(port), threads=args.nthreads)


if __name__ == "__main__":
    main()
