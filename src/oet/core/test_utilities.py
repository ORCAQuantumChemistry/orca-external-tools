"""
Utilities used in the test suite
"""

import multiprocessing as mp
import socket
import subprocess
import time
from collections.abc import Callable
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from multiprocessing.queues import Queue

WATER = [
    ("O", 0.0000, 0.0000, 0.0000),
    ("H", 0.2774, 0.8929, 0.2544),
    ("H", 0.6068, -0.2383, -0.7169),
]

OH = [
    ("O", 0.0000, 0.0000, 0.0000),
    ("H", 0.2774, 0.8929, 0.2544),
]


def read_result_file(filename: str | Path) -> tuple[int, float, list[float]]:
    """
    Reads the engrad file written by the wrapper

    Parameters
    ----------
    filename: str
        Name of the output file

    Returns
    -------
    int: number of atoms
    float: total energy
    grad: gradient list

    Raises
    ------
    OSError: Failing to open file
    ValueError: Failed to convert the input to int/float
    """
    with open(filename) as f:
        lines = f.readlines()

    # Remove comments from '#' to the end of line:
    data_lines = [li for line in lines if (li := line.partition("#")[0].strip())]

    # Extract data
    num_atoms = int(data_lines[0])
    energy = float(data_lines[1])
    gradients = [float(val) for val in data_lines[2:]]

    return num_atoms, energy, gradients


def write_input_file(
    filename: str | Path,
    xyz_filename: Path,
    charge: int,
    multiplicity: int,
    ncores: int,
    do_gradient: int | bool,
    pointcharges_filename: str | None = None,
) -> None:
    """
    Write an input file for the extopt wrapper script with the given parameters.

    Parameters
    ----------
    filename: str
        Output file name
    xyz_filename: str
        Filename of the structure file
    charge:int
        Molecular charge
    multiplicity: int
        Multiplicity
    ncores: int
        Number of cores to use
    do_gradient: int
        Compute gradient (1) or not (0)
    pointcharges_filename: str | None = None
        optional filename of the pointcharges
    """

    # Validate inputs (basic checks)
    if xyz_filename.suffix != ".xyz":
        raise ValueError("xyz_filename did not end with '.xyz'")
    if multiplicity <= 0:
        raise ValueError("multiplicity must be a positive integer")
    if ncores <= 0:
        raise ValueError("ncores must be a positive integer")
    if type(do_gradient) is bool:
        do_gradient = int(do_gradient)
    if do_gradient not in (0, 1):
        raise ValueError("do_gradient must be 0 or 1")

    with open(filename, "w") as f:
        f.write(f"{xyz_filename} # xyz filename: string, ending in '.xyz'\n")
        f.write(f"{charge} # charge: integer\n")
        f.write(f"{multiplicity} # multiplicity: positive integer\n")
        f.write(f"{ncores} # NCores: positive integer\n")
        f.write(f"{do_gradient} # do gradient: 0 or 1\n")
        if pointcharges_filename:
            f.write(f"{pointcharges_filename} # point charge filename: string (optional)\n")
        else:
            f.write("\n")  # Write a blank line if no point charges file given


def write_xyz_file(filename: str | Path, atoms: list[tuple[str, float, float, float]]) -> None:
    """
    Write a file with the given format:

    Parameters
    ----------
    filename: str
        Output file name
    atoms: list[tuple[str, float, float, float]]
        atomic symbols and positions [(symbol, x, y, z), ...]
    """
    with open(filename, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for atom in atoms:
            symbol, x, y, z = atom
            f.write(f"{symbol} {x:.4f} {y:.4f} {z:.4f}\n")


def run_wrapper(
    inputfile: str | Path,
    script_path: str | Path,
    outfile: str | Path,
    args: list[str] | None = None,
    timeout: float | None = 10.0,
) -> None:
    """
    Run the wrapper

    Parameters
    ----------
    inputfile: str | Path
        Inputfile
    script_path: str | Path
        Path to the oet script
    outfile: str | Path
        File to write the output to
    args: list[str] | None, default = None
        Additional arguments
    timeout: float | None, default: 10 s
        Default timeout time (seconds)
    """
    cmd = [script_path, inputfile]
    if args:
        cmd += args

    with open(outfile, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=timeout, check=False)


def add_arguments(args: str | list[str], additions: list[str]) -> list[str]:
    """
    Add arguments

    Parameters
    ----------
    args: str | list[str]
        Arguments that should be extended
    additions: list[str]
        Arguments to add

    Returns
    -------
    list[str]: extended arguments
    """
    if isinstance(args, str):
        args = [args]
    args += additions
    return args


def get_filenames(basename: str) -> tuple[Path, Path, Path, Path]:
    """
    Set the filenames according to how ORCA would do and cleans any input existing
    """
    xyz_file = Path(basename + ".xyz")
    input_file = Path(basename + ".extinp.tmp")
    engrad_out = Path(basename + ".engrad")
    output_file = Path(basename + ".out").resolve()
    clear_files(basename=basename)
    return xyz_file, input_file, engrad_out, output_file


def clear_files(basename: str) -> None:
    """
    Remove every file starting with basename
    """
    dir_path = Path.cwd()
    for f in dir_path.glob(basename + "*"):
        if f.is_file():
            f.unlink()  # remove file


def _worker(
    fn: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any], q: "Queue[bool]"
) -> None:
    """
    Helper for executing a function.

    Parameters
    ----------
    fn: Callable[..., T]
        Callable function that is executed.
    args: tuple[Any]
        Any positional arguments that are given to the function call.
    kwargs: dict[str, Any]
        Any keyword arguments that are given to the function call.
    q: Queue[bool]
        Queue used to put the function call.
    """
    try:
        # Call the function
        _ = fn(*args, **kwargs)
        # Don't check what it did, just return ok, if the function didn't crash
        q.put(True)
    except Exception:
        q.put(False)


def wait_for_server(
    process: subprocess.Popen[str],
    ip_port: str,
    timeout: float = 60.0,
    poll_interval: float = 0.1,
) -> None:
    """
    Wait until a server process accepts TCP connections.

    Parameters
    ----------
    process: subprocess.Popen[str]
        The subprocess used for starting the server.
    id_port: str
        The server address.
    timeout: float, default: 60.0
        The allowed time for waiting.
    poll_interval: float, default: 0.1
        The interval for pinging the server.
    """

    # Get the server address
    host, port_str = ip_port.rsplit(":", 1)
    port = int(port_str)

    # Track the time
    start_time = time.monotonic()
    while time.monotonic() - start_time < timeout:
        # Fail early if the server process already died.
        returncode = process.poll()
        if returncode is not None:
            raise RuntimeError(f"Server terminated unexpectedly with return code {returncode}.")

        # Try pinging the server
        try:
            with socket.create_connection(
                (host, port),
                timeout=poll_interval,
            ):
                return
        except OSError:
            pass

        time.sleep(poll_interval)

    raise TimeoutError(f"Server did not become ready within {timeout:.0f} s.")


class TimeoutCallError(StrEnum):
    """Possible errors that are returned by TimeoutCall"""

    # Function timed out
    TIMEOUT = "timeout"
    # Function crashed
    CRASH = "crash"
    # General error
    ERROR = "error"


class TimeoutCall:
    """
    Class for calling a function with a certain timeout.
    Useful for functions that, e.g., download files.
    Doesn't return the result of the function as it might not be pickled.
    """

    def __init__(self, fn: Callable[..., Any]) -> None:
        """
        Initialization of the class.

        Parameters
        ----------
        fn: Callable[..., T]
            Callable function that is executed.
        """
        self.fn = fn
        self.timed_out = False

    def __call__(
        self, *args: Any, timeout: float = 10, **kwargs: Any
    ) -> tuple[bool, TimeoutCallError | None]:
        """
        Execute the function set in __init__ with the timeout defined there.

        Parameters
        ----------
        args: Any
            Any positional arguments that are given to the function call.
        timeout: float, default: 10 sec.
            Timeout in sec.
        kwargs: Any
            Any keyword arguments that are given to the function call.

        Returns
        -------
        bool
            True, if everything was ok. False otherwise.
        TimeoutCallError | None
            Either the error type if failed or None.
        """

        # Start process and wait the timeout
        q: Queue[bool] = mp.Queue()
        p: mp.Process = mp.Process(target=_worker, args=(self.fn, args, kwargs, q))
        p.start()
        p.join(timeout)

        # Check if there was any error
        if p.exitcode not in (0, None):
            return False, TimeoutCallError.CRASH

        # Check if the process is still alive. If yes, it has timed out.
        if p.is_alive():
            p.terminate()
            p.join()
            return False, TimeoutCallError.TIMEOUT

        # Check if the worker provides the correct result
        try:
            status_ok = q.get(timeout=1)
        except Exception:
            return False, TimeoutCallError.CRASH

        # Check if there was a general error
        if not status_ok:
            return False, TimeoutCallError.ERROR

        # If everything went well, return True and no Error.
        return True, None
