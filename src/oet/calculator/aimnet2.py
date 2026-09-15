#!/usr/bin/env python3
"""
Calculator for using AIMNet2 (https://github.com/isayevlab/aimnetcentral),
compatible with ORCA's ExtTool interface.

Provides
--------
class: Aimnet2Calc(CalcServer)
    Class for performing a AIMNet2 calculation together with ORCA
main: function
    Main function
"""

import sys
import warnings
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any, ClassVar

from requests.exceptions import HTTPError

from oet import ASSETS_DIR
from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import (
    ENERGY_CONVERSION,
    LENGTH_CONVERSION,
    resolve_model_file,
    xyzfile_to_at_coord,
)

try:
    from aimnet.calculators import AIMNet2Calculator
    from aimnet.calculators.model_registry import get_model_path, load_model_registry
except ImportError as err:
    print(
        f"[MISSING] Required module aimnet not found: {err}.\n"
        "Please install the packages in the virtual environment.\n"
    )
    sys.exit(1)
try:
    import torch
except ImportError as e:
    print("[MISSING] torch not found:", e)
    sys.exit(1)


DEFAULT_MODEL_PATH = ASSETS_DIR / "aimnet2"

# Periodic table for symbol→Z conversion.
# Per-model element rejection is delegated to AIMNet2Calculator.eval()
# against model.metadata["implemented_species"].
_PERIODIC_TABLE: tuple[str, ...] = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)
assert len(_PERIODIC_TABLE) == 118, f"expected 118 elements, got {len(_PERIODIC_TABLE)}"
_SYMBOL_TO_Z: dict[str, int] = {sym: i + 1 for i, sym in enumerate(_PERIODIC_TABLE)}


class Aimnet2Calc(BaseCalc):
    """ORCA ExtTool calculator backed by AIMNet2 (aimnet>=0.2,<0.3)."""

    # Wired into extend_parser's choices=, single source of truth.
    _SUPPORTED_DEVICES: tuple[str, ...] = ("cpu", "cuda", "auto")
    _COULOMB_METHODS: tuple[str, ...] = ("simple", "dsf", "ewald")
    # Single source of truth for tri-state choices AND value translation.
    _TRISTATE_MAP: ClassVar[dict[str, bool | None]] = {
        "auto": None,
        "on": True,
        "off": False,
    }

    _calc: AIMNet2Calculator | None = None
    _setup_args: frozenset[tuple[str, Any]] | None = None

    def get_calculator(self) -> AIMNet2Calculator:
        """
        Returns AIMNet2 calculator

        Returns
        -------
        AIMNet2Calculator
            AIMNet2 calculator
        """
        return self._calc

    @staticmethod
    def _aimnet_alias_resolver(name: str) -> str | None:
        """Look up `name` in aimnet's model registry; return canonical
        file name if it is an alias, else None."""
        registry = load_model_registry()
        aliases = registry["aliases"]
        canonical = aliases.get(name) or name
        file = registry["models"].get(canonical, {}).get("file")
        return str(file) if file else None

    @staticmethod
    def _aimnet_subdir_fallback(name: str, exc: Exception) -> str | None:
        """If aimnet's URL fetch 404s on a flat name, retry once under
        the registry's `<family>/<name>` subdirectory layout (e.g.
        `aimnet2_wb97m_1` is actually stored at `aimnet2/aimnet2_wb97m_1`)."""
        if isinstance(exc, HTTPError) and "/" not in name:
            subdir_name = name.split("_")[0] + "/" + name
            print(
                f'Failed to find model "{name}" at URL: {exc.response.url}\n'
                f'Trying again with model name "{subdir_name}"',
                file=sys.stderr,
            )
            return subdir_name
        return None

    @classmethod
    def get_model_file(cls, model: str, model_dir: str) -> Path:
        """
        Make sure model file exists in the correct location.
        If `model` is an absolute path, it must already exist.
        Otherwise, let AIMNet2 download it, then move it to `model_dir`.

        Thin wrapper that supplies aimnet-specific alias / fetch / fallback
        callables to the calculator-agnostic
        :func:`oet.core.misc.resolve_model_file`.

        Parameters
        ----------
        model : str
            Model name, e.g. "aimnet2", or filename, e.g.
            "aimnet2-wb97m-d3_0.pt", or absolute path.
        model_dir : str
            Directory to look for or store the model file.

        Returns
        -------
        Path
            Full path to the model file.

        Raises
        ------
        FileNotFoundError
            If the model file is given by absolute path and does not exist.
        FileExistsError
            If `model_dir` exists but is not a directory, or `cached_path`
            exists but is not a file.
        """
        return resolve_model_file(
            model,
            model_dir,
            fetch=get_model_path,
            alias_resolver=cls._aimnet_alias_resolver,
            fetch_fallback=cls._aimnet_subdir_fallback,
        )

    def setup(
        self,
        model: str,
        model_dir: str,
        device: str | None,
        ncores: int,
        *,
        compile_model: bool = False,
        nb_threshold: int = 120,
        ensemble_member: int = 0,
        coulomb: str = "auto",
        dispersion: str = "auto",
        coulomb_method: str | None = None,
        coulomb_cutoff: float = 15.0,
        dftd3_cutoff: float | None = None,
        dftd3_smoothing_fraction: float | None = None,
    ) -> None:
        """Construct the calculator + apply post-ctor configuration.

        First call wins: subsequent calls with the same args are a no-op;
        calls with different args raise (server cache is responsible for
        keying on args, so a single Aimnet2Calc instance only ever sees
        one configuration). Per-call work belongs in run_aimnet2.

        Parameters
        ----------
        model : str
            Model name (e.g. "aimnet2"), filename, or absolute path. Forwarded
            to :meth:`get_model_file` for resolution and caching.
        model_dir : str
            Local cache directory for downloaded model files.
        device : str | None
            Compute device. One of {"cpu", "cuda", "auto"}; "auto" maps to
            None upstream (auto-detect). When "cuda" is requested, this
            method raises if torch.cuda.is_available() is False.
        ncores : int
            Number of CPU threads. Sets torch.set_num_threads process-wide
            once per worker; does not change per call.
        compile_model : bool, default: False
            If True, wrap the model with torch.compile. Server mode only;
            see readmes/aimnet2.md.
        nb_threshold : int, default: 120
            Adaptive neighbor-list batch size forwarded to the constructor.
        ensemble_member : int, default: 0
            Which ensemble member (0..3) to load. OET runs ONE model per call.
        coulomb : str, default: "auto"
            Tri-state ("auto" / "on" / "off") forcing the long-range Coulomb
            module on or off. "auto" defers to the model's default.
        dispersion : str, default: "auto"
            Tri-state for the DFT-D3 module, same semantics as `coulomb`.
        coulomb_method : str | None, default: None
            One of {"simple", "dsf", "ewald"} or None. If set, applied via
            `set_lrcoulomb_method` after construction.
        coulomb_cutoff : float, default: 15.0
            Coulomb cutoff in Ångström (used by dsf/ewald). For
            `aimnet2-rxn*` models the training-frozen value is 4.6;
            other values trigger a UserWarning.
        dftd3_cutoff : float | None, default: None
            Override the D3 cutoff via `set_dftd3_cutoff` if non-None.
        dftd3_smoothing_fraction : float | None, default: None
            Override the D3 smoothing fraction via `set_dftd3_cutoff` if
            non-None.

        Raises
        ------
        RuntimeError
            If `device` is "cuda" but no CUDA device is available, or if
            this method has already been called with a different argument
            set on this instance.
        """
        # Validate device availability BEFORE the cache short-circuit so
        # a cached CPU calc does not silently service a --device cuda
        # request on a no-CUDA box.
        if device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but not available")

        # Normalize device for storage so "auto" and None compare equal
        # in the args-match check below (both resolve to None upstream).
        device_arg = None if device == "auto" else device

        # Args-match short-circuit (uses normalized device).
        # frozenset of items (not tuple) so adding a new kwarg doesn't
        # require editing positional indices everywhere.
        new_args = frozenset(
            {
                "model": model,
                "model_dir": model_dir,
                "device": device_arg,
                "ncores": ncores,
                "compile_model": compile_model,
                "nb_threshold": nb_threshold,
                "ensemble_member": ensemble_member,
                "coulomb": coulomb,
                "dispersion": dispersion,
                "coulomb_method": coulomb_method,
                "coulomb_cutoff": coulomb_cutoff,
                "dftd3_cutoff": dftd3_cutoff,
                "dftd3_smoothing_fraction": dftd3_smoothing_fraction,
            }.items()
        )
        if self._calc is not None:
            if new_args == self._setup_args:
                return
            raise RuntimeError(
                "Aimnet2Calc.setup() called with different args than the "
                "cached calculator. Server-mode callers must key the cache "
                "on setup args."
            )

        # Validate string-choice args with helpful error messages
        # (argparse normally enforces this for CLI users, but direct-API
        # callers benefit from the early check too).
        for arg_name, arg_val, allowed in (
            ("coulomb", coulomb, tuple(self._TRISTATE_MAP)),
            ("dispersion", dispersion, tuple(self._TRISTATE_MAP)),
        ):
            if arg_val not in allowed:
                raise ValueError(f"setup(): {arg_name}={arg_val!r} not in {allowed}")
        if coulomb_method is not None and coulomb_method not in self._COULOMB_METHODS:
            raise ValueError(
                f"setup(): coulomb_method={coulomb_method!r} not in {self._COULOMB_METHODS}"
            )

        # When coulomb_method=='simple', the cutoff is meaningless (simple has
        # no cutoff in the upstream LR module). Warn loudly so users don't
        # think a custom value applies.
        if coulomb_method == "simple" and coulomb_cutoff != 15.0:
            warnings.warn(
                f"--coulomb-cutoff={coulomb_cutoff} is ignored when "
                "--coulomb-method=simple (simple Coulomb has no cutoff; "
                "cutoff applies only to dsf/ewald methods).",
                UserWarning,
                stacklevel=2,
            )

        model_path = str(self.get_model_file(model, model_dir))

        # Warn if user picked aimnet2-pd: it carries baked-in CPCM/THF
        # solvation, so energies are NOT gas-phase. Easy to miss from the
        # README alone.
        model_path_lower = Path(model_path).stem.lower()
        if "aimnet2-pd" in model_path_lower or "aimnet2_pd" in model_path_lower:
            warnings.warn(
                "aimnet2-pd embeds B97-3c with implicit CPCM/THF solvation; "
                "energies are NOT gas-phase. Do not mix aimnet2-pd energies "
                "with any other family.",
                UserWarning,
                stacklevel=2,
            )

        # Warn if user picked an aimnet2-rxn family model with a
        # non-trained --coulomb-cutoff override. The rxn family was
        # trained with a 4.6 A long-range cutoff; other values produce
        # physically-suspect electrostatics. Only relevant when the user
        # explicitly enabled the long-range method.
        is_rxn_family = "aimnet2-rxn" in model_path_lower or "aimnet2_rxn" in model_path_lower
        if is_rxn_family and coulomb_method is not None and coulomb_cutoff != 4.6:
            warnings.warn(
                f"--coulomb-cutoff={coulomb_cutoff} with -m {model} disagrees with "
                "the aimnet2-rxn training cutoff (4.6 A). Pass --coulomb-cutoff 4.6 "
                "for physically meaningful long-range electrostatics.",
                UserWarning,
                stacklevel=2,
            )

        self._calc = AIMNet2Calculator(
            model_path,
            nb_threshold=nb_threshold,
            needs_coulomb=self._TRISTATE_MAP[coulomb],
            needs_dispersion=self._TRISTATE_MAP[dispersion],
            device=device_arg,
            compile_model=compile_model,
            ensemble_member=ensemble_member,
        )

        # Post-ctor configuration. If any of these raise, roll back via
        # release() so device-resident tensors that the ctor moved to CUDA
        # are walked back to CPU before the reference is dropped (otherwise
        # a repeat-bad-config server loop accumulates VRAM).
        try:
            # Skip cutoff for "simple" mode (no upstream cutoff parameter).
            if coulomb_method is not None:
                if coulomb_method == "simple":
                    self._calc.set_lrcoulomb_method(coulomb_method)
                else:
                    self._calc.set_lrcoulomb_method(coulomb_method, cutoff=coulomb_cutoff)

            if dftd3_cutoff is not None or dftd3_smoothing_fraction is not None:
                self._calc.set_dftd3_cutoff(dftd3_cutoff, dftd3_smoothing_fraction)

            # Process-wide thread setting; one-shot here, never per-call.
            torch.set_num_threads(ncores)
        except Exception:
            self.release()
            raise

        self._setup_args = new_args

    def release(self) -> None:
        """Release device-side resources held by the cached calculator.

        Server-mode callers (sibling PR
        2026-04-26-oet-server-vram-eviction-design.md) invoke this before
        evicting a cached calculator from the worker cache to reclaim GPU
        memory. Until that sibling PR lands, this method has no live caller.

        NOTE: torch.cuda.empty_cache() returns memory to PyTorch's caching
        allocator, not to the OS — nvidia-smi reservation does not drop.
        Compiled models (compile_model=True) keep Inductor kernel artifacts
        in process-global state that .to('cpu') + empty_cache() cannot
        reclaim; only torch._dynamo.reset() can, and that would invalidate
        every other compiled module in the process.
        """
        if self._calc is not None:
            for component_name in ("model", "external_coulomb", "external_dftd3"):
                component = getattr(self._calc, component_name, None)
                if component is None:
                    continue
                try:
                    component.to("cpu")
                except Exception as e:
                    # Don't let a partial-cleanup failure cascade into the
                    # server's eviction logic, but log so silent VRAM leaks
                    # are debuggable.
                    print(
                        f"Aimnet2Calc.release(): failed to move {component_name} to cpu: {e!r}",
                        file=sys.stderr,
                    )
            self._calc = None
            self._setup_args = None
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add AIMNet2 v0.2 options to an argument parser.

        The same parser is used by both the standalone `oet_aimnet2` CLI
        and by `oet_server aimnet2` (server.py:278 calls extend_parser to
        build its CLI), so flags propagate to server mode automatically.

        Parameters
        ----------
        parser : ArgumentParser
            Parser that should be extended.
        """
        # --- model selection ---------------------------------------------
        parser.add_argument(
            "-m",
            "--model",
            type=str,
            dest="model",
            default="aimnet2",
            help=(
                "Model name (registry alias e.g. aimnet2, aimnet2-2025, "
                "aimnet2-nse, aimnet2-rxn, aimnet2-pd), canonical key "
                "(e.g. aimnet2-wb97m-d3_0), HuggingFace repo id, or "
                "absolute local path to a .pt file. "
                "Default: aimnet2 (= aimnet2-wb97m-d3_0). "
                "For non-covalent / screening: try aimnet2-2025. "
                "For open-shell or charged systems: use aimnet2-nse "
                "(other families are closed-shell-trained)."
            ),
        )
        parser.add_argument(
            "-p",
            "--model-path",
            metavar="DIR",
            dest="model_dir",
            type=str,
            default=str(DEFAULT_MODEL_PATH),
            help=f"Local cache directory for downloaded model files. Default: {DEFAULT_MODEL_PATH}.",
        )
        parser.add_argument(
            "-d",
            "--device",
            metavar="DEVICE",
            dest="device",
            type=str,
            choices=cls._SUPPORTED_DEVICES,
            default="cpu",
            help="Compute device. 'auto' lets upstream auto-detect (None). Default: cpu.",
        )

        # --- performance -------------------------------------------------
        parser.add_argument(
            "--compile",
            dest="compile",
            action="store_true",
            help=(
                "Enable torch.compile JIT. SERVER MODE ONLY - standalone "
                "oet_aimnet2 is a fresh process per ORCA call and re-pays "
                "JIT cost every step. First-call latency 10-60s. Recompiles "
                "on shape change. Incompatible with Hessian. Do NOT use with "
                "NEB / OptTS / IRC."
            ),
        )
        parser.add_argument(
            "--nb-threshold",
            dest="nb_threshold",
            type=int,
            default=120,
            help="Adaptive neighbor-list batch size. Default: 120.",
        )
        parser.add_argument(
            "--ensemble-member",
            dest="ensemble_member",
            type=int,
            choices=(0, 1, 2, 3),
            default=0,
            help=(
                "Define the member of the method ensemble to use (default 0). "
                "OET runs ONE model per call; for production-quality accuracy, "
                "run each member (0,1,2,3) separately and average outside ORCA."
            ),
        )

        # --- long-range Coulomb ------------------------------------------
        parser.add_argument(
            "--coulomb",
            dest="coulomb",
            type=str,
            choices=tuple(cls._TRISTATE_MAP),
            default="auto",
            help="Force on/off the model's long-range Coulomb module. Default: auto (model decides).",
        )
        parser.add_argument(
            "--coulomb-method",
            dest="coulomb_method",
            type=str,
            choices=cls._COULOMB_METHODS,
            default=None,
            help=(
                "Override the model's long-range Coulomb method. "
                "If unset, the model default is used (no post-ctor call)."
            ),
        )
        parser.add_argument(
            "--coulomb-cutoff",
            dest="coulomb_cutoff",
            type=float,
            default=15.0,
            help=(
                "Cutoff in Angstrom for dsf/ewald methods. Default: 15.0. "
                "Requires --coulomb-method (rejected otherwise). "
                "NOTE: aimnet2-rxn family was trained with cutoff frozen at "
                "4.6 A; passing other values fires an upstream UserWarning."
            ),
        )

        # --- dispersion (DFT-D3) -----------------------------------------
        parser.add_argument(
            "--dispersion",
            dest="dispersion",
            type=str,
            choices=tuple(cls._TRISTATE_MAP),
            default="auto",
            help="Force on/off the model's D3 dispersion module. Default: auto (model decides).",
        )
        parser.add_argument(
            "--dftd3-cutoff",
            dest="dftd3_cutoff",
            type=float,
            default=None,
            help="Override D3 dispersion cutoff in Angstrom. Default: model's value.",
        )
        parser.add_argument(
            "--dftd3-smoothing-fraction",
            dest="dftd3_smoothing_fraction",
            type=float,
            default=None,
            help="Override D3 cutoff smoothing fraction. Default: model's value.",
        )

        # --- Special options ----------------------------------------------
        parser.add_argument(
            "-do",
            "--download-only",
            action="store_true",
            default=False,
            dest="download_only",
            help="Only download the model files without performing an actual calculation. "
            "Respects the model and cache dir defined via command line. ",
        )

    def handle_special_args(self, input_args: list[str] | None = None) -> bool:
        """
        Handle special arguments that don't require an input file.

        Parameters
        ----------
        input_args: list [str]
            The input arguments.

        Returns
        -------
        bool
            Whether special arguments were handled and the script can exit or not.
        """
        # Create an argument parser.
        parser = ArgumentParser(add_help=False)
        self.extend_parser(parser=parser)
        args, _ = parser.parse_known_args(input_args)

        # Get the calculator specific parameters
        args_parsed = args
        download_only = args_parsed.download_only
        model = args_parsed.model
        model_dir = args_parsed.model_dir

        # Handle download only.
        if download_only:
            print(f"Downloading model files if necessary to {model_dir}.")
            self.get_model_file(model=model, model_dir=str(model_dir))
            print("Done")
            return True

        return False

    def atomic_symbol_to_number(self, symbol: str) -> int:
        """Convert an element symbol to atomic number.

        Per-model element rejection happens upstream in
        AIMNet2Calculator.eval(validate_species=True). For models that
        do NOT populate metadata["implemented_species"] (e.g. legacy raw
        nn.Module .pt files), upstream silently accepts ANY atomic number
        and may produce undefined output for unsupported elements; OET
        relies on the upstream check and does not second-guess.
        """
        try:
            return _SYMBOL_TO_Z[symbol.title()]
        except KeyError:
            raise ValueError(f"Unknown element symbol: {symbol}")

    def serialize_input(
        self,
        atom_types: list[str],
        coordinates: list[tuple[float, float, float]],
        charge: int,
        mult: int,
        dograd: bool,
    ) -> dict[str, Any]:
        """Build kwargs for AIMNet2Calculator.eval().

        `mult` is only included for NSE-class models; non-NSE models reject
        the key in v0.2.

        Parameters
        ----------
        atom_types : list[str]
            List of element symbols (e.g. ["O", "H", "H"]).
        coordinates : list[tuple[float, float, float]]
            List of (x, y, z) coordinates in Ångström.
        charge : int
            Molecular charge.
        mult : int
            Spin multiplicity (2S + 1).
        dograd : bool
            Whether to compute the gradient.

        Returns
        -------
        dict[str, Any]
            kwargs ready to splat into ``AIMNet2Calculator.eval(**...)``.
            Keys:

            - ``data`` : dict with ``coord`` (list of [N, 3] coordinates),
              ``numbers`` (list of [N] atomic numbers), ``charge``
              (list of [B] total charges), and — only when the loaded
              model is NSE — ``mult`` (list of [B] multiplicities).
            - ``forces`` : bool, mirrors `dograd`.
            - ``stress`` : bool, always False (unsupported by the
              ORCA external-tool protocol).
            - ``hessian`` : bool, always False (Hessian sidecar is
              deferred per cplett's RFC reply).
        """
        numbers = [self.atomic_symbol_to_number(sym) for sym in atom_types]
        data: dict[str, Any] = {
            "coord": [coordinates],
            "numbers": [numbers],
            "charge": [charge],
        }
        # run_aimnet2 guarantees _calc is set before this is reached;
        # the previous defensive `_calc is not None` was unreachable.
        assert self._calc is not None  # for mypy; runtime invariant per caller
        if self._calc.is_nse:
            data["mult"] = [mult]
        elif mult != 1:
            # Non-NSE models are closed-shell-trained; passing mult != 1
            # silently yields a spin-restricted energy, not a true UKS-
            # equivalent open-shell value. Warn loudly (once per process
            # per call site, via the default warnings filter) and direct
            # the user to aimnet2-nse.
            warnings.warn(
                f"mult={mult} requested but the loaded model is closed-shell-"
                "trained (non-NSE). The returned energy is spin-restricted, NOT "
                "a UKS-equivalent open-shell value. For genuine open-shell "
                "energetics use -m aimnet2-nse.",
                UserWarning,
                stacklevel=2,
            )
        return {
            "data": data,
            "forces": dograd,
            "stress": False,
            "hessian": False,
        }

    def run_aimnet2(
        self,
        atom_types: list[str],
        coordinates: list[tuple[float, float, float]],
        calc_data: CalculationData,
    ) -> tuple[float, list[float]]:
        """
        Runs an AimNet2 calculation.

        Parameters
        ----------
        atom_types : list[str]
            List of element symbols (e.g., ["O", "H", "H"])
        coordinates : list[tuple[float, float, float]]
            List of (x, y, z) coordinates
        calc_data: CalculationData
            Object with calculation data for the run

        Returns
        -------
        float
            The computed energy (Eh)
        list[float]
            Flattened gradient vector (Eh/Bohr), if computed, otherwise empty.
        """

        # make ase atoms object for calculation
        aimnet2_input = self.serialize_input(
            atom_types=atom_types,
            coordinates=coordinates,
            mult=calc_data.mult,
            charge=calc_data.charge,
            dograd=calc_data.dograd,
        )

        if not self._calc:
            raise RuntimeError("Calculator could not be initialized.")
        results = self._calc(**aimnet2_input)

        energy = float(results["energy"].detach()) / ENERGY_CONVERSION["eV"]
        gradient = []
        if (forces := results.get("forces", None)) is not None:
            # detach() defends against any future create_graph=True default.
            forces = forces.detach()
            # unit conversion & factor of -1 to convert from forces to gradient
            fac = -LENGTH_CONVERSION["Ang"] / ENERGY_CONVERSION["eV"]
            gradient = (forces * fac).flatten().tolist()

        return energy, gradient

    def calc(
        self,
        calc_data: CalculationData,
        args_parsed: Namespace,
        args_not_parsed: list[str],
    ) -> tuple[float, list[float]]:
        """Routine for calculating energy + optional gradient.

        Validates cross-flag constraints, then sets up the calculator and
        runs run_aimnet2.
        """
        model = args_parsed.model
        model_dir = args_parsed.model_dir
        coulomb_cutoff = args_parsed.coulomb_cutoff
        device = str(args_parsed.device)
        coulomb_method = args_parsed.coulomb_method
        compile_model = args_parsed.compile
        nb_threshold = args_parsed.nb_threshold
        ensemble_member = args_parsed.ensemble_member
        coulomb = args_parsed.coulomb
        dispersion = args_parsed.dispersion
        dftd3_cutoff = args_parsed.dftd3_cutoff
        dftd3_smoothing_fraction = args_parsed.dftd3_smoothing_fraction

        # --- cross-flag validation ----------------------------------------
        # --coulomb-cutoff is only meaningful with --coulomb-method.
        # The 15.0 != check has one corner: a user who passes --coulomb-cutoff 15.0
        # explicitly without --coulomb-method silently gets a default-cutoff
        # Coulomb-disabled run rather than a rejection. They should specify
        # --coulomb-method anyway, so the false-negative is acceptable.
        #
        # We raise SystemExit(msg) directly rather than holding a parser ref:
        # argparse's parser.error() is itself just print + SystemExit(2), and
        # stashing the parser in args_parsed used to poison server.py's cache
        # key (frozenset of args_parsed.items()) because ArgumentParser is
        # hashed by identity. SystemExit's message goes to stderr automatically.
        if coulomb_cutoff != 15.0 and coulomb_method is None:
            raise SystemExit("oet_aimnet2: error: --coulomb-cutoff requires --coulomb-method")

        # Mirror the argparse default so mypy sees a guaranteed-str (the
        # parser sets default=str(DEFAULT_MODEL_PATH)).
        device = str(args_parsed.get("device"))
        if device not in self._SUPPORTED_DEVICES:
            raise RuntimeError(
                f"Device {device} not supported. Use one of {self._SUPPORTED_DEVICES}."
            )

        # --- set up calculator (idempotent) ------------------------------
        self.setup(
            model=model,
            model_dir=model_dir,
            device=device,
            ncores=calc_data.ncores,
            compile_model=compile_model,
            nb_threshold=nb_threshold,
            ensemble_member=ensemble_member,
            coulomb=coulomb,
            dispersion=dispersion,
            coulomb_method=coulomb_method,
            coulomb_cutoff=coulomb_cutoff,
            dftd3_cutoff=dftd3_cutoff,
            dftd3_smoothing_fraction=dftd3_smoothing_fraction,
        )

        # --- read XYZ and run --------------------------------------------
        atom_types, coordinates = xyzfile_to_at_coord(calc_data.xyzfile)
        return self.run_aimnet2(
            atom_types=atom_types,
            coordinates=coordinates,
            calc_data=calc_data,
        )


def main(argv: list[str] | None = None) -> None:
    """Main routine for execution."""
    calculator = Aimnet2Calc()
    if calculator.handle_special_args(argv):
        return
    inputfile, args, args_not_parsed = calculator.parse_args(argv)
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
