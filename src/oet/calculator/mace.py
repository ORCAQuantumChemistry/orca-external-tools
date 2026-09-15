#!/usr/bin/env python3
"""
MACE wrapper for ORCA's ExtTool interface.
"""

import sys
import warnings
from argparse import ArgumentParser, Namespace

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import ENERGY_CONVERSION, LENGTH_CONVERSION, xyzfile_to_at_coord

try:
    with warnings.catch_warnings():
        # Catch the warning that the `weights_only` bool is not set during `torch.load`. This will not affect the calculation.
        warnings.filterwarnings(
            "ignore",
            message="Environment variable TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected",
        )
        from ase.calculators.calculator import PropertyNotImplementedError
        from mace.calculators.foundations_models import mace_mp, mace_omol
        from mace.calculators.mace import MACECalculator
except ImportError as e:
    print(
        f"[MISSING] Required module mace not found: {e}.\n"
        "Please install the packages in the virtual environment.\n"
        "Therefore, activate the venv, got to the orca-external-tools "
        "main directory and use pip install -r ./requirements/mace.txt"
    )
    sys.exit(1)

try:
    import torch
except ImportError as e:
    print("[MISSING] torch not found:", e)
    sys.exit(1)

try:
    from ase import Atoms
except ImportError as e:
    print("[MISSING] ase not found:", e)
    sys.exit(1)


class MaceCalc(BaseCalc):
    # MACE calculator used for computations
    _calc: MACECalculator | None = None

    def set_calculator(
        self,
        suite: str,
        model: str,
        dispersion: bool,
        damping: str,
        dispersion_xc: str,
        dispersion_cutoff: float,
        device: str,
        default_dtype: str,
        head: str,
        force: bool = False,
    ) -> None:
        """
        Prepare a MACE calculator object that is compatible with the ASE calculator object to compute energy and gradient, if not done already.

        Parameters
        ----------
        suite: str
            Suite for the calculator.
        model: str
            Model/model path to be used.
        dispersion: bool
            Whether to use a dispersion model or not.
        damping: str
            Dispersion damping.
        dispersion_xc: str
            XC functional for the dispersion model.
        dispersion_cutoff: str
            Cutoff radius for the dispersion model.
        device: str
            Device (cuda or cpu).
        default_dtype: str
            Precision type.
        head: str
            The head of the MACE model.
        force: bool, default = False
            Force re-initialization of the calculator, even if already initialized.
        """
        if self._calc and not force:
            return
        match suite:
            case "mp":
                kwargs: dict[str, object] = {
                    "model": model,
                    "device": device,
                    "default_dtype": default_dtype or "float32",
                    "dispersion": dispersion,
                    "damping": damping,
                    "dispersion_xc": dispersion_xc,
                    "dispersion_cutoff": (dispersion_cutoff * LENGTH_CONVERSION["Ang"]),
                }
                if head:
                    kwargs["head"] = head
                calc = mace_mp(**kwargs)
            # No dispersion applicable for omol
            case "omol":
                calc = mace_omol(
                    model=model,
                    device=device,
                    default_dtype=default_dtype or "float64",
                )
            case _:
                raise RuntimeError(f"The suite {suite} is not supported.")
        self._calc = calc

    def get_calculator(self) -> MACECalculator | None:
        """
        Returns the current MACE calculator
        """
        return self._calc

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """
        Add Mace parsing options.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        parser.add_argument(
            "-s",
            "--suite",
            dest="suite",
            choices=["mp", "omol", "mace-mp", "mace-omol"],
            default="omol",
            help="Select MACE suite: mp/mace-mp or omol/mace-omol. Default: omol",
        )
        # Settings the default model to None uses the MACE default.
        parser.add_argument(
            "-m",
            "--model",
            type=str,
            default=None,
            help="Model spec or local path. MP: small/medium/large/medium-mpa-0/... OMOL: extra_large or path.",
        )
        parser.add_argument(
            "--default-dtype",
            choices=["float32", "float64"],
            default="float64",
            help="Default float precision (recommended: use float32 only for MD).",
        )
        parser.add_argument(
            "-d",
            "--device",
            type=str,
            default="cpu",
            metavar="DEVICE",
            dest="device",
            choices=(device_choices := ["cpu", "cuda"]),
            help="Device to perform the calculation on. "
            "Options: " + ", ".join(device_choices) + ". "
            "Default: cpu. ",
        )
        # MP-specific extras
        parser.add_argument(
            "--dispersion",
            action="store_true",
            help=(
                "Enable D3 dispersion. Only applicable for the MP suite. "
                "This requires the installation of torch_dftd to the virtual "
                "environment and that liblzma was installed when setting up "
                "the environment."
            ),
        )
        parser.add_argument(
            "--damping",
            type=str,
            default="bj",
            choices=(damp_choices := ["zero", "bj", "zerom", "bjm"]),
            help="D3 damping. "
            "Options: " + ", ".join(damp_choices) + "."
            "Only applicable for the MP suite. Default BJ for MP suite. ",
        )
        parser.add_argument(
            "--dispersion-xc",
            type=str,
            default="pbe",
            help="XC functional for D3. Only applicable for the MP suite. Default PBE for MP suite.",
        )
        parser.add_argument(
            "--dispersion-cutoff",
            type=float,
            default=40.0,
            help="Cutoff radius for D3 in Bohr. "
            "Only applicable for the MP suite. "
            "Default: 40 Bohr for MP suite.",
        )
        parser.add_argument(
            "--head",
            type=str,
            default=None,
            help="Advanced option: Select MACE head. Only applicable for the MP suite.",
        )

    def run_mace(
        self,
        atom_types: list[str],
        coordinates: list[tuple[float, float, float]],
        calc_data: CalculationData,
    ) -> tuple[float, list[float]]:
        """
        Runs a Mace calculation.

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
            Flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """

        # set the number of threads
        torch.set_num_threads(calc_data.ncores)

        # make ase atoms object for calculation
        atoms = Atoms(symbols=atom_types, positions=coordinates)
        atoms.info = {"charge": calc_data.charge, "spin": calc_data.mult}
        atoms.calc = self._calc

        # Get the energy and gradients
        energy = atoms.get_potential_energy() / ENERGY_CONVERSION["eV"]
        gradient = []
        try:
            forces = atoms.get_forces()
            # Convert forces to gradient (-1) and unit conversion
            fac = -LENGTH_CONVERSION["Ang"] / ENERGY_CONVERSION["eV"]
            gradient = (fac * forces).flatten().tolist()
        except PropertyNotImplementedError:
            # forces may not be available
            pass

        return energy, gradient

    def calc(
        self,
        calc_data: CalculationData,
        args_parsed: Namespace,
        args_not_parsed: list[str],
    ) -> tuple[float, list[float]]:
        """
        Routine for calculating energy and optional gradient.

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        args_parsed: Namespace
            Arguments parsed as defined in extend_parser
        args_not_parsed: list[str]
            Arguments not parsed so far

        Returns
        -------
        float
            The computed energy (Eh)
        list[float]
            Flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """

        suite = args_parsed.suite
        if suite.startswith("mace"):
            suite = suite.split("-", 1)[1]
        dispersion = args_parsed.dispersion
        dispersion_xc = args_parsed.dispersion_xc
        dispersion_cutoff = args_parsed.dispersion_cutoff
        if dispersion and suite != "mp":
            print(
                "WARNING: Dispersion flag recognized, but MP suite not used. Ignoring all options related to dispersion."
            )

        # Setup calculator if not already set.
        # This is important as usage on a server would otherwise cause an initialization with every call so that nothing is gained.
        # Catch the warning that the `weights_only` bool is not set during `torch.load`. This will not affect the calculation.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Environment variable TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected",
            )
            self.set_calculator(
                suite=suite,
                model=args_parsed.model,
                dispersion=dispersion,
                damping=args_parsed.damping,
                dispersion_xc=dispersion_xc,
                dispersion_cutoff=dispersion_cutoff,
                device=args_parsed.device,
                default_dtype=args_parsed.default_dtype,
                head=args_parsed.head,
            )

        # process the XYZ file
        atom_types, coordinates = xyzfile_to_at_coord(calc_data.xyzfile)

        # run mace
        energy, gradient = self.run_mace(
            atom_types=atom_types, coordinates=coordinates, calc_data=calc_data
        )

        return energy, gradient


def main() -> None:
    """
    Main routine for execution.
    """
    calculator = MaceCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
