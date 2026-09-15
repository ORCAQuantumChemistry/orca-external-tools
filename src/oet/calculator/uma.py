#!/usr/bin/env python3
"""
This package provides [fairchem](https://github.com/facebookresearch/fairchem) wrappers for ORCA's ExtTool interface.
Before starting to use this module, please make sure your have access to the [fairchem repository](https://huggingface.co/facebook/UMA) and logged in with your
huggingface account. For details, please see the GitHub repository or the [respective tutorials](https://fair-chem.github.io/).

Provides
--------
class: UmaCalc(CalcServer)
    Class for performing a UMA calculation together with ORCA
main: function
    Main function
"""

import os
import sys
import warnings
from argparse import ArgumentParser, Namespace
from pathlib import Path

from oet import ASSETS_DIR
from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import ENERGY_CONVERSION, LENGTH_CONVERSION, xyzfile_to_at_coord

try:
    # Suppress pkg_resources deprecated warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from ase.calculators.calculator import PropertyNotImplementedError
        from fairchem.core import FAIRChemCalculator, pretrained_mlip
        from fairchem.core.calculate.pretrained_mlip import available_models
        from fairchem.core.units.mlip_unit.api.inference import UMATask
        from huggingface_hub import hf_hub_download
except ImportError as e:
    print(
        f"[MISSING] Required module fairchem-core not found: {e}.\n"
        "Please install the packages in the virtual environment.\n"
        "Therefore, activate the venv, got to the orca-external-tools "
        "main directory and use pip install -r ./requirements/uma.txt\n"
        "Also, make sure you are logged in with your Hugging Face account:\n"
        "https://fair-chem.github.io/core/install.html#access-to-gated-models-on-huggingface"
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


# Override the default fairchem `CACHE_DIR`, unless the environment variable is set
DEFAULT_CACHE_DIR = str(os.environ.get("FAIRCHEM_CACHE_DIR", ASSETS_DIR / "fairchem"))


class UmaCalc(BaseCalc):
    # Fairchem calculator used to compute energy and grad
    _calc: FAIRChemCalculator | None = None

    def set_calculator(
        self, param: str, basemodel: str, device: str, cache_dir: str, force: bool = False
    ) -> None:
        """
        Prepare the `FAIRChemCalculator` object to compute energy and gradient, if not done already.

        Parameters
        ----------
        param: str
            Parameter set used by fairchem
        basemodel: str
            UMA basemodel
        device: str
            Device that should be used, e.g., cpu or cuda
        cache_dir: str
            Cache directory to read/write downloaded model files to
        force: bool, default = False
            Force re-initialization of the calculator, even if already initialized
        """
        if not self._calc or force:
            # Make sure the cache directory exists
            Path(cache_dir).mkdir(parents=True, exist_ok=True)
            # Monkey-patch the Fairchem CACHE_DIR: the provided one is not always respected.
            # In particular, `pretrained_checkpoint_path_from_name` just uses `CACHE_DIR`
            pretrained_mlip.CACHE_DIR = cache_dir
            # Suppress fairchemcore internal warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                predictor = pretrained_mlip.get_predict_unit(
                    basemodel, device=device, cache_dir=cache_dir
                )
                self._calc = FAIRChemCalculator(predictor, task_name=param)

    def get_calculator(self) -> FAIRChemCalculator:
        """
        Returns UMA calculator
        """
        return self._calc

    def check_for_model_files(self, basemodel: str, cache_dir: str) -> bool:
        """
        Check if model files are available in current cache directory.

        Parameters
        ----------
        basemodel: str
            UMA basemodel
        cache_dir: str
            Cache directory

        Returns
        -------
        bool
            True, if model files were found
        """
        try:
            # First the model parameter
            hf_hub_download(
                filename=basemodel + ".pt",
                repo_id="facebook/UMA",
                subfolder="checkpoints",
                cache_dir=cache_dir,
                local_files_only=True,
            )
            # Then the atomic references
            hf_hub_download(
                filename="iso_atom_elem_refs.yaml",
                repo_id="facebook/UMA",
                subfolder="references",
                cache_dir=cache_dir,
                local_files_only=True,
            )
        except Exception:
            return False
        else:
            return True

    def switch_to_offline_mode(self) -> None:
        """
        Goes into offline mode to prevent HuggingFace from downloading something from the web
        """
        os.environ["HF_HUB_OFFLINE"] = "1"
        # Older versions require a different keyword:
        os.environ["HUGGINGFACE_HUB_OFFLINE"] = "1"

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add Uma parsing options.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        parser.add_argument(
            "-t",
            "--task",
            type=UMATask,
            choices=list(UMATask),
            default=(default_task := UMATask.OMOL),
            metavar="TASK",
            dest="param",
            help="The UMA task/parameter set name. "
            "Options: " + ", ".join(UMATask) + ". "
            f"Default: {default_task}. ",
        )
        parser.add_argument(
            "-m",
            "--model",
            type=str,
            default=(default_model := "uma-s-1p1"),
            metavar="MODEL",
            dest="basemodel",
            choices=available_models,
            help="The UMA base model. "
            "Options: " + ", ".join(available_models) + ". "
            f"Default: {default_model}. ",
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
        parser.add_argument(
            "-c",
            "--cachedir",
            type=str,
            default=str(DEFAULT_CACHE_DIR),
            metavar="DIR",
            dest="cache_dir",
            help="The cache directory to store downloaded model files. "
            "Can also be set via the environment variable FAIRCHEM_CACHE_DIR. "
            f'Default: "{DEFAULT_CACHE_DIR}".',
        )
        parser.add_argument(
            "-do",
            "--download-only",
            action="store_true",
            default=False,
            dest="download_only",
            help="Only download the model files without performing an actual calculation. "
            "Respects the model and cache dir defined via command line. ",
        )
        parser.add_argument(
            "-o",
            "--offline",
            action="store_true",
            default=False,
            dest="offline_mode",
            help="Force into offline mode. Please note that there will be an error if the model parameters are not found.",
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
        # Create an argument parser and parse.
        parser = ArgumentParser(add_help=False)
        self.extend_parser(parser=parser)
        args_parsed, _ = parser.parse_known_args(input_args)

        # Check if the model files should only be downloaded.
        download_only = args_parsed.download_only

        # Handle download only.
        if download_only:
            # Get the calculator specific parameters
            param = args_parsed.param
            basemodel = args_parsed.basemodel
            device = args_parsed.device
            cache_dir = args_parsed.cache_dir
            # Download the files if necessary.
            print(f"Downloading model files if necessary to {DEFAULT_CACHE_DIR}.")
            self.set_calculator(
                param=param, basemodel=basemodel, device=device, cache_dir=cache_dir
            )
            print("Done")
            return True

        return False

    def run_uma(
        self,
        atom_types: list[str],
        coordinates: list[tuple[float, float, float]],
        calc_data: CalculationData,
    ) -> tuple[float, list[float]]:
        """
        Runs an UMA calculation.

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

        # Suppress fairchemcore internal warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        Writes ORCA output


        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        args_parsed: dict[str, Any]
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
        # Get the arguments parsed as defined in extend_parser
        param = args_parsed.param
        basemodel = args_parsed.basemodel
        device = args_parsed.device
        cache_dir = args_parsed.cache_dir
        offline_mode = args_parsed.offline_mode
        # Check if the model files are available
        model_files_available = self.check_for_model_files(basemodel=basemodel, cache_dir=cache_dir)
        # If they are available, switch to offline mode.
        if model_files_available:
            self.switch_to_offline_mode()
            print("Model files detected. Switching to offline mode.")
        # If they are not available, but the user requested the offline mode to prevent online communication,
        # print a warning as this will likely cause subsequent errors.
        elif offline_mode:
            self.switch_to_offline_mode()
            # Check if the model files are locally available. If not, subsequent errors will occur
            # as they cannot be downloaded.
            print(
                "WARNING: Offline mode selected, but no model files were detected. "
                "This will likely cause subsequent errors."
            )

        # setup calculator if not already set
        # this is important as usage on a server would otherwise cause
        # initialization with every call so that nothing is gained
        self.set_calculator(param=param, basemodel=basemodel, device=device, cache_dir=cache_dir)

        # process the XYZ file
        atom_types, coordinates = xyzfile_to_at_coord(calc_data.xyzfile)

        # run uma
        energy, gradient = self.run_uma(
            atom_types=atom_types, coordinates=coordinates, calc_data=calc_data
        )

        return energy, gradient


def main() -> None:
    """
    Main routine for execution
    """
    calculator = UmaCalc()
    if calculator.handle_special_args():
        return
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
