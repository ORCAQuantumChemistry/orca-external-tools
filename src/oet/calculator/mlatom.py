#!/usr/bin/env python3
"""
This is a wrapper for MLatom (http://mlatom.com), compatible with ORCA's ExtTool interface.
Implementation by Pavlo O. Dral (http://dr-dral.com) on 2025.05.04,
based on the example https://github.com/ORCAQuantumChemistry/orca-external-tools/blob/71565ce53837d4d6bdedff71a1ff353f8a289b77/xtb.py
and adapted by FACCTs GmbH on 2026.12.09.

Usage instructions

./mlatom4orca.py <basename_EXT.extinp.tmp> [args for MLatom command line interface]

./mlatom4orca.py -x /path/to/mlatom/shell_cmd.py <basename_EXT.extinp.tmp> [args for MLatom command line interface]

Calculations can be performed with any method or ML model supported by MLatom, which can provide energies (and gradients).

Examples
1. Calculations with AIQM2:

./mlatom4orca.py ethanol_EXT.extinp.tmp AIQM2

2. Calculations with your ML model:

./mlatom4orca.py da_EXT.extinp.tmp useMLmodel MLmodelType=ANI MLmodelIn=da_energy_ani.npz

Note that the calculations via this wrapper calling MLatom for single-point calculations might be significantly slower than using MLatom directly in heavy simulations because of the disk I/O overhead.

Provides
--------
class: MlatomCalc(BaseCalc)
    Class for performing a MLAtom calculation together with ORCA
main: function
    Main function
"""

import os
import shutil
import sys
import tempfile
from argparse import ArgumentParser, Namespace
from pathlib import Path

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import LENGTH_CONVERSION, check_path, run_command


class MlatomCalc(BaseCalc):
    @property
    def PROGRAM_NAMES(self) -> list[str]:
        """
        Program names/paths to search for.

        Returns
        -------
        list[str]
            A list of paths or names to be searched for in PATH.
        """
        exe_names = ["mlatom"]

        candidates: list[str] = []

        # Program names to be searched in the PATH.
        candidates.extend(exe_names)

        # Directories in the current python environments to check
        # Include "Scripts" to catch Windows installations
        bin_dirs = [
            Path(sys.prefix) / "bin",
            Path(sys.prefix) / "Scripts",
            Path(sys.executable).parent,
        ]

        # Add possible full Paths of the program
        for bindir in bin_dirs:
            for name in exe_names:
                candidates.append(str(bindir / name))

        return candidates

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add mlatom parsing options.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        parser.add_argument("-e", "--exe", dest="prog", help="Path to the mlatom executable")

    def run_mlatom(
        self,
        calc_data: CalculationData,
        args: list[str],
    ) -> None:
        """
        Run the MLatom program and redirect its STDOUT and STDERR to a file.

        Parameters
        ----------
        calc_data: CalculationData
            Calculation data
        args : list[str, ...]
            additional arguments to pass to MLatom
        """
        args = list(args)
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            mlatomenergy = os.path.join(cwd, f"{calc_data.basename}.energy")
            mlatomgrad = os.path.join(cwd, f"{calc_data.basename}.gradient")
            shutil.rmtree(tmpdirname)
            shutil.copytree(cwd, tmpdirname)
            args += [
                str(i)
                for i in [
                    f"XYZfile={calc_data.xyzfile}",
                    f"charges={calc_data.charge}",
                    f"multiplicities={calc_data.mult}",
                    f"nthreads={calc_data.ncores}",
                    f"YestFile={mlatomenergy}",
                ]
            ]
            if calc_data.dograd:
                args += [f"YgradXYZestFile={mlatomgrad}"]
            if not calc_data.prog_path:
                raise RuntimeError("Path to program is None.")
            run_command(calc_data.prog_path, calc_data.output_file, args)

    def read_mlatomout(self, calc_data: CalculationData) -> tuple[float, list[float]]:
        """
        Read the output from MLatom

        Parameters
        ----------
        calc_data: CalculationData
            Calculation data

        Returns
        -------
        float
            The computed energy
        list[float] | None
            The gradient (X,Y,Z) for each atom
        """
        mlatomenergy = f"{calc_data.basename}.energy"
        mlatomgrad = f"{calc_data.basename}.gradient"
        energy = None
        gradient = []
        mlatomenergy_path = check_path(mlatomenergy)
        # read the energy from the .energy file
        with mlatomenergy_path.open() as f:
            for line in f:
                energy = float(line)
        # read the gradient from the .gradient file
        if calc_data.dograd:
            mlatomgrad_path = check_path(mlatomgrad)
            icount = 0
            with mlatomgrad_path.open() as f:
                for line in f:
                    icount += 1
                    if icount > 2:
                        gradient += [float(i) * LENGTH_CONVERSION["Ang"] for i in line.split()]
        if not energy:
            raise ValueError(f"Total energy not found in file {calc_data.output_file}")
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
            Calculation data
        args_parsed: Namespace
            Arguments parsed as defined in extend_parser
        args_not_parsed: list[str]
            Arguments not parsed so far

        Returns
        -------
        float
            The computed energy (Eh)
        list[float]
            Flattened gradient vector (Eh/Bohr), if computed, otherwise empty.
        """
        # Get options that were parsed
        prog = args_parsed.prog

        calc_data.set_program_path(prog)
        if calc_data.prog_path:
            print(f"Using executable {calc_data.prog_path}")
        else:
            raise FileNotFoundError(
                f"Could not find a valid executable from standard program names: {self.PROGRAM_NAMES}"
            )

        # run MLatom
        self.run_mlatom(
            calc_data=calc_data,
            args=args_not_parsed,
        )

        # parse the MLatom output
        energy, gradient = self.read_mlatomout(calc_data=calc_data)

        return energy, gradient


def main() -> None:
    """
    Main routine for execution
    """
    calculator = MlatomCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
