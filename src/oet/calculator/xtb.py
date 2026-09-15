#!/usr/bin/env python3
"""
This is a simple wrapper for the xtb binary (github.com/grimme-lab/xtb), compatible with ORCA's ExtTool interface.
It is mostly used for testing purposes, since ORCA has a native interface to xtb.

Provides
--------
class: XtbCalc(BaseCalc)
    Class for performing a xtb calculation together with ORCA
main: function
    Main function
"""

import sys
from argparse import ArgumentParser, Namespace

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import (
    check_path,
    mult_to_nue,
    nat_from_xyzfile,
    run_command,
)


class XtbCalc(BaseCalc):
    @property
    def PROGRAM_NAMES(self) -> list[str]:
        """Program names to search for in PATH"""
        return ["xtb", "otool_xtb"]

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add xtb parsing options.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        parser.add_argument("-e", "--exe", dest="prog", help="Path to the xtb executable")

    def read_xtbout(self, calc_data: CalculationData, natoms: int) -> tuple[float, list[float]]:
        """
        Read the output from XTB

        Parameters
        ----------
        calc_data: CalculationData
            Calculation data
        natoms
            number of atoms in the system

        Returns
        -------
        float
            The computed energy
        list[float]
            The gradient (X,Y,Z) for each atom
        """
        xtbgrad = f"{calc_data.basename}.gradient"
        energy = None
        gradient = []
        # read the energy from the output file
        xtbout = check_path(calc_data.output_file)
        with xtbout.open() as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    energy = float(line.split()[3])
                    break
        # read the gradient from the .gradient file
        if calc_data.dograd:
            xtbgrad_path = check_path(xtbgrad)
            natoms_read = 0
            with xtbgrad_path.open() as f:
                for line in f:
                    if "$grad" in line:
                        break
                for line in f:
                    fields = line.split()
                    if len(fields) == 4:
                        natoms_read += 1
                    elif len(fields) == 3:
                        gradient += [float(i) for i in fields]
                    elif "$end" in line:
                        break
                if natoms_read != natoms:
                    print(
                        f"Number of atoms read: {natoms_read} does not match the expected: {natoms}"
                    )
                    sys.exit(1)
                if len(gradient) != 3 * natoms:
                    print(
                        f"Number of gradient entries: {len(gradient)} does not match 3x number of atoms: {natoms}"
                    )
                    sys.exit(1)
        if not energy:
            raise ValueError(f"Total energy not found in file {calc_data.output_file}")
        return energy, gradient

    def run_xtb(
        self,
        calc_data: CalculationData,
        args: list[str],
    ) -> None:
        """
        Run the xtb program with the given input file and redirect its STDOUT and STDERR to a logfile.

        Parameters
        ----------
        args: list[str]
            Arguments not parsed so far
        calc_data: CalculationData
            Object with calculation data for the run
        """
        # Set up a list of command-line arguments for the xtb program.
        # For the xtb package, not all additional command-line arguments (e.g., specifying a solvation model)
        # can be passed before the structure input, so pass the additional args at the end.
        args = [
            str(i)
            for i in [
                calc_data.xyzfile,
                "-c",
                calc_data.charge,
                "-P",
                calc_data.ncores,
                "--namespace",
                calc_data.basename,
            ]
        ] + args
        nue = mult_to_nue(calc_data.mult)
        if nue:
            args += ["-u", str(nue)]
        if calc_data.dograd:
            args += ["--grad"]
        if not calc_data.prog_path:
            raise RuntimeError("Path to program is None.")
        run_command(calc_data.prog_path, calc_data.output_file, args)

    def calc(
        self, calc_data: CalculationData, args_parsed: Namespace, args_not_parsed: list[str]
    ) -> tuple[float, list[float]]:
        """
        Routine for calculating energy and optional gradient.
        Writes ORCA output

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
        # Get parsed options
        prog = args_parsed.prog
        # Set and check the program path if its executable
        calc_data.set_program_path(prog)
        if calc_data.prog_path:
            print(f"Using executable {calc_data.prog_path}")
        else:
            raise FileNotFoundError(
                f"Could not find a valid executable from standard program names: {self.PROGRAM_NAMES}"
            )

        # run xtb
        self.run_xtb(
            calc_data=calc_data,
            args=args_not_parsed,
        )

        # get the number of atoms from the xyz file
        natoms = nat_from_xyzfile(xyz_file=calc_data.xyzfile)

        # parse the xtb output
        energy, gradient = self.read_xtbout(calc_data=calc_data, natoms=natoms)

        return energy, gradient


def main() -> None:
    """
    Main routine for execution
    """
    calculator = XtbCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
