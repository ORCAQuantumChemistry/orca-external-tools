#!/usr/bin/env python3
"""
Calculator for using the predict.x binary from
aenet (http://ann.atomistic.net), compatible with ORCA's ExtTool interface.

Provides
--------
class: AenetCalc(BaseCalc)
    Class for performing a Aenet calculation together with ORCA
main: function
    Main function
"""

from argparse import ArgumentParser, Namespace
from collections.abc import Mapping
from pathlib import Path

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import (
    ENERGY_CONVERSION,
    LENGTH_CONVERSION,
    check_path,
    get_nns,
    run_command,
    xyz2xsf,
)


class AenetCalc(BaseCalc):
    # directory, containing the NN files - can be hard-coded here for convenience
    NNPATH: str | Path | None = None  # Path('/path/to/nns')
    # extension for the NN files (<Symbol>.<NNEXT>)
    NNEXT: str | None = None

    @property
    def PROGRAM_NAMES(self) -> list[str]:
        """Program names to search for in PATH"""
        return ["predict.x"]

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add aenet parsing options."""
        parser.add_argument("-x", "--exe", dest="prog", help="Path to the aenet executable")
        parser.add_argument(
            "-n",
            "--nnpath",
            metavar="DIR",
            dest="nnpath",
            required=(not cls.NNPATH),
            help="directory containing the NN files <Symbol>.<EXT>"
            + (f' (default: "{cls.NNPATH}")' if cls.NNPATH else ""),
            default=cls.NNPATH,
        )
        parser.add_argument(
            "-e",
            "--nnext",
            metavar="EXT",
            dest="nnext",
            required=False,
            default=cls.NNEXT,
            help="extension of the NN files. "
            + (
                f"(default: {cls.NNEXT})"
                if cls.NNEXT
                else "If not provided, there must be a single file that matches the glob <Symbol>.*"
            ),
        )

    @staticmethod
    def write_predict_input(
        xsfname: str | Path,
        inpname: str | Path,
        dograd: bool,
        nns: Mapping[str, Path],
    ) -> None:
        """Write the input file for predict.x

        Parameters
        ----------
        xsfname : str | Path
            The name of the XSF coordinates file
        inpname : str | Path
            The file to be written
        dograd : bool
            Whether to compute forces
        nns : Mapping[str, Path]
            Keys are element symbols and values are paths to the NN potential files
        """
        inpname = Path(inpname)
        xsfname = Path(xsfname)
        with inpname.open("w") as f:
            # write types
            f.write(f"TYPES\n{len(nns)}\n" + "\n".join(nns) + "\n\n")
            # write NN paths
            f.write("NETWORKS\n" + "\n".join(f'{a}  "{p}"' for a, p in nns.items()) + "\n\n")
            # forces tag if gradient is requested
            if dograd:
                f.write("FORCES\n\n")
            # write XSF file (only one supported)
            f.write(f"FILES\n1\n{xsfname}\n")

    def run_predict(self, predictexe: str | Path, inpname: str, ncores: int, prog_out: str) -> None:
        """
        Run the predict.x program and redirect its STDOUT and STDERR to a file. Exists on a non-zero return code.

        Parameters
        ----------
        predictexe : str | Path
            Path to the predict.x program
        inpname : str
            Path to the input file
        ncores : int
            Number of cores to use # TODO: currently only implemented in serial
        prog_out: str
            Output file of program
        """
        run_command(predictexe, prog_out, [inpname])

    def read_predict_output(
        self, natoms: int, dograd: bool, prog_out: str
    ) -> tuple[float, list[float]]:
        """Read the output from predict.x

        Parameters
        ----------
        natoms : int
            The number of atoms in the system
        dograd : bool
            Whether the gradient was computed
        prog_out : str
            The output file from predict.x

        Returns
        -------
        float
            The computed energy
        list[float]
            The gradient (X,Y,Z) for each atom
        """
        energy = None
        gradient = []
        with open(prog_out) as f:
            for line in f:
                if "Total energy" in line:
                    fields = line.split()
                    unit = fields[-1]
                    if unit not in ENERGY_CONVERSION:
                        raise ValueError(f"Unknown energy unit: {unit}")
                    energy = float(fields[-2]) / ENERGY_CONVERSION[unit]
                elif dograd and "atomic forces" in line:
                    f.readline()  # empty
                    f.readline()  # x,y,z,Fx,Fy,Fz
                    eunit, lunit = f.readline().split()[-1].strip("()").split("/")
                    if eunit not in ENERGY_CONVERSION:
                        raise ValueError(f"Unknown energy unit: {eunit}")
                    if lunit not in LENGTH_CONVERSION:
                        raise ValueError(f"Unknown length unit: {lunit}")
                    # unit conversion & factor of -1 to convert from forces to gradient
                    fac = -LENGTH_CONVERSION[lunit] / ENERGY_CONVERSION[eunit]
                    f.readline()  # ---
                    for i, line2 in enumerate(f):
                        if i >= natoms:
                            break
                        fields = line2.split()
                        gradient += [float(i) * fac for i in fields[-3:]]
        if not energy:
            raise ValueError(f"Total energy not found in file {prog_out}")
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
        # Get the arguments parsed as defined in extend_parser
        prog = args_parsed.prog
        nnpath_name = args_parsed.nnpath
        nnext = args_parsed.nnext

        # Set and check program path
        calc_data.set_program_path(prog)
        if calc_data.prog_path:
            print(f"Using executable {calc_data.prog_path}")
        else:
            raise FileNotFoundError(
                f"Could not find a valid executable from standard program names: {self.PROGRAM_NAMES}"
            )
        # Check other files
        if not isinstance(nnpath_name, str):
            raise TypeError("Problem detecting nnfiles. Please check your input")
        nnpath = check_path(nnpath_name)

        # set filenames
        namespace = calc_data.basename + ".predict"
        xsfname = namespace + ".xsf"
        inpname = namespace + ".in"
        prog_out = namespace + ".out"

        # process the XYZ file
        natoms, atomtypes = xyz2xsf(xyzname=calc_data.xyzfile, xsfname=xsfname)
        # find the NN files
        nns = get_nns(atomtypes=atomtypes, nnpath=nnpath, nnext=nnext)
        # write the input for predict.x
        self.write_predict_input(xsfname=xsfname, inpname=inpname, dograd=calc_data.dograd, nns=nns)
        # run predict.x
        self.run_predict(
            predictexe=calc_data.prog_path,
            inpname=inpname,
            ncores=calc_data.ncores,
            prog_out=prog_out,
        )
        # parse the output
        energy, gradient = self.read_predict_output(natoms, calc_data.dograd, prog_out)

        return energy, gradient


def main() -> None:
    """Main entry point for wrapper"""
    calculator = AenetCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
