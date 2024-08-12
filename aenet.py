#!/usr/bin/env python3
"""
This is a wrapper for the predict.x binary from aenet (http://ann.atomistic.net),
compatible with ORCA's ExtTool interface.
"""
from __future__ import annotations

import sys
import subprocess
from pathlib import Path
from argparse import ArgumentParser
from typing import Iterable, Mapping


# Set the default parameters here (can be overridden from the command line)
# predict.x executable from aenet
PREDICT_EXE: str | Path | None = "predict.x"  # Path('/path/to/aenet/bin/predict.x-2.0.4-gfortran_serial')
# directory, containing the NN files
NNPATH: str | Path | None = None  # Path('/path/to/nns')
# extension for the NN files (<Symbol>.<NNEXT>)
NNEXT: str | None = None  # '20tanh-20tanh.nn_TM'

# Energy and length conversion to atomic units.
ENERGY_CONVERSION = {'eV': 27.21138625}
LENGTH_CONVERSION = {'Ang': 0.529177210903}


# ----------------------------------------------------------------------------------------------------------------------
# Common functions: these are duplicated in all scripts to make them self-contained

def strip_comments(s: str) -> str:
    """Strip comment starting with '#' and continuing until the end of the string. Also strip whitespace."""
    return s.split('#')[0].strip()


def read_input(inpfile: str | Path) -> tuple[str, int, int, int, bool]:
    """ Read the ORCA-generated input file

    Parameters
    ----------
    inpfile
        The input file

    Returns
    -------
    tuple[str, int, int, int, bool]
        xyzname: str
            Name of the XYZ coordinates file
        charge: int
            Total charge
        mult: int
            Spin multiplicity
        ncores: int
            Number of parallel cores available
        dograd: bool
            Whether to compute the gradient
    """
    with open(inpfile) as f:
        xyzname = strip_comments(f.readline())
        charge = int(strip_comments(f.readline()))
        mult = int(strip_comments(f.readline()))
        ncores = int(strip_comments(f.readline()))
        dograd = bool(int(strip_comments(f.readline())))
        # TODO POINT CHARGES
    return xyzname, charge, mult, ncores, dograd


def write_engrad(outfile: str | Path, natoms: int, energy: float,
                 dograd: bool, gradient: Iterable[float] = None) -> None:
    """ Write the energy/gradient file to feed back to ORCA.

    Parameters
    ----------
    outfile
        The engrad file
    natoms
        Number of atoms
    energy
        Total energy
    dograd
        Whether the gradient is computed
    gradient
        The gradient (X,Y,Z) for each atom
    """
    with open(outfile, 'w') as f:
        output = '#\n'
        output += '# Number of atoms\n'
        output += '#\n'
        output += f'{natoms}\n'
        output += '#\n'
        output += '# Total energy [Eh]\n'
        output += '#\n'
        output += f'{energy:.12e}\n'
        if dograd:
            output += '#\n'
            output += '# Gradient [Eh/Bohr] A1X, A1Y, A1Z, A2X, ...\n'
            output += '#\n'
            output += '\n'.join(f'{g: .12e}' for g in gradient) + '\n'
        f.write(output)


def run_command(command: str | Path, outname: str | Path, *args: str) -> None:
    """
    Run the given command and redirect its STDOUT and STDERR to a file. Exists on a non-zero return code.

    Parameters
    ----------
    command
        The command to run or path to an executable
    outname
        The output file to be written to (overwritten!)
    args
        arguments to be passed to the command
    """
    with open(outname, 'w') as of:
        try:
            subprocess.run([command] + list(args), stdout=of, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as err:
            print(err)
            exit(err.returncode)


def clean_output(outfile: str | Path, namespace: str) -> None:
    """
    Print the output file to STDOUT and remove all files starting with `namespace`

    Parameters
    ----------
    outfile
        The output file to print
    namespace
        The starting string of all files to remove.
    """
    # print the output to STDOUT
    with open(outfile) as f:
        for line in f:  # line by line to avoid memory overflow
            print(line, end='')
    # remove all file from the namespace
    for f in Path('.').glob(namespace + '*'):
        f.unlink()

# ----------------------------------------------------------------------------------------------------------------------


def xyz2xsf(xyzname: str | Path, xsfname: str | Path) -> tuple[int, set[str]]:
    """ Convert a XYZ file to XSF format.

    Parameters
    ----------
    xyzname
        The XYZ file to convert
    xsfname
        The output XSF file name

    Returns
    -------
    tuple[int, set[str]]
        natoms: int
            The number of atoms in the XYZ file
        atomtypes: set[str]
            The elements present in the XYZ file
    """
    atomtypes = set()
    with open(xyzname) as xyzf, open(xsfname, 'w') as xsff:
        natoms = int(xyzf.readline())
        xyzf.readline()  # comment line

        xsff.write("#\n\nATOMS\n")
        for i, line in enumerate(xyzf):
            if i >= natoms:
                break
            # add the forces and print
            xsff.write(line.rstrip() + "  0.0  0.0  0.0\n")
            # collect the elements
            atomtypes.add(line.split()[0])
    return natoms, atomtypes


def get_nns(atomtypes: Iterable[str], nnpath: str | Path, nnext: str = None) -> dict[str, Path]:
    """ Find the neural network potential files for each element in `atomtypes`.
    The files must all be in the same directory and be named "<ElementSymbol>.<Extension>" with the same extension.

    Parameters
    ----------
    atomtypes
        The elements needed
    nnpath
        Path to the directory containing the neural network potential files
    nnext
        The extension for each NN file. If none is given '*' is used as a wildcard.
        However, then there must be a single file that matches, otherwise an exception is raised

    Returns
    -------
    dict[str, Path]
        The keys are element symbols and the values are paths to the NN files

    Raises
    ------
    RuntimeError
        If more than one or no NN files are found for a requested element
    """
    nnpath = Path(nnpath)
    if not nnext:
        nnext = '*'
    nns = {}
    for a in atomtypes:
        matches = list(nnpath.glob(a + '.' + nnext))
        if not matches:
            raise RuntimeError(f'No NN files found for {a} in {nnpath}')
        if len(matches) > 1:
            raise RuntimeError(f'Multiple NN files found for {a}: {matches}. Set --nnext to specify the extension')
        nns[a] = matches[0]
    return nns


def write_predict_input(xsfname: str | Path, inpname: str | Path, dograd: bool, nns: Mapping[str, Path]) -> None:
    """Write the input file for predict.x

    Parameters
    ----------
    xsfname
        The name of the XSF coordinates file
    inpname
        The file to be written
    dograd
        Whether to compute forces
    nns
        Keys are element symbols and values are paths to the NN potential files
    """
    with open(inpname, 'w') as f:
        # write types
        f.write(f'TYPES\n{len(nns)}\n' + '\n'.join(nns) + '\n\n')
        # write NN paths
        f.write('NETWORKS\n' + '\n'.join(f'{a}  "{p}"' for a, p in nns.items()) + '\n\n')
        # forces tag if gradient is requested
        if dograd:
            f.write('FORCES\n\n')
        # write XSF file (only one supported)
        f.write(f'FILES\n1\n{xsfname}\n')


def run_predict(predictexe: str | Path, inpname: str | Path, outname: str | Path, ncores: int) -> None:
    """
    Run the predict.x program and redirect its STDOUT and STDERR to a file. Exists on a non-zero return code.

    Parameters
    ----------
    predictexe
        Path to the predict.x program
    inpname
        Path to the input file
    outname
        The output file to be written to (overwritten!)
    ncores
        Number of cores to use # TODO: currently only implemented in serial
    """
    run_command(predictexe, outname, inpname)


def read_predict_output(outname: str | Path, natoms: int, dograd: bool) -> tuple[float, list[float]]:
    """Read the output from predict.x

    Parameters
    ----------
    outname
        The name of the output file
    natoms
        The number of atoms in the system
    dograd
        Whether the gradient was computed

    Returns
    -------
    tuple[float, list[float]]
        energy: float
            The computed energy
        gradient: list[float]
            The gradient (X,Y,Z) for each atom
    """
    energy = None
    gradient = []
    with open(outname) as f:
        for line in f:
            if 'Total energy' in line:
                fields = line.split()
                unit = fields[-1]
                if unit not in ENERGY_CONVERSION:
                    raise ValueError(f'Unknown energy unit: {unit}')
                energy = float(fields[-2]) / ENERGY_CONVERSION[unit]
            elif dograd and 'atomic forces' in line:
                f.readline()  # empty
                f.readline()  # x,y,z,Fx,Fy,Fz
                eunit, lunit = f.readline().split()[-1].strip('()').split('/')
                if eunit not in ENERGY_CONVERSION:
                    raise ValueError(f'Unknown energy unit: {eunit}')
                if lunit not in LENGTH_CONVERSION:
                    raise ValueError(f'Unknown length unit: {lunit}')
                # unit conversion & factor of -1 to convert from forces to gradient
                fac = - LENGTH_CONVERSION[lunit] / ENERGY_CONVERSION[eunit]
                f.readline()  # ---
                for i, line2 in enumerate(f):
                    if i >= natoms:
                        break
                    fields = line2.split()
                    gradient += [float(i) * fac for i in fields[-3:]]
    return energy, gradient


def main(argv: list[str]):
    """Main function to run the script"""
    # parse the CLI arguments
    parser = ArgumentParser(
        prog=argv[0],
        description="Wrapper for aenet's predict.x, compatible with ORCA's "
        "otool_external. Reads the ORCA-generated input <inputfile>, "
        "converts the BaseName.xyz file to XSF, writes an input for predict.x, "
        "runs the latter, parses its output and writes the BaseName.engrad "
        "file for ORCA.",
    )
    parser.add_argument("inputfile")
    parser.add_argument(
        "-n",
        "--nnpath",
        metavar="DIR",
        dest="nnpath",
        required=(not NNPATH),
        help="directory containing the NN files <Symbol>.<EXT>"
        + (f' (default: "{NNPATH}")' if NNPATH else ""),
        default=NNPATH,
    )
    parser.add_argument(
        "-e",
        "--nnext",
        metavar="EXT",
        dest="nnext",
        required=False,
        default=NNEXT,
        help="extension of the NN files. "
        + (
            f"(default: {NNEXT})"
            if NNEXT
            else "If not provided, there must be a single file that matches the glob <Symbol>.*"
        ),
    )
    parser.add_argument(
        "-x",
        "--exe",
        metavar="EXE",
        dest="exe",
        required=(not PREDICT_EXE),
        help="path to the aenet predict.x executable"
        + (f' (default: "{PREDICT_EXE}")' if PREDICT_EXE else ""),
        default=PREDICT_EXE,
    )
    args = parser.parse_args(argv[1:])

    # sanitize paths
    nnpath = Path(args.nnpath).expanduser().resolve()
    predictexe = Path(args.exe).expanduser().resolve()

    # read the ORCA-generated input
    xyzname, charge, mult, ncores, dograd = read_input(args.inputfile)

    # set filenames
    basename = xyzname.rstrip(".xyz")
    orca_engrad = basename + ".engrad"
    namespace = basename + ".predict"
    xsfname = namespace + ".xsf"
    inpname = namespace + ".in"
    outname = namespace + ".out"

    # process the XYZ file
    natoms, atomtypes = xyz2xsf(xyzname, xsfname)
    # find the NN files
    nns = get_nns(atomtypes, nnpath, args.nnext)
    # write the input for predict.x
    write_predict_input(xsfname, inpname, dograd, nns)
    # run predict.x
    run_predict(predictexe, inpname, outname, ncores)
    # parse the output
    energy, gradient = read_predict_output(outname, natoms, dograd)
    # convert to ORCA engrad
    write_engrad(orca_engrad, natoms, energy, dograd, gradient)
    # pipe the output to STDOUT and remove the generated files
    clean_output(outname, namespace)


if __name__ == '__main__':
    main(sys.argv)
