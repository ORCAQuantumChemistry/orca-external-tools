#!/usr/bin/env python3
"""
This is a simple wrapper for the xtb binary (github.com/grimme-lab/xtb), compatible with ORCA's ExtTool interface.
It is mostly used for testing purposes, since ORCA has a native interface to xtb.
"""
from __future__ import annotations

import shutil
import sys
import subprocess
from pathlib import Path
from argparse import ArgumentParser
from typing import Iterable


# path to the xtb executable. If None, will look for all XTB_NAMES in the system PATH
XTB_EXE: str | Path | None = None
XTB_NAMES: list[str] = ['xtb', 'otool_xtb']


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
            subprocess.run([str(command)] + list(args), stdout=of, stderr=subprocess.STDOUT, check=True)
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

def run_xtb(xtbexe: str | Path, xyzname: str, namespace: str,
            charge: int, mult: int, ncores: int, dograd: bool,
            outfile: str | Path, *args: str) -> None:
    """
    Run the XTB program and redirect its STDOUT and STDERR to a file.

    Parameters
    ----------
    xtbexe
        path to the xtb binary
    xyzname
        name of the XYZ file
    namespace
        filename prefix for the xtb output files
    charge
        total charge of the system
    mult
        spin multiplicity of the system
    ncores
        number of threads to use
    dograd
        whether to compute the gradient
    outfile
        the output file
    args
        additional arguments to pass to xtb
    """
    args = list(args)
    args += [str(i) for i in [xyzname, '-c', charge, '-P', ncores, '--namespace', namespace]]
    nue = mult - 1
    if nue:
        args += ['-u', str(nue)]
    if dograd:
        args += ['--grad']
    run_command(xtbexe, outfile, *args)


def read_xtbout(namespace: str, xtbout: str | Path, natoms: int, dograd: bool) -> tuple[float, list[float]]:
    """
    Read the output from XTB

    Parameters
    ----------
    namespace
        filename prefix of the xtb output files
    xtbout
        the main xtb output file
    natoms
        number of atoms in the system
    dograd
        whether to read the gradient

    Returns
    -------
    tuple[float, list[float]]
        energy: float
            The computed energy
        gradient: list[float]
            The gradient (X,Y,Z) for each atom
    """
    xtbgrad = namespace + '.gradient'
    energy = None
    natoms_read = 0
    gradient = []
    # read the energy from the output file
    with open(xtbout) as f:
        for line in f:
            if 'TOTAL ENERGY' in line:
                energy = float(line.split()[3])
                break
    # read the gradient from the .gradient file
    if dograd:
        with open(xtbgrad) as f:
            for line in f:
                if '$grad' in line:
                    break
            for line in f:
                fields = line.split()
                if len(fields) == 4:
                    natoms_read += 1
                elif len(fields) == 3:
                    gradient += [float(i) for i in fields]
                elif '$end' in line:
                    break
            if natoms_read != natoms:
                print(f'Number of atoms read: {natoms_read} does not match the expected: {natoms}')
                exit(1)
            if len(gradient) != 3 * natoms:
                print(f'Number of gradient entries: {len(gradient)} does not match 3x number of atoms: {natoms}')
                exit(1)
    return energy, gradient


def main(argv):
    if not (xtbexe := XTB_EXE):
        for xtb in XTB_NAMES:
            if xtbexe := shutil.which(xtb):
                break

    # parse the CLI arguments
    parser = ArgumentParser(prog=argv[0], allow_abbrev=False,
                            description="Wrapper for xtb, compatible with ORCA's otool_external. "
                                        "Reads the ORCA-generated input <inputfile>, calls xtb, "
                                        "parses its output and writes the BaseName.engrad file for ORCA.")
    parser.add_argument('inputfile')
    parser.add_argument('-x', "--exe", metavar='XTBEXE', dest='xtbexe', required=(not xtbexe),
                        help="path to the xtb executable" +
                             (f' (default: "{xtbexe}")' if xtbexe else ''), default=xtbexe)
    args, xtb_args = parser.parse_known_args(argv[1:])

    # sanitize the path to xtb
    xtbexe = Path(args.xtbexe).expanduser().resolve()

    # read the ORCA-generated input
    xyzname, charge, mult, ncores, dograd = read_input(args.inputfile)

    # set filenames
    basename = xyzname.rstrip(".xyz")
    orca_engrad = basename + ".engrad"
    xtb_namespace = basename + '.xtb'
    xtbout = xtb_namespace + ".out"

    # run xtb
    run_xtb(xtbexe, xyzname, xtb_namespace, charge, mult, ncores, dograd, xtbout, *xtb_args)

    # get the number of atoms from the xyz file
    with open(xyzname) as f:
        natoms = int(f.readline())

    # parse the xtb output
    energy, gradient = read_xtbout(xtb_namespace, xtbout, natoms, dograd)

    # write the ORCA engrad file
    write_engrad(orca_engrad, natoms, energy, dograd, gradient)

    # print the xtb output to STDOUT and remove leftover files
    clean_output(xtbout, xtb_namespace)


if __name__ == '__main__':
    main(sys.argv)


