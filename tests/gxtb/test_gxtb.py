import shutil
import unittest
from pathlib import Path

from oet.core.test_utilities import (
    OH,
    WATER,
    get_filenames,
    read_result_file,
    run_wrapper,
    write_input_file,
    write_xyz_file,
)

# Get the path to the script that should be tested
resolved_gxtb_script = shutil.which("oet_gxtb")
if resolved_gxtb_script is None:
    raise RuntimeError(
        "The 'goet_xtb' script was not found in PATH. "
        "Run the tests with the project's virtual environment activated."
    )
gxtb_script_path = Path(resolved_gxtb_script)

# Leave gxtb_executable_path empty, if gxtb from system path should be called
gxtb_executable_path = ""


def run_gxtb(inputfile: str, output_file: str) -> None:
    if gxtb_executable_path:
        arguments = ["--exe", gxtb_executable_path]
    else:
        arguments = None
    run_wrapper(
        inputfile=inputfile, script_path=gxtb_script_path, outfile=output_file, args=arguments
    )


class GxtbTests(unittest.TestCase):
    def test_H2O_engrad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("H2O")

        write_xyz_file(xyz_file, WATER)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=0,
            multiplicity=1,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb(input_file, output_file)
        expected_num_atoms = 3
        expected_energy = -76.43736490412
        expected_gradients = [
            -8.58374584e-03,
            -6.34732203e-03,
            4.48788670e-03,
            3.68390440e-03,
            5.26218976e-03,
            -4.49684003e-04,
            4.89984144e-03,
            1.08513227e-03,
            -4.03820270e-03,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=9)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=9)

    def test_OH_anion_eng_grad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_anion")
        write_xyz_file(xyz_file, OH)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=-1,
            multiplicity=1,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb(input_file, output_file)
        expected_num_atoms = 2
        expected_energy = -75.80305584316
        expected_gradients = [
            2.28916816e-03,
            7.36155354e-03,
            2.09936121e-03,
            -2.28916816e-03,
            -7.36155354e-03,
            -2.09936121e-03,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=9)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=9)

    def test_OH_rad_eng_grad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_rad")
        write_xyz_file(xyz_file, OH)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=0,
            multiplicity=2,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb(input_file, output_file)
        expected_num_atoms = 2
        expected_energy = -75.74502880794
        expected_gradients = [
            -1.02890363e-04,
            -3.55911885e-04,
            -1.29478984e-04,
            1.02890363e-04,
            3.55911885e-04,
            1.29478984e-04,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=7)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=7)


if __name__ == "__main__":
    unittest.main()
