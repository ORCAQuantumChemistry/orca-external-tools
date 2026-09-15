import shutil
import subprocess
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

# Reference values regenerated against aimnet v0.2 (PyPI). The v0.2
# release ships retrained model files (storage path .../aimnet2v2/...)
# whose energies differ from v0.1.x by ~1e-6 Eh at this geometry.
# v0.2 is bit-exact deterministic across runs and between standalone
# wrapper and server paths, so places=6 holds.

# Model for running the tests
aimnet_model = "aimnet2"

# Default maximum time (in sec) to download the model files if not present
timeout = 600

# Get the path to the script that should be tested
resolved_aimnet2_script = shutil.which("oet_aimnet2")
if resolved_aimnet2_script is None:
    raise RuntimeError(
        "The 'oet_aimnet2' script was not found in PATH. "
        "Run the tests with the project's virtual environment activated."
    )
aimnet2_script_path = Path(resolved_aimnet2_script)


def cache_model_files(model: str) -> None:
    """
    Wrapper to check if the required model files are present. If not, they are downloaded.

    model: str
        Model for computing the test cases.
    """
    subprocess.run(
        [
            aimnet2_script_path,
            "--download-only",
            "--model",
            model,
        ],
        timeout=timeout,
        check=True,
    )


def run_aimnet2(inputfile: str, output_file: str) -> None:
    run_wrapper(
        inputfile=inputfile, script_path=aimnet2_script_path, outfile=output_file, timeout=30
    )


class Aimnet2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Downloading the model files if necessary.
        """
        # Pre-download AIMNet2 model files
        print("Checking the model files and downloading them if necessary.")

        try:
            cache_model_files(aimnet_model)
        except subprocess.TimeoutExpired as e:
            raise TimeoutError(
                "Loading the model files timed out. "
                "Please check your internet connection and consider "
                "increasing the timeout."
            ) from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError("Loading the model files failed.") from e

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
        run_aimnet2(input_file, output_file)
        expected_num_atoms = 3
        expected_energy = -7.647682538153e01
        expected_gradients = [
            -1.020942814648e-02,
            -7.558954879642e-03,
            5.339907482266e-03,
            3.577803261578e-03,
            9.023892693222e-03,
            1.832913840190e-03,
            6.631619296968e-03,
            -1.464935485274e-03,
            -7.172822486609e-03,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=6)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=6)

    def test_OH_anion_eng_grad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_ainion")
        write_xyz_file(xyz_file, OH)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=-1,
            multiplicity=1,
            ncores=2,
            do_gradient=1,
        )
        run_aimnet2(input_file, output_file)
        expected_num_atoms = 2
        expected_energy = -7.582629635076e01
        expected_gradients = [
            -4.858376923949e-04,
            -1.563820987940e-03,
            -4.455552552827e-04,
            4.858376923949e-04,
            1.563823316246e-03,
            4.455552552827e-04,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=6)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=6)

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
        run_aimnet2(input_file, output_file)
        expected_num_atoms = 2
        expected_energy = -7.568258700191e01
        expected_gradients = [
            -3.783945925534e-03,
            -1.217983383685e-02,
            -3.470211755484e-03,
            3.783945692703e-03,
            1.217983569950e-02,
            3.470211755484e-03,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=6)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=6)


if __name__ == "__main__":
    unittest.main()
