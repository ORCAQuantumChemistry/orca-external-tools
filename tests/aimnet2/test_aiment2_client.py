import os
import shutil
import signal
import subprocess
import unittest
from pathlib import Path

from oet.core.test_utilities import (
    OH,
    WATER,
    get_filenames,
    read_result_file,
    run_wrapper,
    wait_for_server,
    write_input_file,
    write_xyz_file,
)

# Get the path to the script that should be tested
resolved_aimnet2_client = shutil.which("oet_client")
if resolved_aimnet2_client is None:
    raise RuntimeError(
        "The 'oet_client' script was not found in PATH. "
        "Run the tests with the project's virtual environment activated."
    )
aimnet2_client_path = Path(resolved_aimnet2_client)

resolved_server_script = shutil.which("oet_server")
if resolved_server_script is None:
    raise RuntimeError(
        "The 'oet_server' script was not found in PATH. "
        "Run the tests with the project's virtual environment activated."
    )
aimnet2_server_path = Path(resolved_server_script)

# Get the path to the script for downloading model files.
resolved_aimnet2_script = shutil.which("oet_aimnet2")
if resolved_aimnet2_script is None:
    raise RuntimeError(
        "The 'oet_aimnet2' script was not found in PATH. "
        "Run the tests with the project's virtual environment activated."
    )
aimnet2_script_path = Path(resolved_aimnet2_script)

# Default ID and port of server. Change if needed
ip_port = "127.0.0.1:9000"

# Model for running the tests
aimnet_model = "aimnet2"

# Default maximum time (in sec) to download the model files if not present
timeout = 600


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
        inputfile=inputfile,
        script_path=aimnet2_client_path,
        outfile=output_file,
        args=["--bind", ip_port],
        timeout=30,
    )


class Aimnet2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Downloading the model files if necessary and starting the server
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

        # Set up the server
        server_out = Path("server.out").resolve()
        print(f"Starting the server. A detailed server log can be found on file {server_out}")
        with open(server_out, "a") as f:
            cls.server = subprocess.Popen(
                [
                    aimnet2_server_path,
                    "aimnet2",
                    "--bind",
                    ip_port,
                    "--nthreads",
                    "2",
                ],
                stdout=f,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        # Wait for the server to be ready.
        wait_for_server(
            process=cls.server,
            ip_port=ip_port,
            timeout=30.0,
        )

    @classmethod
    def tearDownClass(cls):
        """
        Shut the server at the end
        """
        print("Killing the server.")
        os.killpg(os.getpgid(cls.server.pid), signal.SIGTERM)
        cls.server.wait(timeout=10)

    def test_H2O_engrad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("H2O_client")

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
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_anion_client")
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
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_rad_client")
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
