#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse
from collections.abc import Sequence
from pathlib import Path
import shutil
import tomllib

# Set some paths and variables
ROOT = Path(__file__).resolve().parent
PYPROJECT = ROOT / "pyproject.toml"

# Backend extras defined in pyproject.toml that don't come with tests.
NON_BACKEND_EXTRAS = {"dev"}


# Get backends that require a separate venv.
def get_backend_extras() -> list[str]:
    """
    Return backend extras declared in pyproject.toml except the `NON_BACKEND_EXTRAS`.
    """
    with PYPROJECT.open("rb") as handle:
        pyproject = tomllib.load(handle)

    optional_dependencies = (
        pyproject
        .get("project", {})
        .get("optional-dependencies", {})
    )

    return sorted(
        name
        for name in optional_dependencies
        if name not in NON_BACKEND_EXTRAS
    )


EXTRAS = get_backend_extras()

# Minimal python interpreter required by the base class
minimal_python_version = (3, 11)
if sys.version_info < minimal_python_version:
    raise RuntimeError(
        f"Python {minimal_python_version[0]}.{minimal_python_version[1]} or newer is required."
    )


def create_venv(venv_dir: Path) -> None:
    """
    Create virtual environment, if not present.

    Parameters
    ----------
    venv_dir: Path
        Path to the virtual environment
    """
    print(f"Creating virtual environment in '{venv_dir}'...")
    subprocess.check_call([sys.executable, "-m", "venv", str(venv_dir)])
    print("Virtual environment created.")


def get_venv_pip(venv_dir: Path) -> Path:
    """
    Get the path to the `pip` binary inside the virtual environment.

    Parameters
    ----------
    venv_dir: Path
        Path to the virtual environment that should be installed to.

    Returns
    -------
    pip_path: Path
        Path to the `pip` binary inside the virtual environment.

    Raises
    ------
    FileNotFoundError
        if the `pip` binary does not exist or is not executable.
    """
    # Windows pip has a different path
    if os.name == "nt":
            pip_path = (
                venv_dir
                / ("Scripts")
                / ("pip.exe")
            )
    else:
        pip_path = (
            venv_dir
            / ("bin")
            / ("pip")
            )

    if not pip_path.exists():
        raise FileNotFoundError(f"pip not found in venv: {pip_path}")
    return pip_path


def pip_install_target(
    venv_dir: Path,
    script_dir: Path,
    extras: Sequence[str],
    editable: bool,
) -> None:
    """
    Install oet to virtual environment

    Parameters
    ----------
    venv_dir: Path
        Path to the virtual environment
    script_dir: Path
        Path to the final scripts
    extras: Sequence[str]
        Additional extras to be installed. Should match optional dependency groups of 
        the pyproject.toml.
    """

    # Create the directory for storing the final scripts
    script_dir.mkdir(parents=True, exist_ok=True)

    # Get the pip version of the created virtual environment.
    pip_path = get_venv_pip(venv_dir)

    print(f"Installing package to {script_dir} using pip in venv...")

    # Prepare the call for installing the dependencies to the virtual environment.
    target = "."
    if extras:
        target += f"[{','.join(extras)}]"

    command = [pip_path, "install"]

    if editable:
        command.append("-e")
    
    command.append(target)

    # Install the dependencies.
    subprocess.check_call(
        command
    )

    print("Installation complete.")


def copy_oet_scripts(venv_dir: Path, dest_dir: Path, extras: Sequence[str]) -> None:
    """
    Copy all scripts starting with 'oet' from venv/bin to the destination directory.
    Scripts which are "extras" are only copied if actually installed.

    Parameters
    ----------
    venv_dir : Path
        Path to the virtual environment root directory.
    dest_dir : Path
        Directory where the scripts should be copied.
    extras : Sequence[str]
        Installed extras
    """

    # Get the directory with the oet scripts.
    # Exact location depends on the operating system.
    bin_dir = venv_dir / ("Scripts" if os.name == "nt" else "bin")
    if not bin_dir.exists():
        raise FileNotFoundError(f"bin directory not found in venv: {bin_dir}")

    # Make sure the final directory exists.
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Copy the scripts and count how many were copied.
    count = 0
    for script in bin_dir.glob("oet*"):
        if script.is_file():
            # Skip not installed extras
            if (module := script.name.removeprefix("oet_")) in EXTRAS and module not in extras:
                continue
            target = dest_dir / script.name
            # Copy with metadata (so that the scripts remain executable).
            shutil.copy2(script, target)
            print(f"Copied {script.name} → {target}")
            count += 1

    if count == 0:
        print("No oet_ scripts found.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Installation for orca-external-tools package."
    )
    parser.add_argument(
        "--venv-dir",
        "-v",
        type=Path,
        default=Path(".venv"),
        help="Path to the virtual environment directory",
    )
    parser.add_argument(
        "--script-dir",
        "-s",
        type=Path,
        default=Path("./bin"),
        help="Custom directory where bin/packages should be installed",
    )
    parser.add_argument(
        "--extra",
        "-e",
        nargs="+",
        choices=EXTRAS,
        default=[],
        help="Install optional dependency groups defined in pyproject.toml",
    )
    parser.add_argument(
        "--dev",
        "-d",
        action="store_true",
        help="Install optional developer tools.",
    )
    parser.add_argument(
        "--editable",
        action="store_true",
        help="Install in editable mode. Recommended for development.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Create venv
    if not args.venv_dir.exists():
        create_venv(args.venv_dir)
    else:
        print(
            f"Virtual environment already exists in '{args.venv_dir}'.\n"
            "Installing oet to this venv."
        )

    # Setup the extras to be installed
    extras = list(args.extra)
    if args.dev:
        extras.append("dev")

    # Install oet
    pip_install_target(args.venv_dir, args.script_dir, extras, args.editable)

    # Copy scripts for easier usability
    copy_oet_scripts(venv_dir=args.venv_dir, dest_dir=args.script_dir, extras=args.extra)


if __name__ == "__main__":
    main()
