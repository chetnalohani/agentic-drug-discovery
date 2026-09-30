"""
Simple AutoDock Vina Docking Agent
Agentic Drug Discovery Project
"""

from pathlib import Path
import subprocess
import shutil


def find_project_root():
    """Find project root directory."""
    return Path(__file__).resolve().parent.parent


def find_vina(project_root):
    return Path(project_root) / "vina" / "vina.exe"

    vina = project_root / "vina" / "vina.exe"

    if vina.exists():
        return vina

    system_vina = shutil.which("vina")

    if system_vina:
        return Path(system_vina)

    return None


def find_receptor(project_root, receptor):
    """Find a usable receptor PDBQT file."""

    receptor = Path(receptor)

    # If already PDBQT
    if receptor.suffix.lower() == ".pdbqt" and receptor.exists():
        return receptor

    # Check requested path
    if receptor.exists():
        # Look for PDBQT with same name
        pdbqt = receptor.with_suffix(".pdbqt")

        if pdbqt.exists():
            return pdbqt

    # Common prepared receptor locations
    possible = [
        project_root / "1M17.pdbqt",
        project_root / "results" / "docking" / "1M17_receptor.pdbqt",
        project_root / "results" / "docking" / "EGFR_1M17_receptor.pdbqt",
    ]

    for path in possible:
        if path.exists():
            return path

    # Search results/docking
    docking_dir = project_root / "results" / "docking"

    if docking_dir.exists():
        files = list(docking_dir.glob("*receptor*.pdbqt"))

        if files:
            return files[0]

    return None


def run_docking(
    receptor,
    ligand,
    output_file,
    center_x,
    center_y,
    center_z,
    size_x=20,
    size_y=20,
    size_z=20,
):
    """
    Run molecular docking using AutoDock Vina.

    Parameters
    ----------
    receptor : str or Path
        Receptor PDB or PDBQT file.

    ligand : str or Path
        Ligand PDBQT file.

    output_file : str or Path
        Docked output PDBQT file.

    center_x, center_y, center_z : float
        Docking box center.

    size_x, size_y, size_z : float
        Docking box dimensions.

    Returns
    -------
    Path or None
        Output file if docking succeeds.
    """

    print()
    print("=" * 70)
    print("MOLECULAR DOCKING AGENT")
    print("=" * 70)

    project_root = find_project_root()

    print(f"Project root: {project_root}")

    # ---------------------------------------------------------
    # 1. FIND VINA
    # ---------------------------------------------------------

    vina = find_vina(project_root)

    print(f"Vina executable: {vina}")

    if vina is None:
        print()
        print("=" * 70)
        print("ERROR: AUTODOCK VINA NOT FOUND")
        print("=" * 70)
        print()
        print("Expected:")
        print(project_root / "vina" / "vina.exe")
        return None

    # ---------------------------------------------------------
    # 2. FIND RECEPTOR
    # ---------------------------------------------------------

    receptor_pdbqt = find_receptor(project_root, receptor)

    if receptor_pdbqt is None:
        print()
        print("=" * 70)
        print("ERROR: RECEPTOR PDBQT NOT FOUND")
        print("=" * 70)
        print()
        print("Please prepare the receptor first.")
        return None

    print(f"Receptor: {receptor_pdbqt}")

    # ---------------------------------------------------------
    # 3. CHECK LIGAND
    # ---------------------------------------------------------

    ligand = Path(ligand)

    if not ligand.is_absolute():
        ligand = project_root / ligand

    if not ligand.exists():
        print()
        print("=" * 70)
        print("ERROR: LIGAND NOT FOUND")
        print("=" * 70)
        print()
        print(ligand)
        return None

    print(f"Ligand: {ligand}")

    # ---------------------------------------------------------
    # 4. OUTPUT FILE
    # ---------------------------------------------------------

    output_file = Path(output_file)

    if not output_file.is_absolute():
        output_file = project_root / output_file

    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"Output: {output_file}")

    # ---------------------------------------------------------
    # 5. VINA COMMAND
    # ---------------------------------------------------------

    command = [
        str(vina),
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand),
        "--center_x",
        str(center_x),
        "--center_y",
        str(center_y),
        "--center_z",
        str(center_z),
        "--size_x",
        str(size_x),
        "--size_y",
        str(size_y),
        "--size_z",
        str(size_z),
        "--out",
        str(output_file),
    ]

    print()
    print("=" * 70)
    print("STARTING AUTODOCK VINA")
    print("=" * 70)

    print()
    print("Command:")
    print(" ".join(command))

    # ---------------------------------------------------------
    # 6. RUN VINA
    # ---------------------------------------------------------

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            cwd=str(project_root),
        )

    except FileNotFoundError:
        print()
        print("ERROR: vina.exe could not be started.")
        return None

    except PermissionError:
        print()
        print("ERROR: Windows permission error while starting Vina.")
        return None

    except Exception as error:
        print()
        print("DOCKING ERROR")
        print(error)
        return None

    # ---------------------------------------------------------
    # 7. SHOW VINA OUTPUT
    # ---------------------------------------------------------

    print()
    print("=" * 70)
    print("VINA OUTPUT")
    print("=" * 70)

    if result.stdout:
        print(result.stdout)

    if result.stderr:
        print()
        print("=" * 70)
        print("VINA MESSAGE")
        print("=" * 70)
        print(result.stderr)

    # ---------------------------------------------------------
    # 8. CHECK RETURN CODE
    # ---------------------------------------------------------

    print()
    print(f"Vina return code: {result.returncode}")

    if result.returncode != 0:
        print()
        print("=" * 70)
        print("DOCKING FAILED")
        print("=" * 70)
        return None

    # ---------------------------------------------------------
    # 9. CHECK OUTPUT
    # ---------------------------------------------------------

    if not output_file.exists():
        print()
        print("=" * 70)
        print("DOCKING FINISHED BUT OUTPUT WAS NOT CREATED")
        print("=" * 70)
        return None

    if output_file.stat().st_size == 0:
        print()
        print("ERROR: Docking output file is empty.")
        return None

    # ---------------------------------------------------------
    # 10. SUCCESS
    # ---------------------------------------------------------

    print()
    print("=" * 70)
    print("DOCKING SUCCESSFUL")
    print("=" * 70)

    print()
    print(f"Docked file: {output_file}")
    print(f"File size: {output_file.stat().st_size} bytes")

    return output_file


if __name__ == "__main__":

    print()
    print("=" * 70)
    print("DOCKING AGENT TEST")
    print("=" * 70)

    root = find_project_root()

    print(f"Project: {root}")

    vina = find_vina(root)

    if vina:
        print(f"Vina found: {vina}")
    else:
        print("Vina NOT found")

    print()
    print("Docking agent loaded successfully.")