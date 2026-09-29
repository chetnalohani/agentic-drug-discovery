import subprocess
from pathlib import Path


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
        Receptor PDBQT file.

    ligand : str or Path
        Ligand PDBQT file.

    output_file : str or Path
        Output PDBQT file for docked poses.

    center_x, center_y, center_z : float
        Center of the docking box.

    size_x, size_y, size_z : float
        Size of the docking box.
    """

    print()
    print("=" * 70)
    print("MOLECULAR DOCKING AGENT")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. FIND PROJECT ROOT
    # ---------------------------------------------------------

    project_root = Path(__file__).resolve().parent.parent

    print(f"Project root: {project_root}")

    # ---------------------------------------------------------
    # 2. FIND AUTODOCK VINA
    # ---------------------------------------------------------

    vina_executable = project_root / "vina" / "vina.exe"

    print(f"Vina executable: {vina_executable}")

    if not vina_executable.exists():
        print()
        print("ERROR: AutoDock Vina was not found.")
        print(f"Expected location:")
        print(vina_executable)
        print()
        print("Make sure vina.exe is inside:")
        print(project_root / "vina")
        return None

    # ---------------------------------------------------------
    # 3. CONVERT PATHS TO ABSOLUTE PATHS
    # ---------------------------------------------------------

    receptor = Path(receptor)
    ligand = Path(ligand)
    output_file = Path(output_file)

    if not receptor.is_absolute():
        receptor = project_root / receptor

    if not ligand.is_absolute():
        ligand = project_root / ligand

    if not output_file.is_absolute():
        output_file = project_root / output_file

    receptor = receptor.resolve()
    ligand = ligand.resolve()
    output_file = output_file.resolve()

    # ---------------------------------------------------------
    # 4. PRINT INPUT INFORMATION
    # ---------------------------------------------------------

    print()
    print(f"Receptor : {receptor}")
    print(f"Ligand   : {ligand}")
    print(f"Output   : {output_file}")

    # ---------------------------------------------------------
    # 5. CHECK RECEPTOR
    # ---------------------------------------------------------

    if not receptor.exists():
        print()
        print("ERROR: Receptor file not found.")
        print(f"Expected receptor:")
        print(receptor)
        return None

    # ---------------------------------------------------------
    # 6. CHECK LIGAND
    # ---------------------------------------------------------

    if not ligand.exists():
        print()
        print("ERROR: Ligand file not found.")
        print(f"Expected ligand:")
        print(ligand)
        return None

    # ---------------------------------------------------------
    # 7. CHECK FILE TYPES
    # ---------------------------------------------------------

    if receptor.suffix.lower() != ".pdbqt":
        print()
        print("WARNING:")
        print("The receptor file is not a .pdbqt file.")
        print(f"Current receptor: {receptor.name}")
        print()
        print("AutoDock Vina normally expects a prepared PDBQT receptor.")

    if ligand.suffix.lower() != ".pdbqt":
        print()
        print("WARNING:")
        print("The ligand file is not a .pdbqt file.")
        print(f"Current ligand: {ligand.name}")
        print()
        print("AutoDock Vina normally expects a prepared PDBQT ligand.")

    # ---------------------------------------------------------
    # 8. CREATE OUTPUT DIRECTORY
    # ---------------------------------------------------------

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------
    # 9. BUILD VINA COMMAND
    # ---------------------------------------------------------

    command = [
        str(vina_executable),

        "--receptor",
        str(receptor),

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

    # ---------------------------------------------------------
    # 10. SHOW DOCKING PARAMETERS
    # ---------------------------------------------------------

    print()
    print("-" * 70)
    print("DOCKING PARAMETERS")
    print("-" * 70)

    print(f"Center X : {center_x}")
    print(f"Center Y : {center_y}")
    print(f"Center Z : {center_z}")

    print(f"Size X   : {size_x}")
    print(f"Size Y   : {size_y}")
    print(f"Size Z   : {size_z}")

    print()
    print("Starting AutoDock Vina...")
    print()

    # ---------------------------------------------------------
    # 11. RUN AUTODOCK VINA
    # ---------------------------------------------------------

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            cwd=str(project_root),
        )

        # -----------------------------------------------------
        # 12. PRINT VINA OUTPUT
        # -----------------------------------------------------

        if result.stdout:
            print(result.stdout)

        if result.stderr:
            print()
            print("VINA MESSAGE:")
            print(result.stderr)

        # -----------------------------------------------------
        # 13. CHECK RETURN CODE
        # -----------------------------------------------------

        if result.returncode != 0:

            print()
            print("=" * 70)
            print("DOCKING FAILED")
            print("=" * 70)

            print(f"Vina return code: {result.returncode}")

            return None

        # -----------------------------------------------------
        # 14. CHECK OUTPUT FILE
        # -----------------------------------------------------

        if not output_file.exists():

            print()
            print("=" * 70)
            print("DOCKING FINISHED BUT OUTPUT WAS NOT CREATED")
            print("=" * 70)

            print(f"Expected output:")
            print(output_file)

            return None

        # -----------------------------------------------------
        # 15. SUCCESS
        # -----------------------------------------------------

        print()
        print("=" * 70)
        print("DOCKING COMPLETED SUCCESSFULLY")
        print("=" * 70)

        print()
        print(f"Docked output:")
        print(output_file)

        print()
        print(f"Output file size: {output_file.stat().st_size} bytes")

        return output_file

    # ---------------------------------------------------------
    # 16. HANDLE PYTHON ERRORS
    # ---------------------------------------------------------

    except FileNotFoundError:

        print()
        print("=" * 70)
        print("AUTODOCK VINA COULD NOT BE STARTED")
        print("=" * 70)

        print("vina.exe could not be executed.")

        return None

    except Exception as error:

        print()
        print("=" * 70)
        print("DOCKING ERROR")
        print("=" * 70)

        print(f"{type(error).__name__}: {error}")

        return None


# =============================================================
# DIRECT TEST
# =============================================================

if __name__ == "__main__":

    print()
    print("=" * 70)
    print("STAGE 6: MOLECULAR DOCKING")
    print("=" * 70)

    print()
    print("Docking agent is ready.")
    print()
    print("Use run_docking() from main.py to perform docking.")