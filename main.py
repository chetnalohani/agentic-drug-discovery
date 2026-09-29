from pathlib import Path
import re
import subprocess
import sys

from agents.target_validation_agent import target_validation
from agents.literature_agent import literature_research
from agents.compound_discovery_agent import compound_discovery
from agents.compound_filter_agent import check_lipinski
from agents.docking_agent import run_docking


# =========================================================
# PROJECT PATHS
# =========================================================

PROJECT_DIR = Path(__file__).resolve().parent

RECEPTOR_PDB = PROJECT_DIR / "1M17.pdb"

RESULTS_DIR = PROJECT_DIR / "results"
DOCKING_DIR = RESULTS_DIR / "docking"

RECEPTOR_PDBQT = DOCKING_DIR / "EGFR_1M17_receptor.pdbqt"


# =========================================================
# VALIDATED 1M17 DOCKING BOX
# =========================================================

CENTER_X = 22.839
CENTER_Y = 6.163
CENTER_Z = 48.007

SIZE_X = 22.0
SIZE_Y = 16.0
SIZE_Z = 16.0


# =========================================================
# UTILITY FUNCTIONS
# =========================================================

def safe_filename(text):
    text = str(text or "compound")
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)[:80]


def get_value(compound, *keys):
    for key in keys:
        value = compound.get(key)

        if value not in (None, ""):
            return value

    return ""


# =========================================================
# CLEAN RECEPTOR
# =========================================================

def clean_receptor(source, destination):

    print("Cleaning receptor...")

    lines = source.read_text(
        encoding="utf-8",
        errors="ignore"
    ).splitlines(True)

    cleaned = []
    removed_aq4 = 0

    for line in lines:

        if len(line) >= 20:

            record = line[0:6].strip()

            if record in ("ATOM", "HETATM"):

                resname = line[17:20].strip()

                if resname == "AQ4":
                    removed_aq4 += 1
                    continue

        cleaned.append(line)

    destination.write_text(
        "".join(cleaned),
        encoding="utf-8"
    )

    return removed_aq4


# =========================================================
# FIND MEEKO
# =========================================================

def find_meeko():

    commands = [
        "mk_prepare_receptor.exe",
        "mk_prepare_receptor.py",
        "mk_prepare_receptor"
    ]

    for command in commands:

        try:

            result = subprocess.run(
                [command, "--help"],
                capture_output=True,
                text=True
            )

            if result.returncode in (0, 1, 2):
                return command

        except (FileNotFoundError, OSError):
            pass

    scripts_dir = Path(sys.executable).resolve().parent

    for name in [
        "mk_prepare_receptor.exe",
        "mk_prepare_receptor.py"
    ]:

        candidate = scripts_dir / name

        if candidate.exists():
            return str(candidate)

    return None


# =========================================================
# RECEPTOR PREPARATION
# =========================================================

def prepare_receptor():

    print()
    print("=" * 70)
    print("RECEPTOR PREPARATION")
    print("=" * 70)

    DOCKING_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not RECEPTOR_PDB.exists():

        print("ERROR: 1M17.pdb not found.")
        print(RECEPTOR_PDB)

        return None

    print()
    print("Input receptor:")
    print(RECEPTOR_PDB)

    # -----------------------------------------------------
    # USE EXISTING RECEPTOR
    # -----------------------------------------------------

    if RECEPTOR_PDBQT.exists():

        print()
        print("Prepared receptor already exists.")

        return RECEPTOR_PDBQT

    # -----------------------------------------------------
    # CLEAN RECEPTOR
    # -----------------------------------------------------

    cleaned = (
        DOCKING_DIR /
        "EGFR_1M17_receptor_clean.pdb"
    )

    try:

        removed = clean_receptor(
            RECEPTOR_PDB,
            cleaned
        )

    except Exception as error:

        print()
        print("ERROR: Receptor cleaning failed.")
        print(error)

        return None

    print()
    print("Cleaned receptor:")
    print(cleaned)

    print(
        f"Removed AQ4 atoms: {removed}"
    )

    # -----------------------------------------------------
    # FIND MEEKO
    # -----------------------------------------------------

    meeko = find_meeko()

    if meeko is None:

        print()
        print("ERROR: Meeko was not found.")

        print()
        print("Python:")
        print(sys.executable)

        return None

    # -----------------------------------------------------
    # PREPARE RECEPTOR WITH MEEKO
    # -----------------------------------------------------

    output_base = (
        DOCKING_DIR /
        "EGFR_1M17_receptor"
    )

    command = [
        meeko,
        "-i",
        str(cleaned),
        "-o",
        str(output_base),
        "-p",
        "--default_altloc",
        "A"
    ]

    print()
    print("Preparing receptor with Meeko...")

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            cwd=str(PROJECT_DIR)
        )

    except Exception as error:

        print()
        print("ERROR running Meeko:")
        print(error)

        return None

    if result.stdout:
        print(result.stdout)

    if result.stderr:
        print(result.stderr)

    if result.returncode != 0:

        print()
        print("RECEPTOR PREPARATION FAILED")
        print(
            f"Return code: {result.returncode}"
        )

        return None

    if not RECEPTOR_PDBQT.exists():

        print()
        print(
            "ERROR: Receptor PDBQT was not created."
        )

        print(RECEPTOR_PDBQT)

        return None

    print()
    print("RECEPTOR PREPARATION COMPLETED")

    return RECEPTOR_PDBQT


# =========================================================
# LIGAND PREPARATION
# =========================================================

def prepare_ligand(compound, index):

    smiles = get_value(
        compound,
        "smiles",
        "canonical_smiles",
        "isomeric_smiles"
    )

    if not smiles:

        print(
            f"Compound {index}: No SMILES found."
        )

        return None

    cid = get_value(
        compound,
        "cid",
        "pubchem_cid",
        "pubchem_id"
    )

    name = get_value(
        compound,
        "name",
        "compound_name",
        "iupac_name"
    )

    identifier = (
        str(cid)
        if cid
        else safe_filename(name)
    )

    if not identifier:
        identifier = f"compound_{index}"

    ligand_file = (
        DOCKING_DIR /
        f"ligand_{identifier}.pdbqt"
    )

    print()
    print(
        f"Preparing ligand {index}: {name}"
    )

    print(
        f"CID: {cid}"
    )

    if ligand_file.exists():

        print(
            "Using existing ligand:"
        )

        print(ligand_file)

        return ligand_file

    # -----------------------------------------------------
    # RDKit + MEEKO
    # -----------------------------------------------------

    preparation_code = r'''
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

from meeko import (
    MoleculePreparation,
    PDBQTWriterLegacy
)


smiles = sys.argv[1]
output_file = sys.argv[2]


# ---------------------------------------------------------
# CREATE MOLECULE FROM SMILES
# ---------------------------------------------------------

mol = Chem.MolFromSmiles(smiles)

if mol is None:

    raise RuntimeError(
        "RDKit could not parse the SMILES."
    )


# ---------------------------------------------------------
# ADD HYDROGENS
# ---------------------------------------------------------

mol = Chem.AddHs(mol)


# ---------------------------------------------------------
# GENERATE 3D COORDINATES
# ---------------------------------------------------------

print("Generating 3D coordinates...")


status = AllChem.EmbedMolecule(
    mol,
    AllChem.ETKDGv3()
)


if status != 0:

    print(
        "Standard embedding failed."
    )

    status = AllChem.EmbedMolecule(
        mol,
        randomSeed=42
    )


if status != 0:

    raise RuntimeError(
        "RDKit could not generate 3D coordinates."
    )


# ---------------------------------------------------------
# ENERGY MINIMIZATION
# ---------------------------------------------------------

try:

    AllChem.UFFOptimizeMolecule(
        mol,
        maxIters=200
    )

except Exception:

    print(
        "UFF minimization failed; "
        "continuing with embedded structure."
    )


# ---------------------------------------------------------
# MEEKO PREPARATION
# ---------------------------------------------------------

print(
    "Preparing ligand with Meeko..."
)


preparation = MoleculePreparation()


setups = preparation.prepare(mol)


if not setups:

    raise RuntimeError(
        "Meeko produced no ligand setup."
    )


# ---------------------------------------------------------
# WRITE PDBQT
# ---------------------------------------------------------

pdbqt_string, is_ok, error = (
    PDBQTWriterLegacy.write_string(
        setups[0]
    )
)


if not is_ok:

    raise RuntimeError(
        str(error)
    )


with open(
    output_file,
    "w",
    encoding="utf-8"
) as handle:

    handle.write(pdbqt_string)


print(
    "Ligand PDBQT created:"
)

print(output_file)
'''

    # -----------------------------------------------------
    # RUN LIGAND PREPARATION
    # -----------------------------------------------------

    try:

        result = subprocess.run(
            [
                sys.executable,
                "-c",
                preparation_code,
                smiles,
                str(ligand_file)
            ],
            capture_output=True,
            text=True,
            cwd=str(PROJECT_DIR)
        )

    except Exception as error:

        print()
        print(
            "Ligand preparation error:"
        )

        print(error)

        return None

    if result.stdout:
        print(result.stdout)

    if result.stderr:
        print(result.stderr)

    if (
        result.returncode != 0
        or not ligand_file.exists()
    ):

        print()
        print(
            "LIGAND PREPARATION FAILED"
        )

        return None

    print()
    print(
        "Ligand preparation completed."
    )

    return ligand_file


# =========================================================
# DOCKING
# =========================================================

def dock_compound(
    compound,
    index,
    receptor
):

    name = get_value(
        compound,
        "name",
        "compound_name",
        "iupac_name"
    )

    cid = get_value(
        compound,
        "cid",
        "pubchem_cid",
        "pubchem_id"
    )

    ligand = prepare_ligand(
        compound,
        index
    )

    if ligand is None:
        return None

    identifier = (
        str(cid)
        if cid
        else safe_filename(name)
    )

    if not identifier:
        identifier = f"compound_{index}"

    output = (
        DOCKING_DIR /
        f"docked_{identifier}.pdbqt"
    )

    print()
    print(
        f"Docking compound {index}: {name}"
    )

    return run_docking(
        receptor=str(receptor),
        ligand=str(ligand),
        output_file=str(output),
        center_x=CENTER_X,
        center_y=CENTER_Y,
        center_z=CENTER_Z,
        size_x=SIZE_X,
        size_y=SIZE_Y,
        size_z=SIZE_Z
    )


# =========================================================
# MAIN PIPELINE
# =========================================================

def main():

    target = "EGFR"

    print()
    print("=" * 70)
    print("AGENTIC DRUG DISCOVERY SYSTEM")
    print("=" * 70)

    print()
    print(f"Target: {target}")

    # =====================================================
    # STAGE 1
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 1 - TARGET VALIDATION")
    print("=" * 70)

    target_data = target_validation(target)

    if not target_data:

        print(
            "Target validation failed."
        )

        return

    print(
        "Target validation completed."
    )

    # =====================================================
    # STAGE 2
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 2 - LITERATURE RESEARCH")
    print("=" * 70)

    papers = literature_research(
        target,
        max_results=10
    )

    print()
    print(
        f"Literature records: "
        f"{len(papers) if papers else 0}"
    )

    # =====================================================
    # STAGE 3
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 3 - COMPOUND DISCOVERY")
    print("=" * 70)

    compounds = compound_discovery(
        target,
        max_compounds=10
    )

    if not compounds:

        print(
            "No compounds found."
        )

        return

    print()
    print(
        f"Compounds discovered: "
        f"{len(compounds)}"
    )

    # =====================================================
    # STAGE 4
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 4 - COMPOUND FILTERING")
    print("=" * 70)

    filtered = []

    for compound in compounds:

        try:

            if check_lipinski(compound):

                filtered.append(compound)

        except Exception as error:

            print()
            print(
                f"Lipinski error: {error}"
            )

    print()
    print(
        f"Compounds passing Lipinski: "
        f"{len(filtered)}"
    )

    if not filtered:

        print(
            "No compounds passed filtering."
        )

        return

    # =====================================================
    # STAGE 5
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 5 - FILTERED CANDIDATES")
    print("=" * 70)

    for i, compound in enumerate(
        filtered,
        start=1
    ):

        name = get_value(
            compound,
            "name",
            "compound_name"
        )

        cid = get_value(
            compound,
            "cid",
            "pubchem_cid"
        )

        mw = compound.get(
            "molecular_weight",
            ""
        )

        print()
        print(
            f"{i}. {name}"
        )

        print(
            f"   CID: {cid}"
        )

        print(
            f"   MW: {mw}"
        )

    # =====================================================
    # STAGE 6A
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 6A - EGFR RECEPTOR PREPARATION")
    print("=" * 70)

    receptor = prepare_receptor()

    if receptor is None:

        print()

        print(
            "Docking cannot continue because "
            "receptor preparation failed."
        )

        return

    # =====================================================
    # STAGE 6B
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 6B - MOLECULAR DOCKING")
    print("=" * 70)

    docking_results = []

    for index, compound in enumerate(
        filtered,
        start=1
    ):

        try:

            result = dock_compound(
                compound,
                index,
                receptor
            )

            if result is not None:

                docking_results.append(
                    {
                        "compound": compound,
                        "result": result
                    }
                )

        except Exception as error:

            print()

            print(
                f"Docking error for "
                f"compound {index}:"
            )

            print(error)

    # =====================================================
    # STAGE 7
    # =====================================================

    print()
    print("=" * 70)
    print("STAGE 7 - DOCKING SUMMARY")
    print("=" * 70)

    print()

    print(
        f"Filtered compounds: "
        f"{len(filtered)}"
    )

    print(
        f"Successful docking runs: "
        f"{len(docking_results)}"
    )

    for i, item in enumerate(
        docking_results,
        start=1
    ):

        compound = item["compound"]

        name = get_value(
            compound,
            "name",
            "compound_name"
        )

        cid = get_value(
            compound,
            "cid",
            "pubchem_cid"
        )

        print()

        print(
            f"{i}. {name}"
        )

        print(
            f"   CID: {cid}"
        )

        print(
            f"   Output: "
            f"{item['result']}"
        )

    # =====================================================
    # COMPLETE
    # =====================================================

    print()
    print("=" * 70)
    print(
        "AGENTIC DRUG DISCOVERY PIPELINE COMPLETED"
    )
    print("=" * 70)


# =========================================================
# RUN
# =========================================================

if __name__ == "__main__":
    main()