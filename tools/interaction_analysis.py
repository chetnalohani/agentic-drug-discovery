from pathlib import Path
import math
import re
import csv


# =========================================================
# PROJECT PATHS
# =========================================================

PROJECT_DIR = Path(__file__).resolve().parent.parent

DOCKING_DIR = PROJECT_DIR / "results" / "docking"

RECEPTOR_FILE = (
    DOCKING_DIR /
    "EGFR_1M17_receptor_clean.pdb"
)

OUTPUT_DIR = (
    DOCKING_DIR /
    "interaction_analysis"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# =========================================================
# COMPOUND NAMES
# =========================================================

COMPOUND_NAMES = {

    "123631": "gefitinib",

    "176870": "erlotinib",

    "10184653": "afatinib",

    "71496458": "osimertinib",

    "208908": "lapatinib",

    "3081361": "vandetanib",

    "156414": "canertinib",

    "22024915": "icotinib",

    "6445562": "pelitinib",

    "25127713": "poziotinib"
}


# =========================================================
# DISTANCE
# =========================================================

def distance(atom1, atom2):

    return math.sqrt(

        (atom1["x"] - atom2["x"]) ** 2

        +

        (atom1["y"] - atom2["y"]) ** 2

        +

        (atom1["z"] - atom2["z"]) ** 2

    )


# =========================================================
# GET ELEMENT
# =========================================================

def get_element(atom_name):

    letters = re.sub(
        r"[^A-Za-z]",
        "",
        atom_name
    ).upper()

    if not letters:
        return ""

    if letters.startswith("CL"):
        return "CL"

    if letters.startswith("BR"):
        return "BR"

    return letters[0]


# =========================================================
# READ RECEPTOR
# =========================================================

def read_receptor(file_path):

    atoms = []

    with open(
        file_path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as file:

        for line in file:

            if not line.startswith(
                ("ATOM", "HETATM")
            ):
                continue

            try:

                atom_name = (
                    line[12:16].strip()
                )

                element = ""

                if len(line) >= 78:

                    element = (
                        line[76:78]
                        .strip()
                        .upper()
                    )

                if not element:

                    element = get_element(
                        atom_name
                    )

                atom = {

                    "atom": atom_name,

                    "resname":
                        line[17:20].strip(),

                    "chain":
                        line[21].strip()
                        or "-",

                    "resnum":
                        line[22:26].strip(),

                    "x":
                        float(line[30:38]),

                    "y":
                        float(line[38:46]),

                    "z":
                        float(line[46:54]),

                    "element":
                        element
                }

                atoms.append(atom)

            except (
                ValueError,
                IndexError
            ):

                continue

    return atoms


# =========================================================
# READ BEST DOCKING POSE
# =========================================================

def read_best_pose(file_path):

    affinity = None

    ligand_atoms = []

    in_first_model = False

    saw_model = False

    with open(
        file_path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as file:

        for line in file:

            # ---------------------------------------------
            # VINA AFFINITY
            # ---------------------------------------------

            match = re.search(

                r"REMARK VINA RESULT:\s*"
                r"([-+]?\d+(?:\.\d+)?)",

                line
            )

            if (
                match
                and affinity is None
            ):

                affinity = float(
                    match.group(1)
                )

            # ---------------------------------------------
            # MODEL 1
            # ---------------------------------------------

            if line.startswith("MODEL"):

                if not saw_model:

                    saw_model = True

                    in_first_model = True

                else:

                    in_first_model = False

                    break

                continue

            # ---------------------------------------------
            # END MODEL
            # ---------------------------------------------

            if line.startswith("ENDMDL"):

                if in_first_model:

                    break

                continue

            # ---------------------------------------------
            # ATOMS
            # ---------------------------------------------

            if not line.startswith(
                ("ATOM", "HETATM")
            ):

                continue

            if (
                saw_model
                and not in_first_model
            ):

                continue

            try:

                atom_name = (
                    line[12:16].strip()
                )

                element = ""

                if len(line) >= 78:

                    element = (
                        line[76:78]
                        .strip()
                        .upper()
                    )

                if not element:

                    element = get_element(
                        atom_name
                    )

                atom = {

                    "atom": atom_name,

                    "x":
                        float(line[30:38]),

                    "y":
                        float(line[38:46]),

                    "z":
                        float(line[46:54]),

                    "element":
                        element
                }

                ligand_atoms.append(
                    atom
                )

            except (
                ValueError,
                IndexError
            ):

                continue

    return (
        affinity,
        ligand_atoms
    )


# =========================================================
# MAIN
# =========================================================

def main():

    print()
    print("=" * 70)
    print(
        "EGFR - 1M17 DOCKING INTERACTION ANALYSIS"
    )
    print("=" * 70)

    # -----------------------------------------------------
    # CHECK RECEPTOR
    # -----------------------------------------------------

    if not RECEPTOR_FILE.exists():

        print()
        print(
            "ERROR: Receptor file not found:"
        )

        print(
            RECEPTOR_FILE
        )

        return

    # -----------------------------------------------------
    # READ RECEPTOR
    # -----------------------------------------------------

    receptor_atoms = read_receptor(
        RECEPTOR_FILE
    )

    if not receptor_atoms:

        print()
        print(
            "ERROR: No receptor atoms found."
        )

        return

    print()
    print(
        f"Receptor atoms read: "
        f"{len(receptor_atoms)}"
    )

    # -----------------------------------------------------
    # FIND DOCKED FILES
    # -----------------------------------------------------

    docking_files = sorted(
        DOCKING_DIR.glob(
            "docked_*.pdbqt"
        )
    )

    if not docking_files:

        print()
        print(
            "ERROR: No docked PDBQT files found."
        )

        return

    print(
        f"Docked files found: "
        f"{len(docking_files)}"
    )

    print()

    summary_rows = []

    # =====================================================
    # ANALYZE EACH COMPOUND
    # =====================================================

    for docking_file in docking_files:

        # -------------------------------------------------
        # GET CID
        # -------------------------------------------------

        match = re.search(

            r"docked_(\d+)\.pdbqt$",

            docking_file.name,

            re.IGNORECASE
        )

        if match:

            cid = match.group(1)

        else:

            cid = (
                docking_file.stem
                .replace(
                    "docked_",
                    ""
                )
            )

        compound = (
            COMPOUND_NAMES.get(
                cid,
                cid
            )
        )

        print(
            f"Analyzing: "
            f"{compound}"
        )

        # -------------------------------------------------
        # READ DOCKED POSE
        # -------------------------------------------------

        affinity, ligand_atoms = (
            read_best_pose(
                docking_file
            )
        )

        report_file = (
            OUTPUT_DIR /
            f"{compound}_{cid}_interaction_report.txt"
        )

        if not ligand_atoms:

            report_file.write_text(

                f"Compound: {compound}\n"
                f"CID: {cid}\n"
                f"ERROR: No ligand atoms found.\n",

                encoding="utf-8"
            )

            print(
                "  ERROR: No ligand atoms."
            )

            continue

        # -------------------------------------------------
        # FIND CONTACTS
        # -------------------------------------------------

        contacts = []

        for receptor_atom in receptor_atoms:

            closest_distance = min(

                distance(
                    receptor_atom,
                    ligand_atom
                )

                for ligand_atom
                in ligand_atoms

            )

            if closest_distance <= 4.0:

                contacts.append(

                    (
                        closest_distance,
                        receptor_atom
                    )

                )

        contacts.sort(
            key=lambda x: x[0]
        )

        # -------------------------------------------------
        # GROUP BY RESIDUE
        # -------------------------------------------------

        residue_data = {}

        for (
            atom_distance,
            receptor_atom
        ) in contacts:

            key = (

                receptor_atom["chain"],

                receptor_atom["resnum"],

                receptor_atom["resname"]

            )

            if (

                key not in residue_data

                or

                atom_distance
                < residue_data[key]

            ):

                residue_data[key] = (
                    atom_distance
                )

        residues = sorted(

            [

                (
                    atom_distance,
                    chain,
                    resnum,
                    resname
                )

                for (

                    (
                        chain,
                        resnum,
                        resname
                    ),

                    atom_distance

                ) in residue_data.items()

            ],

            key=lambda x: x[0]
        )

        # -------------------------------------------------
        # CLOSE N/O/S PAIRS
        # -------------------------------------------------

        close_pairs = []

        for receptor_atom in receptor_atoms:

            if receptor_atom[
                "element"
            ] not in {

                "N",
                "O",
                "S"

            }:

                continue

            for ligand_atom in ligand_atoms:

                if ligand_atom[
                    "element"
                ] not in {

                    "N",
                    "O",
                    "S"

                }:

                    continue

                atom_distance = distance(

                    receptor_atom,

                    ligand_atom

                )

                if atom_distance <= 3.5:

                    close_pairs.append(

                        (

                            atom_distance,

                            receptor_atom,

                            ligand_atom

                        )

                    )

        close_pairs.sort(
            key=lambda x: x[0]
        )

        # -------------------------------------------------
        # WRITE REPORT
        # -------------------------------------------------

        report = []

        report.append(
            "EGFR - 1M17 DOCKING INTERACTION REPORT"
        )

        report.append(
            "=" * 60
        )

        report.append(
            f"Compound: {compound}"
        )

        report.append(
            f"PubChem CID: {cid}"
        )

        report.append(
            f"Docked file: "
            f"{docking_file.name}"
        )

        if affinity is not None:

            report.append(
                f"Vina best affinity: "
                f"{affinity:.3f} kcal/mol"
            )

        else:

            report.append(
                "Vina best affinity: "
                "Not found"
            )

        report.append("")

        report.append(
            "BINDING-SITE RESIDUES "
            "WITHIN 4.0 A"
        )

        report.append(
            "-" * 60
        )

        if residues:

            for (

                atom_distance,
                chain,
                resnum,
                resname

            ) in residues:

                report.append(

                    f"{resname} "
                    f"{chain}{resnum} "
                    f"- "
                    f"{atom_distance:.2f} A"

                )

        else:

            report.append(
                "No residues found."
            )

        report.append("")

        report.append(
            "CLOSE N/O/S ATOM PAIRS "
            "WITHIN 3.5 A"
        )

        report.append(
            "-" * 60
        )

        report.append(
            "These are proximity-based "
            "interaction candidates."
        )

        report.append(
            "They are NOT confirmed hydrogen bonds."
        )

        report.append("")

        if close_pairs:

            for (

                atom_distance,
                receptor_atom,
                ligand_atom

            ) in close_pairs[:20]:

                report.append(

                    f"{receptor_atom['resname']} "
                    f"{receptor_atom['chain']}"
                    f"{receptor_atom['resnum']} "
                    f"{receptor_atom['atom']} "
                    f"- ligand "
                    f"{ligand_atom['atom']} "
                    f": "
                    f"{atom_distance:.2f} A"

                )

        else:

            report.append(
                "No close N/O/S pairs found."
            )

        # -------------------------------------------------
        # SAVE REPORT
        # -------------------------------------------------

        report_file.write_text(

            "\n".join(report) + "\n",

            encoding="utf-8"
        )

        # -------------------------------------------------
        # SUMMARY DATA
        # -------------------------------------------------

        residue_names = "; ".join(

            (

                f"{resname} "
                f"{chain}{resnum}"

            )

            for (

                atom_distance,
                chain,
                resnum,
                resname

            ) in residues

        )

        summary_rows.append({

            "compound":
                compound,

            "cid":
                cid,

            "vina_affinity_kcal_mol":
                (
                    affinity
                    if affinity is not None
                    else ""
                ),

            "binding_site_residues":
                residue_names,

            "report_file":
                str(report_file)

        })

        if affinity is not None:

            print(

                f"  Vina affinity: "
                f"{affinity:.3f} kcal/mol"

            )

        print(
            f"  Report saved: "
            f"{report_file.name}"
        )

        print()

    # =====================================================
    # SAVE CSV SUMMARY
    # =====================================================

    csv_file = (
        OUTPUT_DIR /
        "interaction_summary.csv"
    )

    with open(

        csv_file,
        "w",
        newline="",
        encoding="utf-8"

    ) as file:

        writer = csv.DictWriter(

            file,

            fieldnames=[

                "compound",

                "cid",

                "vina_affinity_kcal_mol",

                "binding_site_residues",

                "report_file"

            ]

        )

        writer.writeheader()

        writer.writerows(
            summary_rows
        )

    # =====================================================
    # FINISHED
    # =====================================================

    print()
    print("=" * 70)
    print(
        "INTERACTION ANALYSIS COMPLETED"
    )
    print("=" * 70)

    print()
    print(
        f"Reports folder:"
    )

    print(
        OUTPUT_DIR
    )

    print()
    print(
        f"Summary CSV:"
    )

    print(
        csv_file
    )


# =========================================================
# RUN
# =========================================================

if __name__ == "__main__":

    main()