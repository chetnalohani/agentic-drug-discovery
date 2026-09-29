from pathlib import Path
import sys
import pandas as pd

# ============================================================
# PROJECT ROOT
# ============================================================

PROJECT_DIR = Path(__file__).resolve().parent.parent

# Make project root importable when running:
# python tools\admet_analysis.py
sys.path.insert(0, str(PROJECT_DIR))

# ============================================================
# PROJECT IMPORTS
# ============================================================

from agents.compound_discovery_agent import compound_discovery
from agents.compound_filter_agent import check_lipinski

from admet_ai import ADMETModel


# ============================================================
# PATHS
# ============================================================

RESULTS_DIR = PROJECT_DIR / "results"
ADMET_DIR = RESULTS_DIR / "admet"

ADMET_DIR.mkdir(
    parents=True,
    exist_ok=True
)

ADMET_INPUT = ADMET_DIR / "admet_input.csv"
ADMET_OUTPUT = ADMET_DIR / "admet_predictions.csv"


# ============================================================
# SETTINGS
# ============================================================

TARGET = "EGFR"
MAX_COMPOUNDS = 10


# ============================================================
# UTILITY
# ============================================================

def get_value(compound, *keys):
    """
    Return the first non-empty value from a compound dictionary.
    """

    for key in keys:

        value = compound.get(key)

        if value not in (None, ""):
            return value

    return ""


# ============================================================
# CREATE ADMET INPUT
# ============================================================

def create_admet_input():

    print()
    print("=" * 70)
    print("ADMET ANALYSIS")
    print("=" * 70)

    # --------------------------------------------------------
    # 1. COMPOUND DISCOVERY
    # --------------------------------------------------------

    print()
    print("Discovering EGFR compounds...")

    compounds = compound_discovery(
        TARGET,
        max_compounds=MAX_COMPOUNDS
    )

    if not compounds:

        print()
        print("ERROR: No compounds were discovered.")
        return None

    print(
        f"Compounds discovered: {len(compounds)}"
    )

    # --------------------------------------------------------
    # 2. LIPINSKI FILTERING
    # --------------------------------------------------------

    print()
    print("Applying Lipinski filtering...")

    filtered_compounds = []

    for compound in compounds:

        try:

            if check_lipinski(compound):
                filtered_compounds.append(compound)

        except Exception as error:

            print(
                f"Lipinski filtering error: {error}"
            )

    print(
        f"Compounds passing Lipinski: "
        f"{len(filtered_compounds)}"
    )

    if not filtered_compounds:

        print()
        print(
            "ERROR: No compounds passed Lipinski filtering."
        )

        return None

    # --------------------------------------------------------
    # 3. EXTRACT COMPOUND INFORMATION
    # --------------------------------------------------------

    records = []

    for compound in filtered_compounds:

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

        smiles = get_value(
            compound,
            "smiles",
            "canonical_smiles",
            "isomeric_smiles"
        )

        if not smiles:

            print()
            print(
                f"WARNING: No SMILES found for "
                f"{name or cid}. Skipping."
            )

            continue

        records.append(
            {
                "pubchem_cid": cid,
                "compound_name": name,
                "smiles": smiles
            }
        )

    if not records:

        print()
        print(
            "ERROR: No valid SMILES were available."
        )

        return None

    # --------------------------------------------------------
    # 4. SAVE ADMET INPUT
    # --------------------------------------------------------

    input_df = pd.DataFrame(records)

    input_df.to_csv(
        ADMET_INPUT,
        index=False
    )

    print()
    print(
        f"ADMET input saved to:"
    )
    print(ADMET_INPUT)

    print()
    print(
        f"Valid compounds for ADMET: "
        f"{len(input_df)}"
    )

    return input_df


# ============================================================
# RUN ADMET-AI
# ============================================================

def run_admet_prediction(input_df):

    print()
    print("=" * 70)
    print("ADMET-AI PREDICTION")
    print("=" * 70)

    smiles_list = input_df["smiles"].tolist()

    print()
    print(
        f"Running ADMET-AI for "
        f"{len(smiles_list)} compounds..."
    )

    # --------------------------------------------------------
    # Load ADMET-AI model
    # --------------------------------------------------------

    model = ADMETModel()

    # --------------------------------------------------------
    # Generate predictions
    # --------------------------------------------------------

    predictions = model.predict(
        smiles=smiles_list
    )

    # --------------------------------------------------------
    # Convert prediction output to DataFrame
    # --------------------------------------------------------

    if isinstance(predictions, pd.DataFrame):

        prediction_df = predictions.copy()

    else:

        prediction_df = pd.DataFrame(predictions)

    # --------------------------------------------------------
    # Combine compound information with ADMET predictions
    # --------------------------------------------------------

    prediction_df.insert(
        0,
        "pubchem_cid",
        input_df["pubchem_cid"].values
    )

    prediction_df.insert(
        1,
        "compound_name",
        input_df["compound_name"].values
    )

    prediction_df.insert(
        2,
        "smiles",
        input_df["smiles"].values
    )

    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    prediction_df.to_csv(
        ADMET_OUTPUT,
        index=False
    )

    print()
    print(
        "ADMET prediction completed successfully."
    )

    print()
    print(
        "ADMET results saved to:"
    )

    print(ADMET_OUTPUT)

    print()
    print(
        f"Predictions generated for "
        f"{len(prediction_df)} compounds."
    )

    return prediction_df


# ============================================================
# MAIN
# ============================================================

def main():

    input_df = create_admet_input()

    if input_df is None:
        return

    run_admet_prediction(
        input_df
    )

    print()
    print("=" * 70)
    print("ADMET ANALYSIS COMPLETED")
    print("=" * 70)


if __name__ == "__main__":

    main()