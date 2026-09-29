import os
import pandas as pd


# ============================================================
# PATHS
# ============================================================

PROJECT_DIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

ADMET_FILE = os.path.join(
    PROJECT_DIR,
    "results",
    "admet",
    "admet_predictions.csv"
)

OUTPUT_DIR = os.path.join(
    PROJECT_DIR,
    "results",
    "final"
)

OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    "final_candidate_ranking.csv"
)


# ============================================================
# DOCKING SCORES
# AutoDock Vina results already obtained in this project
# More negative = better predicted binding affinity
# ============================================================

DOCKING_SCORES = {
    25127713: -9.012,   # poziotinib
    208908: -8.719,     # lapatinib
    6445562: -8.638,    # pelitinib
    156414: -8.544,     # canertinib
    3081361: -8.362,    # vandetanib
    123631: -8.320,     # gefitinib
    71496458: -8.243,   # osimertinib
    10184653: -8.138,   # afatinib
    176870: -7.015,     # erlotinib
    22024915: -6.156     # icotinib
}


# ============================================================
# MAIN FUNCTION
# ============================================================

def main():

    print()
    print("=" * 70)
    print("FINAL CANDIDATE INTEGRATION")
    print("=" * 70)

    # --------------------------------------------------------
    # Check ADMET file
    # --------------------------------------------------------

    if not os.path.exists(ADMET_FILE):
        raise FileNotFoundError(
            f"ADMET results not found:\n{ADMET_FILE}"
        )

    # --------------------------------------------------------
    # Create output directory
    # --------------------------------------------------------

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Read ADMET predictions
    # --------------------------------------------------------

    df = pd.read_csv(ADMET_FILE)

    if df.empty:
        raise ValueError(
            "ADMET prediction file is empty."
        )

    # --------------------------------------------------------
    # Required column check
    # --------------------------------------------------------

    required_columns = [
        "pubchem_cid",
        "compound_name"
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing_columns)
        )

    # --------------------------------------------------------
    # Convert CID to integer
    # --------------------------------------------------------

    df["pubchem_cid"] = pd.to_numeric(
        df["pubchem_cid"],
        errors="coerce"
    )

    if df["pubchem_cid"].isna().any():
        raise ValueError(
            "One or more PubChem CIDs could not be read."
        )

    df["pubchem_cid"] = df[
        "pubchem_cid"
    ].astype(int)

    # --------------------------------------------------------
    # Add docking scores
    # --------------------------------------------------------

    df["docking_score"] = df[
        "pubchem_cid"
    ].map(DOCKING_SCORES)

    if df["docking_score"].isna().any():

        missing_cids = df.loc[
            df["docking_score"].isna(),
            "pubchem_cid"
        ].tolist()

        raise ValueError(
            "Docking score not available for CID(s): "
            + ", ".join(map(str, missing_cids))
        )

    # --------------------------------------------------------
    # Convert docking score to 0-100 score
    #
    # More negative Vina score = better.
    # Best score receives 100.
    # Worst score receives 0.
    # --------------------------------------------------------

    best_docking = df["docking_score"].min()
    worst_docking = df["docking_score"].max()

    if best_docking == worst_docking:

        df["docking_score_normalized"] = 100.0

    else:

        df["docking_score_normalized"] = (
            (worst_docking - df["docking_score"])
            / (worst_docking - best_docking)
        ) * 100.0

    # --------------------------------------------------------
    # QED score
    # QED is already between 0 and 1.
    # Convert to 0-100.
    # --------------------------------------------------------

    if "QED" in df.columns:

        df["qed_score"] = (
            pd.to_numeric(
                df["QED"],
                errors="coerce"
            )
            .fillna(0)
            * 100.0
        )

    else:

        df["qed_score"] = 0.0

    # --------------------------------------------------------
    # Lipinski score
    # --------------------------------------------------------

    if "Lipinski" in df.columns:

        lipinski_numeric = pd.to_numeric(
            df["Lipinski"],
            errors="coerce"
        )

        df["lipinski_score"] = (
            lipinski_numeric
            .fillna(0)
            .clip(0, 1)
            * 100.0
        )

    else:

        df["lipinski_score"] = 0.0

    # --------------------------------------------------------
    # Integrated score
    #
    # Docking       = 50%
    # QED           = 25%
    # Lipinski      = 25%
    #
    # ADMET predictions are retained in the final table
    # and are NOT incorrectly converted into a single
    # safety score because different ADMET endpoints have
    # different biological meanings and directions.
    # --------------------------------------------------------

    df["final_score"] = (
        0.50 * df["docking_score_normalized"]
        + 0.25 * df["qed_score"]
        + 0.25 * df["lipinski_score"]
    )

    df["final_score"] = df[
        "final_score"
    ].round(2)

    # --------------------------------------------------------
    # Rank compounds
    # --------------------------------------------------------

    df = df.sort_values(
        by="final_score",
        ascending=False
    ).reset_index(drop=True)

    df.insert(
        0,
        "final_rank",
        range(1, len(df) + 1)
    )

    # --------------------------------------------------------
    # Save final ranking
    # --------------------------------------------------------

    df.to_csv(
        OUTPUT_FILE,
        index=False
    )

    # --------------------------------------------------------
    # Display results
    # --------------------------------------------------------

    print()
    print("FINAL CANDIDATE RANKING")
    print("=" * 70)

    for _, row in df.iterrows():

        print()
        print(
            f"Rank {int(row['final_rank'])}: "
            f"{row['compound_name']}"
        )

        print(
            f"PubChem CID: "
            f"{int(row['pubchem_cid'])}"
        )

        print(
            f"Docking score: "
            f"{row['docking_score']:.3f} kcal/mol"
        )

        print(
            f"QED: "
            f"{row['qed_score'] / 100:.3f}"
        )

        print(
            f"Final integrated score: "
            f"{row['final_score']:.2f}/100"
        )

    print()
    print("=" * 70)
    print("FINAL CANDIDATE RANKING COMPLETED")
    print("=" * 70)

    print()
    print("Output saved to:")
    print(OUTPUT_FILE)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()