from pathlib import Path
import pandas as pd


# ============================================================
# LEAD ANALYSIS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

FINAL_FILE = (
    PROJECT_ROOT
    / "results"
    / "final"
    / "final_candidate_ranking.csv"
)

INTERACTION_FILE = (
    PROJECT_ROOT
    / "results"
    / "docking"
    / "interaction_analysis"
    / "interaction_summary.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "final"
)

OUTPUT_CSV = OUTPUT_DIR / "lead_analysis.csv"
OUTPUT_TXT = OUTPUT_DIR / "lead_analysis_report.txt"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def find_column(df, possible_names):
    """
    Find the first matching column from a list of possible names.
    """

    for name in possible_names:
        if name in df.columns:
            return name

    return None


def main():

    print()
    print("=" * 70)
    print("FINAL LEAD ANALYSIS")
    print("=" * 70)

    # --------------------------------------------------------
    # Check input file
    # --------------------------------------------------------

    if not FINAL_FILE.exists():

        print()
        print("ERROR: Final ranking file not found.")
        print(FINAL_FILE)
        return

    # --------------------------------------------------------
    # Load final ranking
    # --------------------------------------------------------

    df = pd.read_csv(FINAL_FILE)

    print()
    print(f"Loaded final ranking:")
    print(FINAL_FILE)

    print()
    print(f"Number of compounds: {len(df)}")

    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    required_columns = [
        "final_rank",
        "pubchem_cid",
        "compound_name",
        "docking_score",
        "docking_score_normalized",
        "qed_score",
        "lipinski_score",
        "final_score",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:

        print()
        print("ERROR: Required columns are missing:")
        print(missing_columns)
        return

    # --------------------------------------------------------
    # Sort according to existing final ranking
    # --------------------------------------------------------

    df = df.sort_values(
        by="final_rank",
        ascending=True
    ).reset_index(drop=True)

    # --------------------------------------------------------
    # Identify ADMET columns
    #
    # ADMET endpoints are retained individually.
    # They are NOT converted into one artificial safety score.
    # --------------------------------------------------------

    base_columns = {
        "final_rank",
        "pubchem_cid",
        "compound_name",
        "smiles",
        "molecular_weight",
        "logP",
        "hydrogen_bond_acceptors",
        "hydrogen_bond_donors",
        "Lipinski",
        "QED",
        "stereo_centers",
        "tpsa",
        "PAINS_alert",
        "docking_score",
        "docking_score_normalized",
        "qed_score",
        "lipinski_score",
        "final_score",
    }

    admet_columns = [
        column
        for column in df.columns
        if column not in base_columns
    ]

    # --------------------------------------------------------
    # Load interaction analysis if available
    # --------------------------------------------------------

    interaction_df = None

    if INTERACTION_FILE.exists():

        try:

            interaction_df = pd.read_csv(
                INTERACTION_FILE
            )

            print()
            print(
                "Interaction analysis loaded:"
            )
            print(INTERACTION_FILE)

        except Exception as error:

            print()
            print(
                "WARNING: Could not load interaction analysis."
            )
            print(error)

    else:

        print()
        print(
            "WARNING: Interaction summary not found."
        )

    # --------------------------------------------------------
    # Create lead-analysis table
    # --------------------------------------------------------

    analysis_columns = [
        "final_rank",
        "pubchem_cid",
        "compound_name",
        "docking_score",
        "docking_score_normalized",
        "qed_score",
        "lipinski_score",
        "final_score",
    ]

    lead_df = df[analysis_columns].copy()

    # --------------------------------------------------------
    # Add ADMET information
    # --------------------------------------------------------

    for column in admet_columns:

        lead_df[column] = df[column]

    # --------------------------------------------------------
    # Add interaction information when possible
    # --------------------------------------------------------

    if interaction_df is not None:

        interaction_cid = find_column(
            interaction_df,
            [
                "pubchem_cid",
                "cid",
                "PubChem_CID",
                "PubChem CID",
            ],
        )

        main_cid = "pubchem_cid"

        if interaction_cid is not None:

            interaction_df = interaction_df.copy()

            interaction_df[interaction_cid] = (
                interaction_df[interaction_cid]
                .astype(str)
            )

            lead_df[main_cid] = (
                lead_df[main_cid]
                .astype(str)
            )

            interaction_columns = [
                column
                for column in interaction_df.columns
                if column != interaction_cid
            ]

            interaction_subset = interaction_df[
                [interaction_cid] + interaction_columns
            ].copy()

            interaction_subset = interaction_subset.rename(
                columns={
                    interaction_cid: main_cid
                }
            )

            lead_df = lead_df.merge(
                interaction_subset,
                on=main_cid,
                how="left",
                suffixes=("", "_interaction"),
            )

            print()
            print(
                "Interaction information merged successfully."
            )

        else:

            print()
            print(
                "Interaction file has no PubChem CID column."
            )
            print(
                "Interaction data will remain separate."
            )

    # --------------------------------------------------------
    # Save CSV
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    lead_df.to_csv(
        OUTPUT_CSV,
        index=False
    )

    # --------------------------------------------------------
    # Create human-readable report
    # --------------------------------------------------------

    top_n = min(5, len(df))

    top_df = df.head(top_n)

    with open(
        OUTPUT_TXT,
        "w",
        encoding="utf-8"
    ) as report:

        report.write(
            "FINAL LEAD ANALYSIS REPORT\n"
        )

        report.write(
            "=" * 70 + "\n\n"
        )

        report.write(
            "Target: EGFR\n"
        )

        report.write(
            f"Total compounds analyzed: {len(df)}\n\n"
        )

        report.write(
            "IMPORTANT:\n"
        )

        report.write(
            "The final integrated score was retained from the "
            "existing ranking workflow.\n"
        )

        report.write(
            "ADMET endpoints are reported individually and are "
            "not collapsed into an artificial safety score.\n\n"
        )

        report.write(
            "TOP LEAD CANDIDATES\n"
        )

        report.write(
            "=" * 70 + "\n\n"
        )

        for _, row in top_df.iterrows():

            report.write(
                f"Rank: {row['final_rank']}\n"
            )

            report.write(
                f"Compound: {row['compound_name']}\n"
            )

            report.write(
                f"PubChem CID: {row['pubchem_cid']}\n"
            )

            report.write(
                f"Docking score: "
                f"{row['docking_score']:.3f} kcal/mol\n"
            )

            report.write(
                f"Docking normalized: "
                f"{row['docking_score_normalized']:.2f}\n"
            )

            report.write(
                f"QED score: "
                f"{row['qed_score']:.2f}\n"
            )

            report.write(
                f"Lipinski score: "
                f"{row['lipinski_score']:.2f}\n"
            )

            report.write(
                f"Final integrated score: "
                f"{row['final_score']:.2f}\n"
            )

            report.write("\n")

            if admet_columns:

                report.write(
                    "ADMET predictions:\n"
                )

                for column in admet_columns:

                    value = row[column]

                    report.write(
                        f"  {column}: {value}\n"
                    )

                report.write("\n")

            report.write(
                "-" * 70 + "\n\n"
            )

    # --------------------------------------------------------
    # Display final summary
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("TOP LEAD CANDIDATES")
    print("=" * 70)

    for _, row in top_df.iterrows():

        print()
        print(
            f"Rank {int(row['final_rank'])}: "
            f"{row['compound_name']}"
        )

        print(
            f"  PubChem CID: {row['pubchem_cid']}"
        )

        print(
            f"  Docking: "
            f"{row['docking_score']:.3f} kcal/mol"
        )

        print(
            f"  QED: "
            f"{row['qed_score']:.2f}"
        )

        print(
            f"  Lipinski: "
            f"{row['lipinski_score']:.2f}"
        )

        print(
            f"  Final score: "
            f"{row['final_score']:.2f}"
        )

    # --------------------------------------------------------
    # Output paths
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("LEAD ANALYSIS COMPLETED")
    print("=" * 70)

    print()
    print("CSV saved to:")
    print(OUTPUT_CSV)

    print()
    print("Report saved to:")
    print(OUTPUT_TXT)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()