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
    """Find the first matching column from a list of possible names."""
    for name in possible_names:
        if name in df.columns:
            return name
    return None


def safe_float(value, decimals=2):
    """Convert value to float safely for display."""
    try:
        return f"{float(value):.{decimals}f}"
    except (ValueError, TypeError):
        return str(value)


def main():

    print()
    print("=" * 70)
    print("FINAL LEAD ANALYSIS")
    print("=" * 70)

    # ========================================================
    # 1. CHECK INPUT FILE
    # ========================================================

    if not FINAL_FILE.exists():
        print()
        print("ERROR: Final ranking file not found:")
        print(FINAL_FILE)
        return

    # ========================================================
    # 2. LOAD FINAL RANKING
    # ========================================================

    try:
        df = pd.read_csv(FINAL_FILE)
    except Exception as error:
        print()
        print("ERROR: Could not read final ranking file.")
        print(error)
        return

    print()
    print("Loaded final ranking:")
    print(FINAL_FILE)
    print()
    print(f"Number of compounds: {len(df)}")

    if df.empty:
        print()
        print("ERROR: Final ranking file is empty.")
        return

    # ========================================================
    # 3. REQUIRED COLUMNS
    # ========================================================

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
        for column in missing_columns:
            print(f"  - {column}")

        print()
        print("Available columns:")
        for column in df.columns:
            print(f"  - {column}")

        return

    # ========================================================
    # 4. CLEAN NUMERIC COLUMNS
    # ========================================================

    numeric_columns = [
        "final_rank",
        "pubchem_cid",
        "docking_score",
        "docking_score_normalized",
        "qed_score",
        "lipinski_score",
        "final_score",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )

    # Remove rows that do not have a final rank
    df = df.dropna(
        subset=["final_rank"]
    ).copy()

    if df.empty:
        print()
        print("ERROR: No valid ranked compounds found.")
        return

    # ========================================================
    # 5. SORT BY EXISTING FINAL RANKING
    # ========================================================

    df = (
        df.sort_values(
            by="final_rank",
            ascending=True
        )
        .reset_index(drop=True)
    )

    # ========================================================
    # 6. IDENTIFY ADMET COLUMNS
    #
    # ADMET endpoints are kept individually.
    # No artificial safety score is created.
    # ========================================================

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

    # ========================================================
    # 7. LOAD INTERACTION ANALYSIS
    # ========================================================

    interaction_df = None

    if INTERACTION_FILE.exists():

        try:
            interaction_df = pd.read_csv(
                INTERACTION_FILE
            )

            print()
            print("Interaction analysis loaded:")
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

    # ========================================================
    # 8. CREATE LEAD ANALYSIS TABLE
    # ========================================================

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

    # ========================================================
    # 9. ADD ADMET INFORMATION
    # ========================================================

    for column in admet_columns:
        lead_df[column] = df[column].values

    # ========================================================
    # 10. ADD INTERACTION INFORMATION
    # ========================================================

    if interaction_df is not None and not interaction_df.empty:

        interaction_cid = find_column(
            interaction_df,
            [
                "pubchem_cid",
                "cid",
                "PubChem_CID",
                "PubChem CID",
            ],
        )

        if interaction_cid is not None:

            interaction_df = interaction_df.copy()

            # Convert both CIDs to strings so that
            # 25127713 and "25127713" match correctly.
            interaction_df[interaction_cid] = (
                interaction_df[interaction_cid]
                .astype(str)
                .str.strip()
            )

            lead_df["pubchem_cid"] = (
                lead_df["pubchem_cid"]
                .astype("Int64")
                .astype(str)
            )

            interaction_columns = [
                column
                for column in interaction_df.columns
                if column != interaction_cid
            ]

            if interaction_columns:

                interaction_subset = interaction_df[
                    [interaction_cid] + interaction_columns
                ].copy()

                interaction_subset = (
                    interaction_subset
                    .rename(
                        columns={
                            interaction_cid: "pubchem_cid"
                        }
                    )
                )

                # Avoid duplicate rows if the interaction file
                # contains multiple records for one CID.
                interaction_subset = (
                    interaction_subset
                    .drop_duplicates(
                        subset=["pubchem_cid"]
                    )
                )

                lead_df = lead_df.merge(
                    interaction_subset,
                    on="pubchem_cid",
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
                    "WARNING: Interaction file contains no "
                    "additional columns."
                )

        else:

            print()
            print(
                "WARNING: Interaction file has no PubChem CID column."
            )
            print(
                "Interaction data was not merged."
            )

    # ========================================================
    # 11. SAVE CSV
    # ========================================================

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    try:

        lead_df.to_csv(
            OUTPUT_CSV,
            index=False
        )

    except Exception as error:

        print()
        print("ERROR: Could not save lead analysis CSV.")
        print(error)
        return

    # ========================================================
    # 12. CREATE HUMAN-READABLE REPORT
    # ========================================================

    top_n = min(5, len(df))
    top_df = df.head(top_n)

    try:

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
                    f"Rank: {safe_float(row['final_rank'], 0)}\n"
                )

                report.write(
                    f"Compound: {row['compound_name']}\n"
                )

                report.write(
                    f"PubChem CID: {row['pubchem_cid']}\n"
                )

                report.write(
                    f"Docking score: "
                    f"{safe_float(row['docking_score'], 3)} "
                    f"kcal/mol\n"
                )

                report.write(
                    f"Docking normalized: "
                    f"{safe_float(row['docking_score_normalized'], 2)}\n"
                )

                report.write(
                    f"QED score: "
                    f"{safe_float(row['qed_score'], 2)}\n"
                )

                report.write(
                    f"Lipinski score: "
                    f"{safe_float(row['lipinski_score'], 2)}\n"
                )

                report.write(
                    f"Final integrated score: "
                    f"{safe_float(row['final_score'], 2)}\n"
                )

                # --------------------------------------------
                # ADMET
                # --------------------------------------------

                if admet_columns:

                    report.write(
                        "\nADMET predictions:\n"
                    )

                    for column in admet_columns:

                        value = row[column]

                        report.write(
                            f"  {column}: {value}\n"
                        )

                report.write(
                    "\n"
                )

                report.write(
                    "-" * 70 + "\n\n"
                )

    except Exception as error:

        print()
        print("ERROR: Could not create report.")
        print(error)
        return

    # ========================================================
    # 13. DISPLAY FINAL SUMMARY
    # ========================================================

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
            f"  PubChem CID: "
            f"{row['pubchem_cid']}"
        )

        print(
            f"  Docking: "
            f"{safe_float(row['docking_score'], 3)} kcal/mol"
        )

        print(
            f"  QED: "
            f"{safe_float(row['qed_score'], 2)}"
        )

        print(
            f"  Lipinski: "
            f"{safe_float(row['lipinski_score'], 2)}"
        )

        print(
            f"  Final score: "
            f"{safe_float(row['final_score'], 2)}"
        )

    # ========================================================
    # 14. OUTPUT PATHS
    # ========================================================

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

    print()
    print("Files created successfully.")


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()