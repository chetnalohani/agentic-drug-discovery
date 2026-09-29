from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

PROJECT_DIR = Path(__file__).resolve().parent.parent

INPUT_FILE = PROJECT_DIR / "results" / "final" / "lead_analysis.csv"

OUTPUT_DIR = PROJECT_DIR / "results" / "final" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD DATA
# ============================================================

def load_data():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE)

    required_columns = [
        "compound_name",
        "docking_score",
        "docking_score_normalized",
        "qed_score",
        "lipinski_score",
        "final_score",
    ]

    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    return df


# ============================================================
# 1. DOCKING SCORE COMPARISON
# ============================================================

def plot_docking_scores(df):

    data = df.sort_values(
        "docking_score",
        ascending=True
    )

    plt.figure(figsize=(11, 6))

    plt.barh(
        data["compound_name"],
        data["docking_score"]
    )

    plt.xlabel("Docking Score (kcal/mol)")
    plt.ylabel("Compound")
    plt.title("EGFR Molecular Docking Scores")

    plt.tight_layout()

    output = OUTPUT_DIR / "01_docking_scores.png"
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 2. FINAL INTEGRATED SCORE
# ============================================================

def plot_final_scores(df):

    data = df.sort_values(
        "final_score",
        ascending=False
    )

    plt.figure(figsize=(11, 6))

    plt.bar(
        data["compound_name"],
        data["final_score"]
    )

    plt.xlabel("Compound")
    plt.ylabel("Final Integrated Score")
    plt.title("Final EGFR Candidate Ranking")

    plt.xticks(
        rotation=45,
        ha="right"
    )

    plt.tight_layout()

    output = OUTPUT_DIR / "02_final_integrated_scores.png"
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 3. QED VS FINAL SCORE
# ============================================================

def plot_qed_vs_final(df):

    plt.figure(figsize=(9, 6))

    plt.scatter(
        df["qed_score"],
        df["final_score"],
        s=80
    )

    for _, row in df.iterrows():

        plt.annotate(
            row["compound_name"],
            (
                row["qed_score"],
                row["final_score"]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8
        )

    plt.xlabel("QED Score")
    plt.ylabel("Final Integrated Score")
    plt.title("QED Score vs Final Integrated Score")

    plt.tight_layout()

    output = OUTPUT_DIR / "03_qed_vs_final_score.png"
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 4. TOP 5 CANDIDATES
# ============================================================

def plot_top_candidates(df):

    data = (
        df.sort_values(
            "final_score",
            ascending=False
        )
        .head(5)
        .sort_values(
            "final_score",
            ascending=True
        )
    )

    plt.figure(figsize=(9, 5))

    plt.barh(
        data["compound_name"],
        data["final_score"]
    )

    plt.xlabel("Final Integrated Score")
    plt.ylabel("Compound")
    plt.title("Top 5 EGFR Candidate Compounds")

    plt.tight_layout()

    output = OUTPUT_DIR / "04_top_5_candidates.png"
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 60)
    print("FINAL VISUALIZATION")
    print("=" * 60)

    df = load_data()

    print()
    print(f"Compounds loaded: {len(df)}")

    plot_docking_scores(df)
    plot_final_scores(df)
    plot_qed_vs_final(df)
    plot_top_candidates(df)

    print()
    print("=" * 60)
    print("FINAL VISUALIZATION COMPLETED")
    print("=" * 60)

    print()
    print("Figures saved to:")
    print(OUTPUT_DIR)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()