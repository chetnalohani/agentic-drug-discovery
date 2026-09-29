from agents.compound_discovery_agent import compound_discovery
from agents.compound_filter_agent import compound_filter


def calculate_candidate_score(compound):
    """
    Rank drug-like compounds using basic chemical properties.

    Higher score = higher priority for further evaluation.

    This is a prioritization score, NOT a docking score.
    """

    score = 0.0
    reasons = []

    # --------------------------------------------------------
    # Extract properties safely
    # --------------------------------------------------------

    try:
        molecular_weight = float(
            compound.get("molecular_weight", 0)
        )
    except (ValueError, TypeError):
        molecular_weight = 0.0

    try:
        xlogp = float(
            compound.get("xlogp", 0)
        )
    except (ValueError, TypeError):
        xlogp = 0.0

    try:
        hbd = int(
            compound.get("h_bond_donors", 0)
        )
    except (ValueError, TypeError):
        hbd = 0

    try:
        hba = int(
            compound.get("h_bond_acceptors", 0)
        )
    except (ValueError, TypeError):
        hba = 0

    try:
        tpsa = float(
            compound.get("tpsa", 0)
        )
    except (ValueError, TypeError):
        tpsa = 0.0

    # --------------------------------------------------------
    # 1. Lipinski compliance
    # --------------------------------------------------------

    if compound.get("lipinski_pass", False):
        score += 30
        reasons.append("Lipinski compliant")

    # --------------------------------------------------------
    # 2. Molecular weight
    # --------------------------------------------------------

    if 250 <= molecular_weight <= 450:
        score += 20
        reasons.append("favorable molecular weight")

    elif 450 < molecular_weight <= 500:
        score += 10
        reasons.append("acceptable molecular weight")

    # --------------------------------------------------------
    # 3. XLogP
    # --------------------------------------------------------

    if 1 <= xlogp <= 4:
        score += 20
        reasons.append("favorable lipophilicity")

    elif 4 < xlogp <= 5:
        score += 10
        reasons.append("acceptable lipophilicity")

    # --------------------------------------------------------
    # 4. Hydrogen bond donors
    # --------------------------------------------------------

    if 0 <= hbd <= 3:
        score += 10

    # --------------------------------------------------------
    # 5. Hydrogen bond acceptors
    # --------------------------------------------------------

    if 2 <= hba <= 8:
        score += 10

    # --------------------------------------------------------
    # 6. TPSA
    # --------------------------------------------------------

    if 40 <= tpsa <= 120:
        score += 10
        reasons.append("favorable TPSA")

    # --------------------------------------------------------
    # Return ranking information
    # --------------------------------------------------------

    return {
        "candidate_score": round(score, 2),
        "ranking_reasons": reasons,
        "molecular_weight": molecular_weight,
        "xlogp": xlogp,
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "tpsa": tpsa,
    }


def rank_candidates(compounds):
    """
    Rank filtered compounds according to their
    calculated candidate-prioritization score.
    """

    print()
    print("=" * 60)
    print("CANDIDATE RANKING AGENT")
    print("=" * 60)

    if not compounds:
        print("No compounds available for ranking.")
        return []

    ranked_compounds = []

    for compound in compounds:

        evaluation = calculate_candidate_score(compound)

        ranked_compound = compound.copy()
        ranked_compound.update(evaluation)

        ranked_compounds.append(ranked_compound)

    # Highest score first
    ranked_compounds.sort(
        key=lambda x: x.get("candidate_score", 0),
        reverse=True
    )

    print()
    print("CANDIDATE RANKING RESULTS")
    print("=" * 60)

    for index, compound in enumerate(
        ranked_compounds,
        start=1
    ):

        print()
        print(f"Rank {index}")
        print(
            f"Name: "
            f"{compound.get('compound_name', '')}"
        )
        print(
            f"PubChem CID: "
            f"{compound.get('pubchem_cid', '')}"
        )
        print(
            f"Candidate Score: "
            f"{compound.get('candidate_score', 0)}"
        )
        print(
            f"MW: "
            f"{compound.get('molecular_weight', '')}"
        )
        print(
            f"XLogP: "
            f"{compound.get('xlogp', '')}"
        )
        print(
            f"HBD: "
            f"{compound.get('h_bond_donors', '')}"
        )
        print(
            f"HBA: "
            f"{compound.get('h_bond_acceptors', '')}"
        )
        print(
            f"TPSA: "
            f"{compound.get('tpsa', '')}"
        )

    print()
    print("=" * 60)
    print("CANDIDATE RANKING STAGE COMPLETED")
    print("=" * 60)

    return ranked_compounds


# ============================================================
# TEST STAGE 5
# ============================================================

if __name__ == "__main__":

    print()
    print("=" * 60)
    print("STAGE 5: CANDIDATE RANKING")
    print("=" * 60)

    # --------------------------------------------------------
    # Stage 3: Compound discovery
    # --------------------------------------------------------

    compounds = compound_discovery(
        "EGFR",
        max_compounds=10
    )

    # --------------------------------------------------------
    # Stage 4: Lipinski filtering
    # --------------------------------------------------------

    filtered_compounds = compound_filter(
        compounds
    )

    # --------------------------------------------------------
    # Stage 5: Candidate ranking
    # --------------------------------------------------------

    ranked_compounds = rank_candidates(
        filtered_compounds
    )

    # --------------------------------------------------------
    # Final candidates
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("TOP EGFR CANDIDATES")
    print("=" * 60)

    if not ranked_compounds:

        print("No candidates available.")

    else:

        for index, compound in enumerate(
            ranked_compounds,
            start=1
        ):

            print()
            print(f"Rank {index}")
            print(
                f"Compound: "
                f"{compound.get('compound_name', '')}"
            )
            print(
                f"PubChem CID: "
                f"{compound.get('pubchem_cid', '')}"
            )
            print(
                f"Candidate Score: "
                f"{compound.get('candidate_score', 0)}"
            )

    print()
    print("=" * 60)
    print("STAGE 5 COMPLETED")
    print("=" * 60)