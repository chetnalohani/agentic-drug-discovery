from agents.compound_discovery_agent import compound_discovery


def check_lipinski(compound):
    """
    Evaluate a compound using Lipinski Rule of Five.

    Rules:
    - Molecular weight <= 500 Da
    - XLogP <= 5
    - H-bond donors <= 5
    - H-bond acceptors <= 10
    """

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

    violations = []

    # Lipinski Rule of Five checks
    if molecular_weight > 500:
        violations.append("MW > 500")

    if xlogp > 5:
        violations.append("XLogP > 5")

    if hbd > 5:
        violations.append("HBD > 5")

    if hba > 10:
        violations.append("HBA > 10")

    passed = len(violations) == 0

    return {
        "lipinski_pass": passed,
        "violations": violations,
        "molecular_weight": molecular_weight,
        "xlogp": xlogp,
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba
    }


def compound_filter(compounds):
    """
    Filter discovered compounds according to
    Lipinski drug-likeness criteria.

    The original compound information is preserved.
    """

    print()
    print("=" * 60)
    print("COMPOUND FILTERING AGENT")
    print("=" * 60)

    filtered_compounds = []

    for compound in compounds:

        # ------------------------------------------
        # Evaluate Lipinski properties
        # ------------------------------------------

        evaluation = check_lipinski(compound)

        # ------------------------------------------
        # IMPORTANT:
        # Preserve ALL original compound information
        # and add the Lipinski evaluation to it.
        # ------------------------------------------

        compound_result = compound.copy()

        compound_result.update(evaluation)

        # ------------------------------------------
        # PASS / FAIL
        # ------------------------------------------

        if evaluation["lipinski_pass"]:

            compound_result["drug_like"] = True

            filtered_compounds.append(
                compound_result
            )

            print()
            print(
                f"PASS: "
                f"{compound_result.get('compound_name', 'Unknown')}"
            )

        else:

            compound_result["drug_like"] = False

            print()
            print(
                f"FAIL: "
                f"{compound_result.get('compound_name', 'Unknown')}"
            )

            if evaluation["violations"]:

                print(
                    "Reason: "
                    + ", ".join(
                        evaluation["violations"]
                    )
                )

    # ------------------------------------------
    # FILTERING SUMMARY
    # ------------------------------------------

    print()
    print("=" * 60)
    print(
        f"Total compounds evaluated: "
        f"{len(compounds)}"
    )

    print(
        f"Drug-like compounds: "
        f"{len(filtered_compounds)}"
    )

    print("=" * 60)

    return filtered_compounds


# ============================================================
# DIRECT TEST
# ============================================================

if __name__ == "__main__":

    print()
    print("=" * 60)
    print("STAGE 4: COMPOUND FILTERING")
    print("=" * 60)

    # ------------------------------------------
    # Stage 3:
    # Discover compounds
    # ------------------------------------------

    compounds = compound_discovery(
        "EGFR",
        max_compounds=10
    )

    if not compounds:

        print()
        print("No compounds were discovered.")
        print("Filtering stage cannot continue.")

    else:

        # ------------------------------------------
        # Stage 4:
        # Apply Lipinski filtering
        # ------------------------------------------

        filtered = compound_filter(
            compounds
        )

        # ------------------------------------------
        # FINAL FILTERED RESULTS
        # ------------------------------------------

        print()
        print("=" * 60)
        print("FINAL FILTERED COMPOUNDS")
        print("=" * 60)

        if not filtered:

            print()
            print(
                "No compounds passed "
                "the Lipinski filtering stage."
            )

        else:

            for index, compound in enumerate(
                filtered,
                start=1
            ):

                print()
                print(
                    f"Compound {index}"
                )

                print(
                    f"Name: "
                    f"{compound.get('compound_name', '')}"
                )

                print(
                    f"PubChem CID: "
                    f"{compound.get('pubchem_cid', '')}"
                )

                print(
                    f"Formula: "
                    f"{compound.get('molecular_formula', '')}"
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
                    f"Lipinski Pass: "
                    f"{compound.get('lipinski_pass', False)}"
                )

                print(
                    f"Drug-like: "
                    f"{compound.get('drug_like', False)}"
                )

                violations = compound.get(
                    "violations",
                    []
                )

                if violations:

                    print(
                        "Violations: "
                        + ", ".join(violations)
                    )

                else:

                    print(
                        "Violations: None"
                    )

        print()
        print("=" * 60)
        print(
            "COMPOUND FILTERING STAGE COMPLETED"
        )
        print("=" * 60)