import requests
import time


PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def compound_discovery(target, max_compounds=10):
    """
    Compound Discovery Agent

    Searches PubChem for compounds associated with the target
    and retrieves basic chemical properties.
    """

    print()
    print("=" * 60)
    print("COMPOUND DISCOVERY AGENT")
    print("=" * 60)
    print(f"Target: {target}")

    # Known EGFR compounds for the first prototype.
    target_compounds = {
        "EGFR": [
            "gefitinib",
            "erlotinib",
            "afatinib",
            "osimertinib",
            "lapatinib",
            "vandetanib",
            "canertinib",
            "icotinib",
            "pelitinib",
            "poziotinib",
        ]
    }

    compounds_to_search = target_compounds.get(
        target.upper(),
        []
    )

    if not compounds_to_search:
        print(
            f"No compound search list configured for target: {target}"
        )
        return []

    compounds = []

    for compound_name in compounds_to_search[:max_compounds]:

        try:
            # PubChem properties required by the pipeline.
            properties = (
                "MolecularFormula,"
                "MolecularWeight,"
                "ConnectivitySMILES,"
                "SMILES,"
                "InChIKey,"
                "IUPACName,"
                "XLogP,"
                "TPSA,"
                "HBondDonorCount,"
                "HBondAcceptorCount,"
                "RotatableBondCount"
            )

            url = (
                f"{PUBCHEM_URL}/compound/name/"
                f"{compound_name}/property/"
                f"{properties}/JSON"
            )

            response = requests.get(
                url,
                timeout=30
            )

            response.raise_for_status()

            data = response.json()

            property_list = (
                data
                .get("PropertyTable", {})
                .get("Properties", [])
            )

            if not property_list:
                print(
                    f"Could not retrieve data for "
                    f"{compound_name}"
                )
                continue

            pubchem_data = property_list[0]

            # IMPORTANT:
            # Keep the original compound identity fields.
            # These must remain available for later agents.
            compound_data = {
                "target": target.upper(),

                # Compound identity
                "compound_name": compound_name,
                "pubchem_cid": pubchem_data.get("CID", ""),

                # Chemical properties
                "molecular_formula": pubchem_data.get(
                    "MolecularFormula",
                    ""
                ),

                "molecular_weight": pubchem_data.get(
                    "MolecularWeight",
                    ""
                ),

                "canonical_smiles": pubchem_data.get(
                    "ConnectivitySMILES",
                    ""
                ),

                "isomeric_smiles": pubchem_data.get(
                    "SMILES",
                    ""
                ),

                "inchikey": pubchem_data.get(
                    "InChIKey",
                    ""
                ),

                "iupac_name": pubchem_data.get(
                    "IUPACName",
                    ""
                ),

                "xlogp": pubchem_data.get(
                    "XLogP",
                    ""
                ),

                "tpsa": pubchem_data.get(
                    "TPSA",
                    ""
                ),

                "h_bond_donors": pubchem_data.get(
                    "HBondDonorCount",
                    ""
                ),

                "h_bond_acceptors": pubchem_data.get(
                    "HBondAcceptorCount",
                    ""
                ),

                "rotatable_bonds": pubchem_data.get(
                    "RotatableBondCount",
                    ""
                ),

                # Provenance
                "source": "PubChem"
            }

            compounds.append(compound_data)

            # Display retrieved compound
            print()
            print(f"Compound: {compound_data['compound_name']}")
            print(
                f"PubChem CID: "
                f"{compound_data['pubchem_cid']}"
            )
            print(
                f"Molecular weight: "
                f"{compound_data['molecular_weight']}"
            )
            print(
                f"Molecular formula: "
                f"{compound_data['molecular_formula']}"
            )
            print(
                f"XLogP: "
                f"{compound_data['xlogp']}"
            )
            print(
                f"HBD: "
                f"{compound_data['h_bond_donors']}"
            )
            print(
                f"HBA: "
                f"{compound_data['h_bond_acceptors']}"
            )

            # Stay below PubChem request rate.
            time.sleep(0.25)

        except requests.RequestException as error:

            print(
                f"PubChem request failed for "
                f"{compound_name}: {error}"
            )

        except (
            ValueError,
            KeyError,
            TypeError
        ) as error:

            print(
                f"Data processing failed for "
                f"{compound_name}: {error}"
            )

    print()
    print("=" * 60)
    print(
        f"Compound discovery completed. "
        f"Compounds retrieved: {len(compounds)}"
    )
    print("=" * 60)

    return compounds


# ============================================================
# TEST THE COMPOUND DISCOVERY AGENT DIRECTLY
# ============================================================

if __name__ == "__main__":

    results = compound_discovery(
        "EGFR",
        max_compounds=10
    )

    print()
    print("=" * 60)
    print("FINAL COMPOUND RESULTS")
    print("=" * 60)

    if not results:

        print("No compounds were retrieved.")

    else:

        for index, compound in enumerate(
            results,
            start=1
        ):

            print()
            print(f"Compound {index}")

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

    print()
    print("=" * 60)
    print("COMPOUND DISCOVERY STAGE COMPLETED")
    print("=" * 60)