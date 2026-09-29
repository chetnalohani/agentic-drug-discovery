import requests


def target_validation(target):
    """
    Target Validation Agent

    Retrieves basic biological information about a
    drug-discovery target from UniProt.
    """

    print()
    print("========================================")
    print("TARGET VALIDATION AGENT")
    print("========================================")
    print(f"Target: {target}")

    # UniProt REST API
    url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": f"gene_exact:{target} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 10
    }

    try:
        response = requests.get(
            url,
            params=params,
            timeout=30
        )

        response.raise_for_status()

        data = response.json()

    except requests.exceptions.RequestException as e:
        print()
        print("ERROR: Unable to connect to UniProt.")
        print(f"Details: {e}")
        return {}

    except ValueError:
        print()
        print("ERROR: UniProt returned invalid JSON data.")
        return {}

    results = data.get("results", [])

    if not results:
        print()
        print("No UniProt target information found.")
        return {}

    # Select the first matching reviewed human protein
    result = results[0]

    accession = result.get(
        "primaryAccession",
        ""
    )

    protein_description = (
        result.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )

    # Fallback if recommended name is unavailable
    if not protein_description:
        protein_description = (
            result.get("proteinDescription", {})
            .get("submissionNames", [{}])[0]
            .get("fullName", {})
            .get("value", "")
        )

    organism = (
        result.get("organism", {})
        .get("scientificName", "")
    )

    target_data = {
        "target": target,
        "uniprot_accession": accession,
        "protein_name": protein_description,
        "organism": organism
    }

    print()
    print("Target information retrieved successfully.")
    print(f"UniProt accession: {accession}")
    print(f"Protein: {protein_description}")
    print(f"Organism: {organism}")

    return target_data