import requests
import xml.etree.ElementTree as ET
import json


def literature_research(target, max_results=10):

    print(f"Searching PubMed for relevant literature on: {target}")

    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # More specific search than simply searching the target name
    search_term = (
        f'"{target}"[Title/Abstract] AND '
        f'(inhibitor OR inhibition OR kinase OR mutation OR '
        f'cancer OR therapy OR drug)'
    )

    params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "json",
        "retmax": max_results
    }

    response = requests.get(
        search_url,
        params=params,
        timeout=30
    )

    response.raise_for_status()

    data = response.json()

    pmids = data["esearchresult"]["idlist"]

    print(f"Found {len(pmids)} potentially relevant papers.")

    if not pmids:
        return []

    # ---------------------------------------------------------
    # Retrieve detailed PubMed information
    # ---------------------------------------------------------

    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml"
    }

    fetch_response = requests.get(
        fetch_url,
        params=fetch_params,
        timeout=30
    )

    fetch_response.raise_for_status()

    root = ET.fromstring(fetch_response.text)

    papers = []

    for article in root.findall(".//PubmedArticle"):

        pmid_element = article.find(".//PMID")
        title_element = article.find(".//ArticleTitle")

        pmid = (
            pmid_element.text
            if pmid_element is not None
            else ""
        )

        title = (
            "".join(title_element.itertext())
            if title_element is not None
            else ""
        )

        # -----------------------------------------------------
        # Authors
        # -----------------------------------------------------

        authors = []

        for author in article.findall(".//Author"):

            name_parts = []

            lastname = author.find("LastName")
            initials = author.find("Initials")

            if lastname is not None:
                name_parts.append(lastname.text)

            if initials is not None:
                name_parts.append(initials.text)

            if name_parts:
                authors.append(" ".join(name_parts))

        # -----------------------------------------------------
        # Journal
        # -----------------------------------------------------

        journal_element = article.find(".//Journal/Title")

        journal = (
            journal_element.text
            if journal_element is not None
            else ""
        )

        # -----------------------------------------------------
        # Abstract
        # -----------------------------------------------------

        abstract_parts = []

        for abstract_text in article.findall(".//AbstractText"):

            text = "".join(abstract_text.itertext())

            if text:
                abstract_parts.append(text)

        abstract = " ".join(abstract_parts)

        # -----------------------------------------------------
        
        # -----------------------------------------------------

        combined_text = (
            title + " " + abstract
        ).lower()

        target_lower = target.lower()

        relevance_score = 0

        if target_lower in title.lower():
            relevance_score += 5

        if target_lower in abstract.lower():
            relevance_score += 3

        important_terms = [
            "inhibitor",
            "inhibition",
            "kinase",
            "mutation",
            "cancer",
            "therapy",
            "drug"
        ]

        for term in important_terms:

            if term in combined_text:
                relevance_score += 1

        # -----------------------------------------------------
        # Store structured paper information
        # -----------------------------------------------------

        paper = {
            "pmid": pmid,
            "title": title,
            "authors": authors[:10],
            "journal": journal,
            "abstract": abstract,
            "relevance_score": relevance_score
        }

        papers.append(paper)

    # ---------------------------------------------------------
    # Rank papers by relevance
    # ---------------------------------------------------------

    papers.sort(
        key=lambda paper: paper["relevance_score"],
        reverse=True
    )

    # Keep only the requested number
    papers = papers[:max_results]

    print()
    print("Literature agent completed.")
    print(f"Relevant papers retrieved: {len(papers)}")

    # ---------------------------------------------------------
    # Save results
    # ---------------------------------------------------------

    with open(
        "results/egfr_literature.json",
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            papers,
            f,
            indent=4,
            ensure_ascii=False
        )

    print()
    print("Research results saved to:")
    print("results/egfr_literature.json")

    return papers


# -------------------------------------------------------------
# Test the literature agent
# -------------------------------------------------------------

if __name__ == "__main__":

    results = literature_research("EGFR", max_results=10)

    print()
    print("TOP RELEVANT PAPERS")
    print("====================")

    for i, paper in enumerate(results, start=1):

        print()
        print(f"PAPER {i}")
        print(f"PMID: {paper['pmid']}")
        print(f"Title: {paper['title']}")
        print(f"Journal: {paper['journal']}")
        print(f"Relevance score: {paper['relevance_score']}")