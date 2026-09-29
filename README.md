# Agentic Drug Discovery System

An agentic computational drug-discovery pipeline designed to automate target validation, literature research, compound discovery, drug-likeness filtering, molecular docking, interaction analysis, ADMET prediction, candidate ranking, lead analysis, and visualization.

## Project Overview

This project integrates multiple computational drug-discovery stages into a single automated workflow.

The pipeline uses EGFR as the target protein and evaluates candidate compounds through sequential computational analysis.

> **Important:** This is a computational drug-discovery project. Docking, ADMET, and ranking results are predictions and require experimental validation before any biological or therapeutic conclusions can be made.

## Target

- **Target:** EGFR (Epidermal Growth Factor Receptor)
- **PDB:** 1M17
- **Docking Software:** AutoDock Vina v1.2.7

## Workflow

The complete computational workflow is:

Target Validation
↓
Literature Research
↓
Compound Discovery
↓
Lipinski Filtering
↓
Molecular Docking
↓
Interaction Analysis
↓
ADMET Prediction
↓
Candidate Ranking
↓
Final Integrated Ranking
↓
Lead Analysis
↓
Visualization

## Key Features

- Automated target validation
- Literature research for the selected target
- Compound discovery using PubChem
- Drug-likeness filtering using Lipinski's Rule of Five
- Molecular docking using AutoDock Vina
- Protein-ligand interaction analysis
- Computational ADMET prediction
- Integrated candidate scoring
- Lead candidate analysis
- Automated visualization and result generation

## Technologies

- Python
- RDKit
- Meeko
- ProDy
- AutoDock Vina
- ADMET-AI
- Pandas
- NumPy
- Matplotlib
- PubChem
- RCSB PDB
- PubMed / NCBI
- UCSF ChimeraX
- Git / GitHub

## Project Structure

```text
agentic-drug-discovery/
│
├── agents/
│   ├── target_validation_agent.py
│   ├── literature_agent.py
│   ├── compound_discovery_agent.py
│   └── compound_filter_agent.py
│
├── tools/
│   ├── admet_analysis.py
│   ├── final_candidate_ranking.py
│   ├── final_visualization.py
│   ├── interaction_analysis.py
│   └── lead_analysis.py
│
├── results/
│   ├── docking/
│   ├── final/
│   │   ├── figures/
│   │   ├── final_candidate_ranking.csv
│   │   ├── lead_analysis.csv
│   │   └── lead_analysis_report.txt
│   └── egfr_literature.json
│
├── vina/
│
├── 1M17.pdb
├── main.py
├── requirements.txt
├── README.md
└── .gitignore
