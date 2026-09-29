# Agentic Drug Discovery System

An agentic computational drug-discovery pipeline integrating target validation, literature research, compound discovery, drug-likeness filtering, molecular docking, interaction analysis, ADMET prediction, candidate ranking, lead analysis, and visualization.

## Target

- Target: EGFR (Epidermal Growth Factor Receptor)
- PDB: 1M17
- Docking software: AutoDock Vina v1.2.7

## Workflow

Target Validation
→ Literature Research
→ Compound Discovery
→ Lipinski Filtering
→ Molecular Docking
→ Interaction Analysis
→ Candidate Ranking
→ ADMET Prediction
→ Final Integrated Ranking
→ Lead Analysis
→ Visualization

## Technologies

- Python 3.11
- RDKit
- Meeko
- ProDy
- Gemmi
- AutoDock Vina
- ADMET-AI
- Pandas
- NumPy
- Matplotlib
- UCSF ChimeraX
- PubChem
- RCSB PDB
- PubMed / NCBI
- Anaconda
- Visual Studio Code
- Git

## Project Structure

```text
agentic - drug - discovery/
│
├── agents/
│   ├── __init__.py
│   ├── target_validation_agent.py
│   ├── literature_agent.py
│   ├── compound_discovery_agent.py
│   ├── compound_filter_agent.py
│   ├── docking_agent.py
│   └── candidate_ranking_agent.py
│
├── tools/
│   ├── __init__.py
│   ├── admet_analysis.py
│   ├── interaction_analysis.py
│   ├── final_candidate_ranking.py
│   ├── lead_analysis.py
│   └── final_visualization.py
│
├── config/
├── data/
│
├── results/
│   ├── admet/
│   ├── docking/
│   └── final/
│       ├── figures/
│       ├── final_candidate_ranking.csv
│       ├── lead_analysis.csv
│       ├── lead_analysis_report.txt
│       └── egfr_literature.json
│
├── vina/
│   └── vina.exe
│
├── 1M17.pdb
├── main.py
├── requirements.txt
├── README.md
└── .env