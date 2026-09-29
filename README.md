# Agentic Drug Discovery System

An agentic computational drug-discovery pipeline integrating target validation, literature research, compound discovery, drug-likeness filtering, molecular docking, interaction analysis, ADMET prediction, candidate ranking, lead analysis, and visualization.

---

## Target

- **Target:** EGFR (Epidermal Growth Factor Receptor)
- **PDB:** 1M17
- **Docking Software:** AutoDock Vina v1.2.7

---

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
Candidate Ranking  
↓  
ADMET Prediction  
↓  
Final Integrated Ranking  
↓  
Lead Analysis  
↓  
Visualization

---

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

---

## Project Structure

```text
agentic-drug-discovery/
│
├── tools/
│   ├── __init__.py
│   ├── admet_analysis.py
│   ├── final_candidate_ranking.py
│   ├── final_visualization.py
│   ├── interaction_analysis.py
│   └── lead_analysis.py
│
├── results/
│   ├── docking/
│   │   └── interaction_analysis/
│   │
│   └── final/
│       ├── figures/
│       │   ├── 01_docking_scores.png
│       │   ├── 02_final_integrated_scores.png
│       │   ├── 03_qed_vs_final_score.png
│       │   └── 04_top_5_candidates.png
│       │
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
├── .env
└── .gitignore