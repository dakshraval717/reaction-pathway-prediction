# Reaction Pathway Prediction

A machine learning pipeline that predicts organic reaction mechanisms by analyzing chemical structure changes. This project automates the extraction of reaction data from NIST and classifies reaction pathways using Morgan fingerprinting and Random Forest classifiers.

## Project Structure

```bash
reaction-pathway-prediction/
├── data/       # Dataset storage (Git LFS tracked)
├── scripts/    # Data processing and ML pipelines
├── docs/       # Documentation and visuals
└── requirements.txt
```

## Key Components

### 1. Data Acquisition (`pandas_scraper.py`)
- Automated scraper for the NIST Chemical Kinetics Database.
- Extracts activation energies, pre-exponential factors, and reaction orders using regex.
- Converts raw HTML reaction strings into structured CSV format.

### 2. Feature Engineering (`morgan_processing.py`)
- Converts chemical reaction SMILES into 4096-bit Morgan difference fingerprints.
- Captures structural changes (bond breaking/forming) as vector inputs for ML models.

### 3. Mechanism Classification (`mechanism_classifier.py`)
- **Model:** Random Forest Classifier (scikit-learn).
- **Task:** Predicts reaction mechanism steps based on structural fingerprint differences.
- **Input:** 4096-bit difference vectors.

## Installation & Usage

1. **Environment Setup:**
   ```bash
   ./setup_conda_env.sh
   conda activate thermo-ml
   ```

2. **Run the Pipeline:**
   ```bash
   cd scripts
   python morgan_processing.py  # Generate features
   python mechanism_classifier.py  # Train model
   ```

## Dependencies
- **RDKit**: Chemical informatics and fingerprint generation.
- **scikit-learn**: Model training and evaluation.
- **pandas/numpy**: Data manipulation.
