conda create -n reaction_prediction_env python=3.10 -y
conda activate reaction-prediction_env
conda install -c conda-forge rdkit pandas numpy scikit-learn requests regex os typing -y
pip install cirpy

