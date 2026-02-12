# Predicting Chemical Reaction Pathways

### Language Model Method

Trained a language model to decompose 'black-box' organic chemical reactions into likely pathways/steps by encoding 3D geometric structure, bond resonance, and aromatic rings using unique strings.

![SMILES Chemical String](SMILES.webp)

### Bit Vector Method

Trained a Random Forest model with scikit-learn to classify single-step reaction types by encoding 'reaction changes' as a 4096-bit vector.

![Morgan Fingerprint Bit Vector](morgan.png)
