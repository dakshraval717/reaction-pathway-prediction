import pandas
import numpy

import sklearn

from rdkit import Chem
from rdkit.Chem import rdmolops # for graph generation and operations/manipulations on molecules

smiles = 'CCO'  # Example SMILES string for ethanol
molecule = Chem.MolFromSmiles(smiles)

# nodes of graph neural network are atoms in the molecule
atoms = []
for atom in molecule.GetAtoms():
    atoms.append({
        'Atom Index': atom.GetIdx(),
        'Atomic Number': atom.GetAtomicNum(),
        'Periodic Table Symbol': atom.GetSymbol(),
        'Bonding Degree': atom.GetDegree(), # number of bonds to neighbours 
        'Formal Charge': atom.GetFormalCharge(), # Lewis structure formal charge
        'Hybridization': atom.GetHybridization().name, # valence bond hybridization state (sp, sp2, tetrahedral sp3, planar sp3, etc.)
        'Aromatic?': atom.GetIsAromatic(), # whether the atom is part of an aromatic ring
        'In Benzene Ring?': atom.GetIsAromatic(), # whether the atom is part of a benzene ring
        'In Ring?': atom.IsInRing(), # whether the atom is part of any ring
    }
    )

# edges of graph neural network are bonds in the molecule
bonds = []
for bond in molecule.GetBonds():
    bonds.append({
        '1st Atom Index': bond.GetBeginAtomIdx(),
        '2nd Atom Index': bond.GetEndAtomIdx(),
        'Bond Type': bond.GetBondType().name, # type of bond (single, double, triple, aromatic)
        'Aromatic?': bond.GetIsAromatic(), 
        'Conjugated?': bond.GetIsConjugated(), # part of chain of double/triple-single-double/triple bonds
        'In Ring?': bond.IsInRing()
    })


# specifies which nodes (atoms) are connected; chemical with N atoms has N x N adjacency matrix
adj_matrix = rdmolops.GetAdjacencyMatrix(molecule)
print("Adjacency Matrix:")
print(adj_matrix)

print(atoms)
print(bonds)

