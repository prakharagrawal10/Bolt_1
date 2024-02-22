from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors

def get_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Generate 2D coordinates for visualization
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)

        # Molecular properties
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        molecular_weight = Descriptors.MolWt(mol)
        logP = Descriptors.MolLogP(mol)
        num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        hydrogen_bond_donors = Descriptors.NumHDonors(mol)
        hydrogen_bond_acceptors = Descriptors.NumHAcceptors(mol)
        polar_surface_area = Descriptors.TPSA(mol)
        molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
        aromatic_atoms = Descriptors.NumAromaticRings(mol)
        ring_count = Descriptors.RingCount(mol)

        img = Draw.MolToImage(mol, size=(300, 300))

        # Print properties
        print(f"Number of atoms: {num_atoms}")
        print(f"Number of bonds: {num_bonds}")
        print(f"Molecular weight: {molecular_weight:.2f} g/mol")
        print(f"CLogP (lipophilicity): {logP:.2f}")
        print(f"Number of rotatable bonds: {num_rotatable_bonds}")
        print(f"Hydrogen bond donors: {hydrogen_bond_donors}")
        print(f"Hydrogen bond acceptors: {hydrogen_bond_acceptors}")
        print(f"Polar surface area: {polar_surface_area:.2f} Å^2")
        print(f"Molecular formula: {molecular_formula}")
        print(f"Aromatic atoms: {aromatic_atoms}")
        print(f"Ring count: {ring_count}")
    else:
        print(f"Invalid SMILES string: {smiles}")

# Get SMILES input from the user
user_smiles = input("Enter a SMILES string: ")
get_molecular_properties(user_smiles)
