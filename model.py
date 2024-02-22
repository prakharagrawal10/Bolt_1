from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors

def get_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    result = {}

    if mol is not None:
        # Generate 2D coordinates for visualization
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)

        # Molecular properties
        result['num_atoms'] = mol.GetNumAtoms()
        result['num_bonds'] = mol.GetNumBonds()
        result['molecular_weight'] = Descriptors.MolWt(mol)
        result['logP'] = Descriptors.MolLogP(mol)
        result['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
        result['hydrogen_bond_donors'] = Descriptors.NumHDonors(mol)
        result['hydrogen_bond_acceptors'] = Descriptors.NumHAcceptors(mol)
        result['polar_surface_area'] = Descriptors.TPSA(mol)
        result['molecular_formula'] = rdMolDescriptors.CalcMolFormula(mol)
        result['aromatic_atoms'] = Descriptors.NumAromaticRings(mol)
        result['ring_count'] = Descriptors.RingCount(mol)

        # Save the image to a file
        img_path = "static/image.png"
        Draw.MolToFile(mol, img_path, size=(300, 300))
        result['img_path'] = img_path
    else:
        result['error'] = 'Invalid SMILES string'

    return result
