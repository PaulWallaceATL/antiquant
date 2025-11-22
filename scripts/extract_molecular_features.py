import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import numpy as np

def calculate_molecular_descriptors(smiles):
    """Calculate molecular descriptors from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        descriptors = {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'num_h_donors': Lipinski.NumHDonors(mol),
            'num_h_acceptors': Lipinski.NumHAcceptors(mol),
            'num_rotatable_bonds': Lipinski.NumRotatableBonds(mol),
            'num_aromatic_rings': Lipinski.NumAromaticRings(mol),
            'tpsa': Descriptors.TPSA(mol),
            'num_atoms': mol.GetNumAtoms()
        }
        
        return descriptors
        
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return None

def extract_features_from_csv(input_csv='data/molecules_sample.csv', 
                               output_csv='data/molecular_features.csv'):
    """Extract molecular features from CSV of molecules"""
    
    print(f"Reading molecules from {input_csv}...")
    df = pd.read_csv(input_csv)
    
    features_list = []
    
    for idx, row in df.iterrows():
        mol_id = row['id']
        smiles = row['smiles']
        
        descriptors = calculate_molecular_descriptors(smiles)
        
        if descriptors:
            descriptors['id'] = mol_id
            features_list.append(descriptors)
        else:
            print(f"Skipping invalid molecule: {mol_id}")
            
        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1} molecules...")
    
    # Create DataFrame
    features_df = pd.DataFrame(features_list)
    
    # Reorder columns to have id first
    cols = ['id'] + [col for col in features_df.columns if col != 'id']
    features_df = features_df[cols]
    
    # Save
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    features_df.to_csv(output_csv, index=False)
    
    print(f"✓ Extracted features for {len(features_df)} molecules")
    print(f"✓ Saved to {output_csv}")
    print(f"Features: {list(features_df.columns)}")
    
    return features_df

if __name__ == "__main__":
    extract_features_from_csv()
