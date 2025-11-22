import os
import pandas as pd
import requests
from io import StringIO

def download_esol():
    """Download ESOL (solubility) dataset"""
    print("Downloading ESOL dataset...")
    
    # ESOL dataset from DeepChem
    url = "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        df = pd.read_csv(StringIO(response.text))
        
        # Rename columns to match our format
        # ESOL has: SMILES, measured log solubility in mols per litre
        df = df.rename(columns={
            'smiles': 'smiles',
            'measured log solubility in mols per litre': 'target_property'
        })
        
        # Add ID column
        df['id'] = ['mol_' + str(i) for i in range(len(df))]
        
        # Select relevant columns
        df = df[['id', 'smiles', 'target_property']]
        
        # Take subset (first 1000 molecules)
        df = df.head(1000)
        
        # Save to data directory
        os.makedirs('data', exist_ok=True)
        output_path = 'data/molecules_sample.csv'
        df.to_csv(output_path, index=False)
        
        print(f"✓ Downloaded {len(df)} molecules")
        print(f"✓ Saved to {output_path}")
        print(f"Property: Log Solubility (mol/L)")
        
        return df
        
    except Exception as e:
        print(f"Error downloading ESOL: {e}")
        print("Creating sample dataset instead...")
        
        # Fallback: Create a small sample dataset
        sample_data = {
            'id': ['mol_0', 'mol_1', 'mol_2'],
            'smiles': [
                'CCO',  # Ethanol
                'CC(=O)O',  # Acetic acid
                'c1ccccc1'  # Benzene
            ],
            'target_property': [-0.77, -0.17, -2.13]  # Log solubility values
        }
        df = pd.DataFrame(sample_data)
        
        os.makedirs('data', exist_ok=True)
        output_path = 'data/molecules_sample.csv'
        df.to_csv(output_path, index=False)
        
        print(f"✓ Created {len(df)} sample molecules")
        print(f"✓ Saved to {output_path}")
        
        return df

if __name__ == "__main__":
    download_esol()
