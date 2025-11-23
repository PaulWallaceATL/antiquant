import sys
import json
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

DATA_DIR = os.path.join(os.path.dirname(__file__), "../data")
MOLECULES_PATH = os.path.join(DATA_DIR, "molecules_sample.csv")

def validate_smiles(smiles):
    """Validate SMILES string using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def find_similar(smiles, database_path=None, top_k=5):
    """
    Find similar molecules in the database using Tanimoto similarity
    
    Args:
        smiles: Input SMILES string
        database_path: Path to CSV file with molecules (default: data/molecules_sample.csv)
        top_k: Number of top similar molecules to return
    
    Returns:
        List of dictionaries with similar molecule info
    """
    try:
        # Load database
        if database_path is None:
            database_path = MOLECULES_PATH
        
        if not os.path.exists(database_path):
            return []
        
        df = pd.read_csv(database_path)
        
        # Validate input SMILES
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            return []
        
        input_fp = FingerprintMols.FingerprintMol(input_mol)
        
        # Calculate similarities
        similarities = []
        for idx, row in df.iterrows():
            try:
                db_smiles = row.get('smiles', '')
                if not db_smiles:
                    continue
                
                db_mol = Chem.MolFromSmiles(db_smiles)
                if db_mol is None:
                    continue
                
                db_fp = FingerprintMols.FingerprintMol(db_mol)
                similarity = DataStructs.TanimotoSimilarity(input_fp, db_fp)
                
                similarities.append({
                    'id': row.get('id', f'mol_{idx}'),
                    'smiles': db_smiles,
                    'similarity': float(similarity),
                    'target_property': row.get('target_property', None)
                })
            except Exception as e:
                continue
        
        # Sort by similarity and return top K
        similarities.sort(key=lambda x: x['similarity'], reverse=True)
        return similarities[:top_k]
        
    except Exception as e:
        print(f"Error in find_similar: {e}", file=sys.stderr)
        return []

def main():
    """CLI interface for testing"""
    if len(sys.argv) < 3:
        print("Usage: python chat_helpers.py <command> <smiles> [top_k]")
        print("Commands: validate, find_similar")
        sys.exit(1)
    
    command = sys.argv[1]
    smiles = sys.argv[2]
    
    if command == "validate":
        result = validate_smiles(smiles)
        print(json.dumps({"valid": result}))
    elif command == "find_similar":
        top_k = int(sys.argv[3]) if len(sys.argv) > 3 else 5
        result = find_similar(smiles, top_k=top_k)
        print(json.dumps(result))
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

if __name__ == "__main__":
    main()

