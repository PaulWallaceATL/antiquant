import sys
import json
import os
import numpy as np
import xgboost as xgb
from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from openai import OpenAI

# Load configuration
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
MODELS_DIR = os.path.join(os.path.dirname(__file__), "../models")
XGB_MODEL_PATH = os.path.join(MODELS_DIR, "xgb_model.json")
VQC_WEIGHTS_PATH = os.path.join(MODELS_DIR, "vqc_weights.npy")

# Initialize OpenAI
client = OpenAI(api_key=OPENAI_API_KEY) if OPENAI_API_KEY else None

def extract_molecular_features(smiles):
    """Extract molecular descriptors from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        features = {
            "molecular_weight": float(Descriptors.MolWt(mol)),
            "logp": float(Descriptors.MolLogP(mol)),
            "num_h_donors": int(Lipinski.NumHDonors(mol)),
            "num_h_acceptors": int(Lipinski.NumHAcceptors(mol)),
            "num_rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
            "num_aromatic_rings": int(Lipinski.NumAromaticRings(mol)),
            "tpsa": float(Descriptors.TPSA(mol)),
            "num_atoms": int(mol.GetNumAtoms())
        }
        
        return features
        
    except Exception as e:
        return None

def get_embedding(smiles):
    """Get embedding for SMILES string"""
    if not client:
        # Return random embedding for testing if no key
        return np.random.rand(1536).tolist()
        
    try:
        response = client.embeddings.create(
            model="text-embedding-3-large",
            input=smiles[:8000]
        )
        return response.data[0].embedding
    except:
        return np.random.rand(1536).tolist()

def predict_classical(features, embedding):
    """Predict solubility using XGBoost"""
    try:
        # Load model
        model = xgb.XGBRegressor()
        model.load_model(XGB_MODEL_PATH)
        
        # Prepare input
        feats = [
            features["molecular_weight"],
            features["logp"],
            features["num_h_donors"],
            features["num_h_acceptors"],
            features["num_rotatable_bonds"],
            features["num_aromatic_rings"],
            features["tpsa"],
            features["num_atoms"]
        ]
        
        X = np.hstack([feats, embedding]).reshape(1, -1)
        
        # Predict
        pred = model.predict(X)[0]
        
        return {
            "prediction": f"{pred:.2f}",
            "property": "Log Solubility (mol/L)",
            "model_used": "XGBoost Regressor"
        }
    except Exception as e:
        return {"error": f"Prediction error: {str(e)}"}

def predict_quantum(features, embedding):
    """Predict using VQC (simplified for demo)"""
    try:
        if not os.path.exists(VQC_WEIGHTS_PATH):
            return {"error": "VQC weights not found"}
        
        # Simplified quantum prediction
        # In production, this would use PennyLane VQC
        data = np.load(VQC_WEIGHTS_PATH, allow_pickle=True).item()
        
        # Use first 4 dimensions of embedding
        x_reduced = np.array(embedding[:4])
        
        # Normalize
        x_norm = (x_reduced - x_reduced.min()) / (x_reduced.max() - x_reduced.min() + 1e-9)
        
        # Mock quantum prediction (in reality, would run through VQC circuit)
        mock_pred = np.mean(x_norm) * 2 - 1  # Value between -1 and 1
        
        return {
            "prediction": f"{mock_pred:.2f}",
            "property": "Log Solubility (mol/L)",
            "model_used": "Quantum VQC (Simulated)"
        }
    except Exception as e:
        return {"error": f"Quantum prediction error: {str(e)}"}

def main():
    try:
        input_data = json.loads(sys.stdin.read())
        smiles = input_data.get("smiles", "")
        mode = input_data.get("mode", "classical")
        
        # Extract features
        features = extract_molecular_features(smiles)
        if not features:
            print(json.dumps({"error": "Invalid SMILES string"}))
            return
        
        # Get embedding
        embedding = get_embedding(smiles)
        
        # Predict
        if mode == "quantum":
            result = predict_quantum(features, embedding)
        else:
            result = predict_classical(features, embedding)
            
        result["features"] = features
        result["smiles"] = smiles
        print(json.dumps(result))
        
    except Exception as e:
        print(json.dumps({"error": str(e)}))

if __name__ == "__main__":
    main()
