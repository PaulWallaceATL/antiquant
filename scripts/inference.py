import sys
import json
import os
import numpy as np
import xgboost as xgb
import pennylane as qml
from sklearn.decomposition import PCA
import lizard
from openai import OpenAI

# Load configuration
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
MODELS_DIR = os.path.join(os.path.dirname(__file__), "../models")
XGB_MODEL_PATH = os.path.join(MODELS_DIR, "xgb_model.json")
VQC_WEIGHTS_PATH = os.path.join(MODELS_DIR, "vqc_weights.npy")

# Initialize OpenAI
client = OpenAI(api_key=OPENAI_API_KEY)

def extract_features(code):
    # Use lizard to extract metrics
    # Create a temp file because lizard works best with files
    # Or use lizard.analyze_file.analyze_source_code
    
    analysis = lizard.analyze_file.analyze_source_code("snippet.py", code)
    
    loc = analysis.nloc
    num_functions = len(analysis.function_list)
    
    ccn = 0
    if num_functions > 0:
        ccn = max([f.cyclomatic_complexity for f in analysis.function_list])
    else:
        ccn = 1 # Default
        
    # Nesting depth approximation
    max_nesting = 0
    for line in code.split('\n'):
        stripped = line.lstrip()
        if stripped and not stripped.startswith('#'):
            indent = len(line) - len(stripped)
            depth = indent // 4
            if depth > max_nesting:
                max_nesting = depth
                
    return {
        "lines_of_code": loc,
        "num_functions": num_functions,
        "cyclomatic_complexity": ccn,
        "max_nesting_depth": max_nesting
    }

def get_embedding(code):
    if not OPENAI_API_KEY:
        # Return random embedding for testing if no key
        return np.random.rand(1536).tolist()
        
    response = client.embeddings.create(
        model="text-embedding-3-large",
        input=code[:8000]
    )
    return response.data[0].embedding

def predict_classical(features, embedding):
    # Load model
    model = xgb.XGBClassifier()
    model.load_model(XGB_MODEL_PATH)
    
    # Prepare input: [lines_of_code, num_functions, max_nesting_depth] + embedding
    # Note: We excluded cyclomatic_complexity from training in train_xgb.py
    
    feats = [
        features["lines_of_code"],
        features["num_functions"],
        features["max_nesting_depth"]
    ]
    
    X = np.hstack([feats, embedding]).reshape(1, -1)
    
    # Predict
    pred = model.predict(X)[0]
    probs = model.predict_proba(X)[0]
    
    labels = ["Low", "Medium", "High"]
    return {
        "prediction": labels[pred],
        "confidence": float(max(probs)),
        "model_used": "XGBoost"
    }

def predict_quantum(features, embedding):
    # Load weights
    if not os.path.exists(VQC_WEIGHTS_PATH):
        return {"error": "VQC weights not found"}
        
    data = np.load(VQC_WEIGHTS_PATH, allow_pickle=True).item()
    params = data['params']
    bias = data['bias']
    
    # PCA reduction (mocking the PCA transform from training)
    # In production, we should save the PCA model too.
    # For now, we'll just slice the first 4 dimensions or project randomly
    # to match the N_QUBITS=4.
    # Ideally: pca = load(PCA_PATH); x_pca = pca.transform(embedding)
    
    # Simple slice for prototype
    x_pca = np.array(embedding[:4]) 
    
    # Normalize
    x_norm = (x_pca - x_pca.min()) / (x_pca.max() - x_pca.min() + 1e-9)
    x_norm = x_norm * 2 * np.pi - np.pi
    
    # Quantum Circuit
    dev = qml.device("default.qubit", wires=4)
    
    @qml.qnode(dev)
    def circuit(params, x):
        for i in range(4):
            qml.RX(x[i], wires=i)
        for l in range(len(params)):
            for i in range(4):
                qml.Rot(params[l, i, 0], params[l, i, 1], params[l, i, 2], wires=i)
            for i in range(3):
                qml.CNOT(wires=[i, i + 1])
            if 4 > 1:
                qml.CNOT(wires=[3, 0])
        return qml.expval(qml.PauliZ(0))
        
    def variational_classifier(params, bias, x):
        return circuit(params, x) + bias
        
    # Predict
    raw_pred = variational_classifier(params, bias, x_norm)
    label = "High" if raw_pred > 0 else "Low" # Binary mapping
    
    return {
        "prediction": label,
        "confidence": float(abs(raw_pred)), # Pseudo-confidence
        "model_used": "PennyLane VQC"
    }

def main():
    try:
        input_data = json.loads(sys.stdin.read())
        code = input_data.get("code", "")
        mode = input_data.get("mode", "classical")
        
        features = extract_features(code)
        embedding = get_embedding(code)
        
        if mode == "quantum":
            result = predict_quantum(features, embedding)
        else:
            result = predict_classical(features, embedding)
            
        result["features"] = features
        print(json.dumps(result))
        
    except Exception as e:
        print(json.dumps({"error": str(e)}))

if __name__ == "__main__":
    main()
