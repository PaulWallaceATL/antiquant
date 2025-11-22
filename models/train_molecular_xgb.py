import pandas as pd
import numpy as np
import json
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import os

# Paths
DATA_DIR = "data"
FEATURES_PATH = os.path.join(DATA_DIR, "molecular_features.csv")
MOLECULES_PATH = os.path.join(DATA_DIR, "molecules_sample.csv")
EMBEDDINGS_PATH = os.path.join(DATA_DIR, "embeddings.json")
EMBEDDINGS_INDEX_PATH = os.path.join(DATA_DIR, "embeddings_index.csv")
MODEL_PATH = os.path.join("models", "xgb_model.json")

def load_data():
    print("Loading data...")
    features_df = pd.read_csv(FEATURES_PATH)
    molecules_df = pd.read_csv(MOLECULES_PATH)
    
    # Load embeddings
    with open(EMBEDDINGS_PATH, 'r') as f:
        embeddings_list = json.load(f)
    embeddings = np.array(embeddings_list)
    
    # Load index
    index_df = pd.read_csv(EMBEDDINGS_INDEX_PATH)
    
    # Merge features with molecules to get target
    merged = pd.merge(features_df, molecules_df[['id', 'target_property']], on='id')
    merged = pd.merge(merged, index_df, on='id')
    
    # Features: molecular descriptors
    feature_cols = ['molecular_weight', 'logp', 'num_h_donors', 'num_h_acceptors', 
                    'num_rotatable_bonds', 'num_aromatic_rings', 'tpsa', 'num_atoms']
    X_features = merged[feature_cols].values
    
    # Get corresponding embeddings
    indices = merged['index'].values
    X_embeddings = embeddings[indices]
    
    # Concatenate features + embeddings
    X = np.hstack([X_features, X_embeddings])
    
    # Target: solubility (regression task now)
    y = merged['target_property'].values
    
    return X, y

def train():
    X, y = load_data()
    
    # Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    print(f"Training on {len(X_train)} molecules, testing on {len(X_test)}...")
    
    # Train XGBoost Regressor (changed from Classifier)
    model = xgb.XGBRegressor(n_estimators=100, max_depth=5, learning_rate=0.1)
    model.fit(X_train, y_train)
    
    # Evaluate
    preds = model.predict(X_test)
    mse = mean_squared_error(y_test, preds)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_test, preds)
    r2 = r2_score(y_test, preds)
    
    print(f"RMSE: {rmse:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"R²: {r2:.4f}")
    
    # Save model
    os.makedirs("models", exist_ok=True)
    model.save_model(MODEL_PATH)
    print(f"Model saved to {MODEL_PATH}")

if __name__ == "__main__":
    train()
