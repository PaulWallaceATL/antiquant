import pandas as pd
import numpy as np
import json
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import os

# Paths
DATA_DIR = "data"
FEATURES_PATH = os.path.join(DATA_DIR, "features.csv")
EMBEDDINGS_PATH = os.path.join(DATA_DIR, "embeddings.json")
EMBEDDINGS_INDEX_PATH = os.path.join(DATA_DIR, "embeddings_index.csv")
MODEL_PATH = os.path.join("models", "xgb_model.json")

def load_data():
    print("Loading data...")
    features_df = pd.read_csv(FEATURES_PATH)
    
    # Load embeddings
    with open(EMBEDDINGS_PATH, 'r') as f:
        embeddings_list = json.load(f)
    embeddings = np.array(embeddings_list)
    
    # Load index to map embeddings to features
    index_df = pd.read_csv(EMBEDDINGS_INDEX_PATH)
    
    # Merge
    # features_df has 'id', index_df has 'id', 'index'
    merged = pd.merge(features_df, index_df, on='id')
    
    # Create X (features + embeddings) and y (complexity)
    # For complexity, we'll define classes based on Cyclomatic Complexity
    # Low (1-5), Medium (6-10), High (>10)
    
    X_features = merged[['lines_of_code', 'num_functions', 'cyclomatic_complexity', 'max_nesting_depth']].values
    
    # Get corresponding embeddings
    indices = merged['index'].values
    X_embeddings = embeddings[indices]
    
    # Concatenate
    X = np.hstack([X_features, X_embeddings])
    
    # Define target
    # We want to predict complexity class.
    # Let's use cyclomatic_complexity as the ground truth source, but bin it.
    # Note: This is a bit circular if we include cyclomatic_complexity in features, 
    # so usually we'd predict it from *code embeddings* alone or other features.
    # But the prompt says "Load features.csv + embeddings ... Train ... predicting complexity".
    # If we include 'cyclomatic_complexity' in input, accuracy will be 100%.
    # So we should EXCLUDE 'cyclomatic_complexity' from X if that's the target.
    # Let's exclude it from X.
    
    X_features_clean = merged[['lines_of_code', 'num_functions', 'max_nesting_depth']].values
    X = np.hstack([X_features_clean, X_embeddings])
    
    y_raw = merged['cyclomatic_complexity'].values
    y = []
    for val in y_raw:
        if val <= 5:
            y.append(0) # Low
        elif val <= 10:
            y.append(1) # Medium
        else:
            y.append(2) # High
            
    return X, np.array(y)

def train():
    X, y = load_data()
    
    # Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    print(f"Training on {len(X_train)} samples, testing on {len(X_test)}...")
    
    # Train XGBoost
    model = xgb.XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
    model.fit(X_train, y_train)
    
    # Evaluate
    preds = model.predict(X_test)
    acc = accuracy_score(y_test, preds)
    print(f"Accuracy: {acc:.4f}")
    print("Classification Report:")
    print(classification_report(y_test, preds))
    
    # Save model
    os.makedirs("models", exist_ok=True)
    model.save_model(MODEL_PATH)
    print(f"Model saved to {MODEL_PATH}")

if __name__ == "__main__":
    train()
