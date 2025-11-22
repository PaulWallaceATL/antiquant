import numpy as np
import os

MODELS_DIR = "models"
WEIGHTS_PATH = os.path.join(MODELS_DIR, "vqc_weights.npy")

def generate_weights():
    # Matches the shape in the notebook/inference script
    # N_LAYERS = 2, N_QUBITS = 4, 3 params per rot
    params = np.random.randn(2, 4, 3)
    bias = np.array(0.0)
    
    os.makedirs(MODELS_DIR, exist_ok=True)
    np.save(WEIGHTS_PATH, {'params': params, 'bias': bias})
    print(f"Saved mock VQC weights to {WEIGHTS_PATH}")

if __name__ == "__main__":
    generate_weights()
