"""
FastAPI wrapper for molecular inference API
Deploy to Google Cloud Run
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, Literal
import sys
import os

# Add parent directory to path to import scripts
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from scripts.molecular_inference import (
    extract_molecular_features,
    calculate_drug_likeness,
    get_3d_structure_enhanced,
    get_embedding,
    predict_classical,
    predict_quantum
)
from scripts.market_analysis import get_market_analysis
from rdkit import Chem

app = FastAPI(title="MoleculeAI API", version="1.0.0")

# CORS middleware - allow Vercel frontend to call this
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your Vercel domain
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class AnalyzeRequest(BaseModel):
    smiles: str
    mode: Optional[Literal["classical", "quantum", "both"]] = "both"

@app.get("/")
def root():
    return {"status": "ok", "service": "MoleculeAI API"}

@app.get("/health")
def health():
    return {"status": "healthy"}

@app.post("/analyze")
async def analyze(request: AnalyzeRequest):
    """Main analysis endpoint - supports classical, quantum, or both"""
    try:
        smiles = request.smiles.strip()
        mode = request.mode or "both"
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Extract features
        features = extract_molecular_features(smiles)
        if not features:
            raise HTTPException(status_code=500, detail="Feature extraction failed")
        
        # Calculate drug-likeness
        drug_likeness = calculate_drug_likeness(mol)
        if not drug_likeness:
            raise HTTPException(status_code=500, detail="Drug-likeness calculation failed")
        
        # Get market analysis
        market_analysis = get_market_analysis(smiles)
        
        # Get enhanced 3D structure
        structure_3d = get_3d_structure_enhanced(smiles)
        
        # Get embedding
        embedding = get_embedding(smiles)
        
        # Predict based on mode
        if mode == "both":
            classical_result = predict_classical(features, embedding)
            quantum_result = predict_quantum(features, embedding)
            
            result = {
                "features": features,
                "smiles": smiles,
                "structure_3d": structure_3d,
                "drug_likeness": drug_likeness,
                "marketAnalysis": market_analysis,
                "classical": classical_result,
                "quantum": quantum_result
            }
        elif mode == "quantum":
            quantum_result = predict_quantum(features, embedding)
            result = {
                "features": features,
                "smiles": smiles,
                "structure_3d": structure_3d,
                "drug_likeness": drug_likeness,
                "marketAnalysis": market_analysis,
                **quantum_result
            }
        else:  # classical
            classical_result = predict_classical(features, embedding)
            result = {
                "features": features,
                "smiles": smiles,
                "structure_3d": structure_3d,
                "drug_likeness": drug_likeness,
                "marketAnalysis": market_analysis,
                **classical_result
            }
        
        return result
        
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        print(f"Error in analyze: {e}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/analyze/classical")
async def analyze_classical(request: AnalyzeRequest):
    """Classical XGBoost prediction only"""
    request.mode = "classical"
    return await analyze(request)

@app.post("/analyze/quantum")
async def analyze_quantum(request: AnalyzeRequest):
    """Quantum VQC prediction only"""
    request.mode = "quantum"
    return await analyze(request)

if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", 8080))
    uvicorn.run(app, host="0.0.0.0", port=port)

