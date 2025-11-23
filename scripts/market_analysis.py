import sys
import json
import os
import random
from rdkit import Chem
from rdkit.Chem import Descriptors, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

# Known drug database for competitor matching (mock data)
KNOWN_DRUGS = [
    {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "indication": "Pain relief, Anti-inflammatory", "status": "Approved", "company": "Bayer"},
    {"name": "Ibuprofen", "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "indication": "Pain relief, Anti-inflammatory", "status": "Approved", "company": "Various"},
    {"name": "Acetaminophen", "smiles": "CC(=O)Nc1ccc(O)cc1", "indication": "Pain relief, Fever reduction", "status": "Approved", "company": "Various"},
    {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "indication": "Stimulant", "status": "Approved", "company": "Various"},
    {"name": "Metformin", "smiles": "CN(C)C(=N)N=C(N)N", "indication": "Type 2 Diabetes", "status": "Approved", "company": "Various"},
    {"name": "Atorvastatin", "smiles": "CC(C)C(=O)Nc1c(C(C)C(=O)O)cccc1c1ccc(O)cc1", "indication": "Cholesterol reduction", "status": "Approved", "company": "Pfizer"},
    {"name": "Lisinopril", "smiles": "CCCCN1CCCC1C(=O)N2CCCC2C(=O)N3CCCC3C(=O)O", "indication": "Hypertension", "status": "Approved", "company": "Various"},
    {"name": "Amlodipine", "smiles": "CCOC(=O)C1=C(C)NC(=C(C1C2=CC=CC=C2Cl)C(=O)OCC)C", "indication": "Hypertension, Angina", "status": "Approved", "company": "Pfizer"},
]

def calculate_similarity(smiles1, smiles2):
    """Calculate Tanimoto similarity between two SMILES strings"""
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 is None or mol2 is None:
            return 0.0
        
        fp1 = FingerprintMols.FingerprintMol(mol1)
        fp2 = FingerprintMols.FingerprintMol(mol2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except:
        return 0.0

def get_market_analysis(smiles, molecule_name=None):
    """Generate market analysis data based on molecular properties"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Calculate molecular properties for heuristics
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        num_rings = Descriptors.RingCount(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        
        # Find similar competitors
        competitors = []
        for drug in KNOWN_DRUGS:
            similarity = calculate_similarity(smiles, drug["smiles"])
            if similarity > 0.3:  # Threshold for relevance
                competitors.append({
                    "name": drug["name"],
                    "smiles": drug["smiles"],
                    "similarity": round(similarity, 3),
                    "indication": drug["indication"],
                    "status": drug["status"],
                    "company": drug["company"]
                })
        
        # Sort by similarity
        competitors.sort(key=lambda x: x["similarity"], reverse=True)
        competitors = competitors[:5]  # Top 5
        
        # Generate mock patent data based on molecular properties
        patents = []
        if mw < 300:
            patents.append({
                "patentNumber": f"US{random.randint(9000000, 9999999)}B2",
                "title": "Small molecule therapeutic compounds and methods of use",
                "filingDate": "2020-03-15",
                "status": "Granted",
                "assignee": "PharmaCorp Inc.",
                "relevance": "High"
            })
        if num_aromatic_rings > 0:
            patents.append({
                "patentNumber": f"US{random.randint(9000000, 9999999)}B2",
                "title": "Aromatic compound derivatives for pharmaceutical applications",
                "filingDate": "2019-07-22",
                "status": "Granted",
                "assignee": "BioTech Solutions LLC",
                "relevance": "Medium"
            })
        if logp > 3:
            patents.append({
                "patentNumber": f"US{random.randint(9000000, 9999999)}B2",
                "title": "Lipophilic compounds and pharmaceutical compositions",
                "filingDate": "2021-11-08",
                "status": "Pending",
                "assignee": "MedPharm Industries",
                "relevance": "Medium"
            })
        
        # Determine therapeutic area based on properties (heuristic)
        therapeutic_areas = []
        if mw < 400 and logp < 2:
            therapeutic_areas.append("Cardiovascular")
        if num_aromatic_rings > 1:
            therapeutic_areas.append("Oncology")
        if logp > 3:
            therapeutic_areas.append("CNS Disorders")
        if mw < 300:
            therapeutic_areas.append("Metabolic Disorders")
        
        therapeutic_area = therapeutic_areas[0] if therapeutic_areas else "General Therapeutics"
        
        # Estimate market size based on therapeutic area (mock data)
        market_sizes = {
            "Cardiovascular": {"size": "45.2", "growth": "5.2%"},
            "Oncology": {"size": "185.3", "growth": "8.7%"},
            "CNS Disorders": {"size": "98.5", "growth": "6.1%"},
            "Metabolic Disorders": {"size": "67.8", "growth": "4.9%"},
            "General Therapeutics": {"size": "32.1", "growth": "3.5%"}
        }
        
        market_info = market_sizes.get(therapeutic_area, market_sizes["General Therapeutics"])
        
        # Regulatory status based on drug-likeness
        if mw <= 500 and logp <= 5:
            regulatory_status = "Favorable - Meets key drug-likeness criteria"
        elif mw <= 600 and logp <= 6:
            regulatory_status = "Moderate - Some optimization may be needed"
        else:
            regulatory_status = "Challenging - Significant optimization required"
        
        return {
            "competitors": competitors,
            "patents": patents[:3],  # Top 3 most relevant
            "marketSize": {
                "therapeuticArea": therapeutic_area,
                "estimatedMarketSize": market_info["size"],
                "currency": "USD",
                "year": 2024,
                "growthRate": market_info["growth"]
            },
            "regulatoryStatus": regulatory_status
        }
        
    except Exception as e:
        print(f"Market analysis error: {e}", file=sys.stderr)
        return None

def main():
    """CLI interface for testing"""
    if len(sys.argv) < 2:
        print("Usage: python market_analysis.py <SMILES>")
        sys.exit(1)
    
    smiles = sys.argv[1]
    result = get_market_analysis(smiles)
    
    if result:
        print(json.dumps(result, indent=2))
    else:
        print(json.dumps({"error": "Failed to generate market analysis"}))

if __name__ == "__main__":
    main()

