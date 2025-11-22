# MoleculeAI - Hybrid Classical + Quantum Molecular Property Predictor

A production-ready web application that predicts molecular properties (solubility) from SMILES strings using a hybrid classical ML + quantum-inspired approach.

## 🧬 Features

- **SMILES Input**: Paste molecular SMILES notation for instant property prediction
- **Classical ML**: XGBoost regressor trained on molecular descriptors + embeddings
- **Quantum VQC**: Variational Quantum Classifier using PennyLane's quantum simulator
- **Dual Predictions**: Compare classical vs quantum model results
- **Molecular Descriptors**: RDKit-powered feature extraction (MW, LogP, H-bonds, TPSA, etc.)
- **Beautiful UI**: Modern Next.js interface with real-time visualization

## 🏗️ Architecture

```
├── web/                    # Next.js App Router frontend
│   ├── app/
│   │   ├── api/           # API routes (/analyze, /analyze-quantum)
│   │   └── page.tsx       # Main UI
│   └── lib/               # Inference wrappers
├── scripts/               # Data & ML pipelines
│   ├── ingest_molecules.py          # Download ESOL dataset
│   ├── extract_molecular_features.py # RDKit descriptors
│   ├── mock_embed.ts                # Mock embeddings (or use OpenAI)
│   └── molecular_inference.py       # Python inference engine
├── models/                # Trained models
│   ├── train_molecular_xgb.py      # XGBoost training
│   └── xgb_model.json              # Saved model
├── notebooks/             # Jupyter experiments
│   └── 02_vqc_pennylane.ipynb      # Quantum VQC training
└── data/                  # Datasets
    ├── molecules_sample.csv         # ESOL solubility data
    ├── molecular_features.csv       # RDKit descriptors
    └── embeddings.json             # SMILES embeddings
```

## 🚀 Getting Started

### Prerequisites

- **Node.js** 18+
- **Python** 3.9+
- **OpenAI API Key** (optional, for production embeddings)

### Installation

1. **Install JavaScript dependencies:**
   ```bash
   npm install
   cd web && npm install && cd ..
   ```

2. **Install Python dependencies:**
   ```bash
   pip install rdkit-pypi xgboost pennylane scikit-learn pandas numpy openai requests
   ```

3. **Set up environment variables:**
   ```bash
   cp .env.example .env
   # Add your OPENAI_API_KEY (optional for testing)
   ```

### Running Locally

1. **Download and prepare dataset** (automatic):
   ```bash
   # Download ESOL dataset (~1000 molecules)
   python3 scripts/ingest_molecules.py
   
   # Extract RDKit molecular descriptors
   python3 scripts/extract_molecular_features.py
   
   # Generate embeddings (mock or OpenAI)
   npx ts-node scripts/mock_embed.ts
   ```

2. **Train the XGBoost model:**
   ```bash
   python3 models/train_molecular_xgb.py
   ```

3. **Run the web app:**
   ```bash
   cd web
   npm run dev
   ```
   
   Visit **http://localhost:3000**

## 🎯 Usage

1. **Enter a SMILES string** (e.g., `CCO` for ethanol, `c1ccccc1` for benzene)
2. **Select model**: Classical (XGBoost) or Quantum (VQC)
3. **Click "Predict Solubility"**
4. **View results**: Log solubility prediction + molecular descriptors

### Example SMILES

- `CCO` - Ethanol
- `CC(=O)Oc1ccccc1C(=O)O` - Aspirin
- `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` - Caffeine
- `c1ccccc1` - Benzene

## ☁️ Deployment on Vercel

1. **Push to GitHub**
   
2. **Import to Vercel**:
   - Set **Root Directory** to `web`
   - Framework Preset: **Next.js**

3. **Add Environment Variables** in Vercel dashboard:
   ```
   OPENAI_API_KEY=sk-...
   APP_SECRET=your-secret
   JWT_SECRET=your-jwt
   ```
   
   (Supabase vars are optional for production data storage)

4. **Deploy!**

### ⚠️ Production Notes

For serverless deployment (Vercel), the Python inference script must be compatible with the runtime. Heavy dependencies like `rdkit` and `pennylane` may exceed bundle limits. Consider:

- **Option A**: Deploy Python API separately (Google Cloud Run, AWS Lambda with containers)
- **Option B**: Use Vercel's Edge Runtime with lighter alternatives
- **Option C**: Pre-compute predictions and serve from Supabase

## 📊 Dataset

**ESOL (Delaney Solubility Dataset)**
- ~1,128 molecules with measured aqueous solubility
- Property: Log Solubility in mol/L
- Source: [DeepChem/ESOL](https://github.com/deepchem/deepchem/tree/master/datasets)

### Citation

```
Delaney, J. S. "ESOL: Estimating Aqueous Solubility Directly from Molecular Structure." 
Journal of Chemical Information and Computer Sciences, 2004.
```

## 🧪 Tech Stack

- **Frontend**: Next.js 14, TailwindCSS, Recharts, Lucide Icons
- **Backend**: Node.js API Routes (Vercel Functions)
- **ML Classical**: XGBoost, RDKit, scikit-learn
- **ML Quantum**: PennyLane (default.qubit simulator)
- **Embeddings**: OpenAI text-embedding-3-large
- **Storage**: Supabase (optional)

## 📝 License

MIT

## 🙏 Acknowledgments

- ESOL dataset from DeepChem
- RDKit for molecular descriptor calculation
- PennyLane for quantum computing framework
