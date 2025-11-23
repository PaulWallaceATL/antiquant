# 🚀 Google Cloud Run Deployment

Deploy your Python backend separately on Google Cloud Run while keeping your Next.js frontend on Vercel.

## ✅ Free Tier (You'll Likely Pay $0!)

Google Cloud Run has a generous free tier:
- **2 million requests/month** - FREE
- **180,000 vCPU-seconds/month** (~50 hours CPU time) - FREE  
- **360,000 GB-seconds/month** memory - FREE
- **1 GB egress/month** (North America) - FREE

**Your estimated usage:** ~15,000 requests/month = **$0/month** 🎉

---

## 📁 Project Structure

```
quant-code-complexity/
├── api/                          # Cloud Run API
│   ├── main.py                   # FastAPI app
│   ├── Dockerfile                # Container config
│   ├── requirements.txt          # Python dependencies
│   └── cloud-run-deploy.md      # Detailed deployment guide
├── scripts/                      # Python inference scripts
├── models/                       # ML model files
├── web/                          # Next.js frontend (Vercel)
└── CLOUD_RUN_INTEGRATION.md     # Frontend integration guide
```

---

## 🚀 Quick Start

### 1. Install Google Cloud CLI

```bash
# macOS
brew install google-cloud-sdk

# Or download from: https://cloud.google.com/sdk/docs/install
```

### 2. Authenticate & Create Project

```bash
gcloud auth login
gcloud projects create your-project-id
gcloud config set project your-project-id
```

### 3. Enable APIs

```bash
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com
```

### 4. Deploy to Cloud Run

From the project root:

```bash
gcloud run deploy moleculeai-api \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10 \
  --set-env-vars OPENAI_API_KEY=sk-your-key-here
```

**You'll get a URL like:**
```
https://moleculeai-api-xxxxx-uc.a.run.app
```

### 5. Update Frontend

In **Vercel Dashboard** → Settings → Environment Variables:
```
MOLECULEAI_API_URL=https://your-api-url-uc.a.run.app
```

The frontend code is already updated to use Cloud Run when this env var is set!

---

## 📚 Detailed Guides

- **`api/cloud-run-deploy.md`** - Full deployment instructions
- **`CLOUD_RUN_INTEGRATION.md`** - Frontend integration details

---

## 🧪 Test Locally

### Option 1: Test FastAPI Directly

```bash
cd api
pip install -r requirements.txt
python main.py

# Test
curl http://localhost:8080/health
curl -X POST http://localhost:8080/analyze \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "mode": "both"}'
```

### Option 2: Test with Docker

```bash
docker build -t moleculeai-api -f api/Dockerfile .
docker run -p 8080:8080 -e OPENAI_API_KEY=sk-... moleculeai-api
```

---

## 💰 Cost Estimate

**Free tier usage:**
- Requests: 15,000/month → **$0** (under 2M limit)
- CPU: ~125,000 vCPU-seconds/month → **$0** (under 180K limit)
- Memory: ~250,000 GB-seconds/month → **$0** (under 360K limit)

**Total: $0/month** ✅

Even if you exceed free tier, typical costs are $0-5/month for small apps.

---

## 🔧 Troubleshooting

### Build fails
- Check `api/Dockerfile` and `api/requirements.txt`
- Verify RDKit dependencies are installed

### Deployment fails  
- Ensure billing is enabled (required even for free tier)
- Check APIs are enabled: `run.googleapis.com`, `cloudbuild.googleapis.com`

### Runtime errors
- Check logs: `gcloud run logs read moleculeai-api`
- Verify `OPENAI_API_KEY` is set
- Ensure model files exist in `/models` directory

### CORS errors
- Update CORS settings in `api/main.py`
- Add your Vercel domain to allowed origins

---

## 📊 Monitoring

View logs and metrics:
```bash
# View logs
gcloud run logs read moleculeai-api --limit 50

# Stream logs
gcloud run logs tail moleculeai-api

# View in console
https://console.cloud.google.com/run
```

---

## 🎯 Next Steps

1. ✅ Deploy to Cloud Run
2. ✅ Get API URL
3. ✅ Set `MOLECULEAI_API_URL` in Vercel
4. ✅ Test integration
5. ✅ Monitor usage

**You're all set!** 🎉

