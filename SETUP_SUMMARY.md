# 🎉 Cloud Run Setup Complete!

## ✅ What Was Created

### 1. **FastAPI Backend** (`api/main.py`)
- REST API wrapper around your Python inference scripts
- Endpoints: `/analyze`, `/analyze/classical`, `/analyze/quantum`
- CORS enabled for Vercel frontend
- Health check endpoint: `/health`

### 2. **Docker Configuration**
- `api/Dockerfile` - Optimized container for Cloud Run
- `api/requirements.txt` - All Python dependencies
- `api/.dockerignore` - Excludes unnecessary files
- `api/.gcloudignore` - Excludes files from Cloud Run deployment

### 3. **Frontend Integration**
- Updated `web/lib/molecular_inference.ts` to call Cloud Run API
- Automatically uses Cloud Run if `MOLECULEAI_API_URL` is set
- Falls back to local Python for development

### 4. **Documentation**
- `api/cloud-run-deploy.md` - Detailed deployment guide
- `CLOUD_RUN_INTEGRATION.md` - Frontend integration guide  
- `README_CLOUD_RUN.md` - Quick start guide

---

## 🚀 Next Steps

### Step 1: Deploy to Cloud Run

```bash
# Install Google Cloud CLI (if not already installed)
brew install google-cloud-sdk

# Authenticate
gcloud auth login

# Create project (or use existing)
gcloud projects create YOUR_PROJECT_ID
gcloud config set project YOUR_PROJECT_ID

# Enable APIs
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com

# Set up billing (required even for free tier)
# Go to: https://console.cloud.google.com/billing

# Deploy!
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

**You'll get a URL like:** `https://moleculeai-api-xxxxx-uc.a.run.app`

### Step 2: Update Vercel

1. Go to **Vercel Dashboard** → Your Project → **Settings** → **Environment Variables**
2. Add:
   ```
   MOLECULEAI_API_URL=https://your-api-url-uc.a.run.app
   ```
3. **Redeploy** your Next.js app (or push to trigger auto-deploy)

### Step 3: Test

1. Test Cloud Run API:
   ```bash
   curl https://your-api-url-uc.a.run.app/health
   ```

2. Test your Vercel app - it should now call Cloud Run!

---

## 💰 Cost Estimate

**Free Tier Coverage:**
- ✅ 2M requests/month (you'll use ~15K)
- ✅ 180K vCPU-seconds/month (you'll use ~125K)
- ✅ 360K GB-seconds/month (you'll use ~250K)

**Estimated monthly cost: $0** 🎉

---

## 📚 Documentation

- **Full deployment guide:** `api/cloud-run-deploy.md`
- **Frontend integration:** `CLOUD_RUN_INTEGRATION.md`
- **Quick reference:** `README_CLOUD_RUN.md`

---

## 🐛 Troubleshooting

### Build fails
- Check `api/Dockerfile` syntax
- Verify all dependencies in `api/requirements.txt`

### Deployment fails
- Ensure billing is enabled (required even for free tier)
- Check APIs are enabled: `run.googleapis.com`, `cloudbuild.googleapis.com`

### API doesn't work
- Check logs: `gcloud run logs read moleculeai-api`
- Verify `OPENAI_API_KEY` is set
- Test health endpoint: `curl https://your-url/health`

### Frontend still uses local Python
- Verify `MOLECULEAI_API_URL` is set in Vercel
- Check environment variable name matches exactly
- Redeploy after setting env var

---

## ✨ Benefits

1. ✅ **Separate scaling** - Backend scales independently
2. ✅ **Free tier** - Likely $0/month for MVP
3. ✅ **Fast deployment** - Deploy in minutes
4. ✅ **Auto-scaling** - Scales to zero when idle
5. ✅ **HTTPS** - Automatic SSL certificate
6. ✅ **Monitoring** - Built-in logs and metrics

---

## 🎯 You're Ready!

Follow the steps above and you'll have your Python backend running on Cloud Run in minutes!

**Questions?** Check the detailed docs or Google Cloud Run documentation.

