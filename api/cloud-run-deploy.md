# Deploy to Google Cloud Run

## ✅ Free Tier Limits (Plenty for MVP)

- **2 million requests/month** - FREE
- **180,000 vCPU-seconds/month** (~50 hours CPU time) - FREE
- **360,000 GB-seconds/month** memory - FREE
- **1 GB egress/month** (North America) - FREE

**Estimated usage for your app:**
- ~100-500 requests/day = ~15,000/month (well under 2M limit)
- ~5-10 seconds per request = ~125,000 vCPU-seconds/month (under 180K limit)

**You'll likely pay $0/month!** 🎉

---

## 📋 Prerequisites

1. **Google Cloud Account** (free $300 credit for new accounts)
2. **Google Cloud CLI** installed
3. **Docker** (for local testing)

---

## 🚀 Step-by-Step Deployment

### 1. Install Google Cloud CLI

**macOS:**
```bash
brew install google-cloud-sdk
```

**Or download from:** https://cloud.google.com/sdk/docs/install

### 2. Authenticate & Set Up Project

**First, login:**
```bash
gcloud auth login
```

This will open your browser to authenticate with your Google Cloud account.

**Then, check your existing projects:**
```bash
gcloud projects list
```

**Option A: Use an existing project**
```bash
gcloud config set project YOUR_EXISTING_PROJECT_ID
```

**Option B: Create a new project (recommended for isolation)**
```bash
# Project IDs must be globally unique and lowercase
gcloud projects create moleculeai-api --name="MoleculeAI API"

# Set it as your active project
gcloud config set project moleculeai-api
```

**Verify your active project:**
```bash
gcloud config get-value project
```

### 3. Enable Required APIs

```bash
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com
```

### 4. Set Up Billing (Required, but Free Tier Applied)

Even with free tier, Google Cloud requires billing enabled (you won't be charged if you stay under free limits):

1. Go to: https://console.cloud.google.com/billing
2. Create billing account (or use existing)
3. Link to your project

### 5. Build and Deploy

From the project root directory:

```bash
# Build and deploy to Cloud Run
gcloud run deploy moleculeai-api \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10 \
  --set-env-vars OPENAI_API_KEY=$OPENAI_API_KEY
```

**Note:** Replace `$OPENAI_API_KEY` with your actual key, or set it in Cloud Run console after deployment.

### 6. Get Your API URL

After deployment, you'll get a URL like:
```
https://moleculeai-api-xxxxx-uc.a.run.app
```

Save this URL - you'll need it for your Next.js frontend!

---

## 🔧 Configuration Options

### Environment Variables

Set in Cloud Run console or via CLI:

```bash
gcloud run services update moleculeai-api \
  --update-env-vars OPENAI_API_KEY=sk-...
```

### Resource Limits

- **Memory:** 2GB recommended (RDKit needs memory)
- **CPU:** 2 vCPU recommended (faster processing)
- **Timeout:** 300 seconds (5 minutes max)
- **Concurrency:** 10 requests per instance

### Cost Optimization

- **Min instances:** 0 (scale to zero when idle = no cost)
- **Max instances:** 10 (auto-scales under load)
- **CPU:** Only allocated during requests (pay per use)

---

## 🧪 Test Locally (Before Deploying)

### Option 1: Test with Docker

```bash
# Build image
docker build -t moleculeai-api -f api/Dockerfile .

# Run container
docker run -p 8080:8080 \
  -e OPENAI_API_KEY=sk-... \
  moleculeai-api

# Test
curl http://localhost:8080/health
curl -X POST http://localhost:8080/analyze \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "mode": "both"}'
```

### Option 2: Test FastAPI Directly

```bash
cd api
pip install -r requirements.txt
python main.py

# Test in another terminal
curl http://localhost:8080/health
```

---

## 📊 Monitoring

After deployment:

1. **View logs:**
   ```bash
   gcloud run logs read moleculeai-api --limit 50
   ```

2. **View in console:**
   https://console.cloud.google.com/run

3. **Check metrics:**
   - Request count
   - Latency
   - Error rate
   - CPU/Memory usage

---

## 🔒 Security

1. **Add authentication** (optional, for production):
   ```bash
   gcloud run deploy moleculeai-api --no-allow-unauthenticated
   ```

2. **Restrict CORS** in `main.py`:
   ```python
   allow_origins=["https://antiquant.vercel.app"]
   ```

3. **Use Cloud Run's built-in HTTPS** (automatic)

---

## 💰 Estimated Costs

**With free tier:**
- 0-15,000 requests/month: **$0**
- Under resource limits: **$0**

**If you exceed free tier:**
- Additional requests: $0.40 per million
- Additional vCPU: $0.00002400 per vCPU-second
- Additional memory: $0.00000250 per GB-second

**Typical monthly cost for small app: $0-5**

---

## 🔄 Update Deployment

To update after code changes:

```bash
gcloud run deploy moleculeai-api --source .
```

Or update specific settings:

```bash
gcloud run services update moleculeai-api \
  --update-env-vars KEY=value \
  --update-memory 4Gi
```

---

## 🐛 Troubleshooting

### Build fails
- Check `api/Dockerfile` and `api/requirements.txt`
- Ensure all dependencies are listed

### Deployment fails
- Check billing is enabled
- Ensure APIs are enabled (run, cloudbuild)
- Check quota limits in console

### Runtime errors
- Check logs: `gcloud run logs read moleculeai-api`
- Verify environment variables are set
- Check model files exist in `/models` directory

### Slow cold starts
- Set min instances to 1 (costs more but faster)
- Optimize Docker image size
- Use Cloud Run's concurrency settings

---

## 📝 Next Steps

After deployment:

1. **Save your Cloud Run URL**
2. **Update Next.js API routes** to call Cloud Run instead of Python spawn
3. **Test the integration**
4. **Monitor usage** in Cloud Console

See `CLOUD_RUN_INTEGRATION.md` for updating Next.js routes.

