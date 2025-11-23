# 🚀 Quick Start - Deploy to Cloud Run

## Step 1: Login & Set Up Project

```bash
# 1. Authenticate (opens browser)
gcloud auth login

# 2. Check existing projects (optional)
gcloud projects list

# 3a. Use existing project (if you have one)
gcloud config set project YOUR_PROJECT_ID

# OR 3b. Create new project (recommended)
gcloud projects create moleculeai-api --name="MoleculeAI API"
gcloud config set project moleculeai-api

# 4. Verify active project
gcloud config get-value project
```

## Step 2: Enable APIs & Billing

```bash
# Enable required APIs
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com

# Set up billing (go to browser)
# https://console.cloud.google.com/billing
# Link billing account to your project
```

## Step 3: Deploy!

```bash
# From project root directory
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

**Copy the URL you get!** (e.g., `https://moleculeai-api-xxxxx-uc.a.run.app`)

## Step 4: Update Vercel

1. Go to **Vercel Dashboard** → Your Project → **Settings** → **Environment Variables**
2. Add: `MOLECULEAI_API_URL=https://your-url-from-step-3`
3. Redeploy your Next.js app

## Done! 🎉

See `cloud-run-deploy.md` for detailed troubleshooting and configuration.

