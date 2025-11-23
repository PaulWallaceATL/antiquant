# ✅ Complete Setup Checklist

Follow these steps **in order** to get everything working.

---

## 🔧 STEP 1: Fix Vercel Build

### A. Check Root Directory Setting

1. **Go to:** Vercel Dashboard → `antiquant2` → **Settings** → **General**
2. **Find:** "Root Directory" section
3. **Verify:** It says `web` (not `./` or empty)
4. **If wrong:** Click "Edit", type `web`, click "Save"

### B. Remove Custom Build Commands

1. **Still in Settings → General**
2. **Expand:** "Build and Output Settings" section
3. **Check each field:**

   **Install Command:**
   - ❌ Should NOT have: `cd web && npm install`
   - ✅ Should be: **EMPTY** (or just `npm install`)
   - **Action:** Delete any `cd web &&` part, leave empty

   **Build Command:**
   - ❌ Should NOT have: `cd web && npm run build`
   - ✅ Should be: **EMPTY** (or just `npm run build`)
   - **Action:** Delete any `cd web &&` part, leave empty

   **Output Directory:**
   - ✅ Should be: **EMPTY** (Vercel auto-detects `.next`)

4. **Click "Save"**

### C. Redeploy

1. Go to **Deployments** tab
2. Click **"Redeploy"** button (top right)
3. Click **"Redeploy"** in the popup
4. Wait for build to complete

**✅ Success indicator:** Build should complete without "cd: web: No such file or directory" error

---

## ☁️ STEP 2: Deploy Python Backend to Cloud Run

### A. Verify You're Ready

Check your terminal - you should have already run:
```bash
✅ gcloud auth login
✅ gcloud config set project moleculeai-api
```

### B. Enable Billing & APIs

1. **Set up billing** (if not already):
   - Go to: https://console.cloud.google.com/billing
   - Link a billing account to `moleculeai-api` project
   - **Note:** Free tier applies, you won't be charged under limits

2. **Enable required APIs:**
   ```bash
   gcloud services enable run.googleapis.com
   gcloud services enable cloudbuild.googleapis.com
   ```

### C. Deploy to Cloud Run

**From your project root directory** (`/Users/paulwallace/Desktop/quant-code-complexity`):

```bash
# Deploy (this will take 5-10 minutes first time)
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

**⚠️ Replace `sk-your-key-here` with your actual OpenAI API key!**

### D. Save Your Cloud Run URL

After deployment completes, you'll see:
```
Service [moleculeai-api] revision [moleculeai-api-00001-xxx] has been deployed...
Service URL: https://moleculeai-api-xxxxx-uc.a.run.app
```

**📋 COPY THIS URL!** You'll need it for Step 3.

### E. Test Cloud Run API

```bash
# Test health endpoint
curl https://your-url-from-above/health

# Should return: {"status":"healthy"}
```

---

## 🔗 STEP 3: Connect Vercel to Cloud Run

### A. Set Environment Variables in Vercel

1. **Go to:** Vercel Dashboard → `antiquant2` → **Settings** → **Environment Variables**
2. **Add these variables:**

   **Variable 1:**
   - Name: `OPENAI_API_KEY`
   - Value: `sk-...` (your OpenAI API key)
   - Environments: ✅ Production, ✅ Preview, ✅ Development

   **Variable 2:**
   - Name: `MOLECULEAI_API_URL`
   - Value: `https://moleculeai-api-xxxxx-uc.a.run.app` (URL from Step 2D)
   - Environments: ✅ Production, ✅ Preview, ✅ Development

3. **Click "Save"** for each variable

### B. Redeploy Vercel

After adding environment variables:
1. Go to **Deployments** tab
2. Click **"Redeploy"** on the latest deployment
3. **Or** push a new commit to trigger auto-deploy

---

## ✅ STEP 4: Test Everything

### A. Test Vercel Site

1. Go to your Vercel deployment URL
2. Try analyzing a molecule (e.g., SMILES: `CCO` for ethanol)
3. Should see predictions from Cloud Run API!

### B. Check Cloud Run Logs (if issues)

```bash
# View recent logs
gcloud run logs read moleculeai-api --limit 50

# Stream logs in real-time
gcloud run logs tail moleculeai-api
```

### C. Check Vercel Logs

1. Go to Vercel Dashboard → `antiquant2` → **Logs**
2. Check for any API errors

---

## 🐛 Troubleshooting

### Build Still Failing?

1. **Double-check Root Directory is `web`** (not `./web` or `/web`)
2. **Make sure all build commands are empty** in Settings
3. **Delete the project and re-import:**
   - Delete `antiquant2` in Vercel
   - Re-import from GitHub
   - **Set Root Directory to `web` during import**

### Cloud Run Deployment Failing?

1. **Check billing is enabled:** https://console.cloud.google.com/billing
2. **Verify APIs are enabled:**
   ```bash
   gcloud services list --enabled
   ```
3. **Check build logs:**
   ```bash
   gcloud run services describe moleculeai-api
   ```

### API Not Working on Site?

1. **Verify environment variables are set** in Vercel
2. **Check Cloud Run is running:**
   ```bash
   gcloud run services list
   ```
3. **Test Cloud Run directly:**
   ```bash
   curl -X POST https://your-cloud-run-url/analyze \
     -H "Content-Type: application/json" \
     -d '{"smiles": "CCO", "mode": "both"}'
   ```

---

## 📊 Quick Status Check

**✅ Vercel:**
- [ ] Root Directory set to `web`
- [ ] Build commands are empty
- [ ] Build passes
- [ ] Environment variables set
- [ ] Site loads

**✅ Cloud Run:**
- [ ] Project created (`moleculeai-api`)
- [ ] Billing enabled
- [ ] APIs enabled
- [ ] Service deployed
- [ ] URL saved
- [ ] Health check works

**✅ Integration:**
- [ ] `MOLECULEAI_API_URL` set in Vercel
- [ ] `OPENAI_API_KEY` set in Vercel
- [ ] Site can call Cloud Run API
- [ ] Molecules analyze successfully

---

## 🎉 You're Done!

Once all checkboxes are ✅, your app should work like it does locally!

**Need help?** Check the logs in Vercel and Cloud Run to see what's failing.

