# ✅ Vercel Import Setup - Right Now

## On the "New Project" Page:

### 1. Change Root Directory:
- Click "Edit" next to Root Directory
- Change from `./` to: `web`
- Press Enter or click away

### 2. Expand "Build and Output Settings":
- Click the arrow `>` to expand
- Make sure all fields are EMPTY:
  - Install Command: **EMPTY**
  - Build Command: **EMPTY**
  - Output Directory: **EMPTY**

### 3. Add Environment Variables:

Click the arrow `>` next to "Environment Variables" to expand it.

**Add Variable 1:**
- Key: `OPENAI_API_KEY`
- Value: `sk-...` (your OpenAI API key)
- Environments: ✅ Production, ✅ Preview, ✅ Development

**Add Variable 2:**
- Key: `MOLECULEAI_API_URL`
- Value: `https://moleculeai-api-xxxxx-uc.a.run.app` 
  - **Note:** Get this after Cloud Run deploys (check console in ~5-10 minutes)
  - **OR** you can add it later in Settings → Environment Variables
- Environments: ✅ Production, ✅ Preview, ✅ Development

### 4. Click "Deploy"!

---

## ⏳ After Deployment:

1. **Wait for Cloud Run to finish** (5-10 minutes)
2. **Get Cloud Run URL:**
   - Go to: https://console.cloud.google.com/run?project=moleculeai-api
   - Click `moleculeai-api` service
   - Copy the URL shown at top
3. **If you didn't add `MOLECULEAI_API_URL` yet:**
   - Go to: Vercel Dashboard → `antiquant` → Settings → Environment Variables
   - Add: `MOLECULEAI_API_URL` = `https://your-cloud-run-url`
   - Redeploy

---

## ✅ That's It!

Once Root Directory is `web` and environment variables are set, everything should work!

