# 🔍 Troubleshooting - Why Analysis Fails

## ❌ The Problem:

Your app at https://antiquant.vercel.app/ is loading, but **"Analysis failed"** error occurs.

This means:
1. **Either:** Cloud Run isn't deployed/running
2. **Or:** `MOLECULEAI_API_URL` environment variable isn't set in Vercel

---

## ✅ Fix Step-by-Step:

### STEP 1: Check Cloud Run Status

Go to: https://console.cloud.google.com/run?project=moleculeai-api

1. **Click** on `moleculeai-api` service
2. **Check:**
   - Is there a green checkmark? (service is running)
   - Is there a URL shown at the top?
   - Any errors in the status?

### STEP 2: Get Cloud Run URL

If the service is running:
- **Copy the URL** from the top of the page (looks like `https://moleculeai-api-xxxxx-uc.a.run.app`)

If it's NOT running or failed:
- Check the **"Revisions"** tab
- Look at the latest revision logs
- See what error occurred

### STEP 3: Add Environment Variable to Vercel

1. **Go to:** https://vercel.com/dashboard → `antiquant` → **Settings** → **Environment Variables**
2. **Check if `MOLECULEAI_API_URL` exists:**
   - If NOT, add it:
     - Key: `MOLECULEAI_API_URL`
     - Value: `https://your-cloud-run-url-from-step-2`
     - Environments: ✅ Production, ✅ Preview, ✅ Development
   - If YES, verify the URL is correct
3. **Click "Save"**
4. **Redeploy:**
   - Go to **Deployments** tab
   - Click **"Redeploy"** on the latest deployment

---

## 🐛 Common Issues:

### Issue 1: Cloud Run Not Running

**Solution:** Check Cloud Run console for errors. The container might be failing to start.

### Issue 2: Environment Variable Not Set

**Symptom:** App loads but API calls fail immediately
**Solution:** Add `MOLECULEAI_API_URL` to Vercel environment variables

### Issue 3: Wrong URL in Environment Variable

**Symptom:** API calls timeout or return 404
**Solution:** Double-check the Cloud Run URL is correct

### Issue 4: CORS Error

**Symptom:** Network errors in browser console
**Solution:** Cloud Run API should already have CORS enabled, but check if your Vercel domain is allowed

---

## 🧪 Test Cloud Run Directly:

Once you have the Cloud Run URL, test it:

```bash
curl -X POST https://your-cloud-run-url/health
# Should return: {"status":"healthy"}

curl -X POST https://your-cloud-run-url/analyze \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "mode": "both"}'
# Should return JSON with analysis results
```

If these fail, Cloud Run isn't working properly.

---

**Start with STEP 1 - Check Cloud Run console to see what's wrong!**

