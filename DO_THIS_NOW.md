# 🚨 DO THIS NOW - Step by Step

## ✅ STEP 1: Fix Vercel Build (Do This First!)

**This will get your site working.**

1. **Go to:** https://vercel.com/dashboard
2. **Click:** `antiquant2` project
3. **Click:** **Settings** (top nav)
4. **Click:** **General** tab
5. **Scroll down to:** **"Build and Output Settings"**
6. **Click the arrow** to expand it
7. **For each field, make it EMPTY or remove `cd web &&`:**

   **Install Command:**
   - Delete everything, leave **EMPTY**
   
   **Build Command:**
   - Delete everything, leave **EMPTY**
   
   **Output Directory:**
   - Delete everything, leave **EMPTY**

8. **Scroll up to:** **"Root Directory"**
9. **Verify it says:** `web` (not `./`)
10. **Click:** **Save** button
11. **Go to:** **Deployments** tab
12. **Click:** **Redeploy** button
13. **Wait** for build to complete

**✅ If build succeeds, your site should work!** (It will use local Python for now, which might not work on Vercel, but the site will load)

---

## ☁️ STEP 2: Deploy Cloud Run (After Vercel Works)

**Only do this after Vercel build passes!**

### Enable APIs First

Run this command:
```bash
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com
```

Say **"y"** when it asks to enable APIs.

### Set Up Billing

1. **Go to:** https://console.cloud.google.com/billing?project=moleculeai-api
2. **Link billing account** to `moleculeai-api` project
3. **Note:** Free tier applies - you won't pay under limits!

### Deploy

Run this (replace `sk-...` with your real OpenAI key):
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
  --set-env-vars OPENAI_API_KEY=sk-your-actual-key-here
```

**⏰ This takes 5-10 minutes the first time!**

### Copy the URL

After it finishes, you'll see:
```
Service URL: https://moleculeai-api-xxxxx-uc.a.run.app
```

**📋 Copy that URL!**

---

## 🔗 STEP 3: Connect Everything

1. **Go to:** Vercel Dashboard → `antiquant2` → **Settings** → **Environment Variables**
2. **Add:**
   - `MOLECULEAI_API_URL` = `https://your-cloud-run-url-from-step-2`
   - `OPENAI_API_KEY` = `sk-...` (if not already set)
3. **Save**
4. **Redeploy** Vercel

---

## ✅ DONE!

After all steps, your app will work like it does locally!

**Need help?** Check `COMPLETE_SETUP_CHECKLIST.md` for detailed troubleshooting.

