# ✅ FINAL STEPS - Almost Done!

## 🎯 What To Do RIGHT NOW:

### STEP 1: Fix Vercel Build (2 minutes) 🔴

**The issue:** Vercel has a custom Install Command with `cd web && npm install`.

**Fix it:**
1. On the **Build and Deployment** page you're viewing:
   - Scroll down below "Root Directory"
   - Look for **"Build Command"** or **"Install Command"** sections
   - **DELETE everything** in those fields (make them EMPTY)
   - Click **"Save"**
2. Go to **"Deployments"** tab
3. Click **"Redeploy"**

**This will fix your Vercel build!**

---

### STEP 2: Get Cloud Run URL (1 minute)

Once Cloud Run deployment is complete:

1. **Go to:** https://console.cloud.google.com/run?project=moleculeai-api
2. **Click** on `moleculeai-api` service
3. **Copy the URL** shown at the top (looks like `https://moleculeai-api-xxxxx-uc.a.run.app`)

**OR check via CLI:**
```bash
gcloud run services list --project=moleculeai-api
```

---

### STEP 3: Connect Cloud Run to Vercel (2 minutes)

1. **Go to:** Vercel Dashboard → `antiquant2` → **Settings** → **Environment Variables**
2. **Add these variables:**

   **Variable 1:**
   - Name: `MOLECULEAI_API_URL`
   - Value: `https://your-cloud-run-url-from-step-2`
   - Environments: ✅ Production, ✅ Preview, ✅ Development

   **Variable 2** (if not already set):
   - Name: `OPENAI_API_KEY`
   - Value: `sk-...` (your OpenAI API key)
   - Environments: ✅ Production, ✅ Preview, ✅ Development

3. **Click "Save"** for each

4. **Redeploy Vercel** (will auto-deploy or click Redeploy)

---

### STEP 4: Test! 🎉

1. Go to your Vercel deployment URL
2. Try analyzing a molecule (e.g., SMILES: `CCO`)
3. Should work! ✅

---

## 📊 Checklist:

- [ ] Fix Vercel build settings (clear Install/Build commands)
- [ ] Redeploy Vercel successfully
- [ ] Get Cloud Run URL
- [ ] Add `MOLECULEAI_API_URL` to Vercel environment variables
- [ ] Add `OPENAI_API_KEY` to Vercel (if not already set)
- [ ] Test the app!

---

**Start with STEP 1 (Fix Vercel) - that's blocking everything right now!**

