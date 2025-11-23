# ✅ FINAL STEPS - What I Need You To Do

I've automated everything I can. Here's what you need to do:

---

## 🌐 STEP 1: Enable Billing (2 minutes)

**This is required for Cloud Run (free tier applies - you won't be charged under limits!)**

1. **Click this link:** https://console.cloud.google.com/billing?project=moleculeai-api
2. **Click:** "Link a billing account" (or use existing)
3. **Follow the prompts**
4. **Done!**

**After enabling billing, run:**
```bash
bash deploy.sh
```

**OR** tell me "billing enabled" and I'll continue the deployment!

---

## 🔧 STEP 2: Fix Vercel Settings (3 minutes)

1. **Go to:** https://vercel.com/dashboard → `antiquant2` → **Settings** → **General**

2. **Check Root Directory:**
   - Should say: `web`
   - If not, click "Edit" and change to `web`

3. **Expand "Build and Output Settings":**
   - **Install Command:** Delete everything, leave **EMPTY**
   - **Build Command:** Delete everything, leave **EMPTY**  
   - **Output Directory:** Delete everything, leave **EMPTY**

4. **Click "Save"**

5. **Go to "Deployments" tab**

6. **Click "Redeploy"**

---

## ⏳ STEP 3: Wait for Cloud Run Deployment

After you enable billing and run `bash deploy.sh` (or tell me), I'll deploy Cloud Run.

**It takes 5-10 minutes.** You'll see the Service URL when it's done.

---

## 🔗 STEP 4: Connect Everything (2 minutes)

After Cloud Run deploys, you'll get a URL like:
```
https://moleculeai-api-xxxxx-uc.a.run.app
```

**Then:**

1. **Go to:** Vercel Dashboard → `antiquant2` → **Settings** → **Environment Variables**

2. **Add these variables:**
   
   **Variable 1:**
   - Name: `MOLECULEAI_API_URL`
   - Value: `https://your-cloud-run-url-from-above`
   - Check all environments: ✅ Production, ✅ Preview, ✅ Development
   
   **Variable 2** (if not already set):
   - Name: `OPENAI_API_KEY`
   - Value: `sk-...` (your OpenAI API key)
   - Check all environments: ✅ Production, ✅ Preview, ✅ Development

3. **Click "Save"** for each

4. **Redeploy Vercel** (Deployments tab → Redeploy)

---

## 🎉 DONE!

After all steps:
- ✅ Vercel build should pass
- ✅ Site should load
- ✅ Molecules should analyze using Cloud Run API

---

## 🐛 Troubleshooting

**If Vercel build still fails:**
- Double-check Root Directory is exactly `web` (not `./web`)
- Make sure ALL build commands are EMPTY
- Try deleting project and re-importing

**If Cloud Run deployment fails:**
- Make sure billing is enabled
- Check APIs are enabled: `gcloud services list --enabled`
- Check logs: `gcloud run logs read moleculeai-api`

---

**Start with STEP 1 (Enable Billing) and let me know when done!** 🚀

