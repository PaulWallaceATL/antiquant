# 🔗 How to Get Your Cloud Run URL

## ✅ Quick Way - Cloud Console:

1. **Go to:** https://console.cloud.google.com/run?project=moleculeai-api

2. **Click** on the service name: `moleculeai-api`

3. **At the top of the page**, you'll see a section that says:
   - **Service URL:** `https://moleculeai-api-xxxxx-uc.a.run.app`
   - **OR** it might say "Creating revision..." if still deploying

4. **Copy the URL** shown (looks like: `https://moleculeai-api-xxxxx-uc.a.run.app`)

---

## 📋 What Value to Add in Vercel:

**In Vercel Environment Variables:**

1. **Key:** `MOLECULEAI_API_URL` (already exists, just needs value)

2. **Value:** `https://moleculeai-api-xxxxx-uc.a.run.app`
   - This is your Cloud Run URL from the console
   - Replace `xxxxx` with your actual service hash

3. **Environments:** ✅ Production, ✅ Preview, ✅ Development (check all)

4. **Click "Save"**

5. **Redeploy Vercel** (Deployments tab → Redeploy)

---

## ⏳ If Cloud Run is Still Deploying:

If you see "Creating revision..." in the Cloud Console:
- **Wait 5-10 minutes** for deployment to complete
- Then get the URL from the top of the page
- Add it to Vercel environment variables

---

## 🧪 Test the URL Works:

Once you have the URL, test it:

```bash
curl https://your-cloud-run-url/health
```

Should return: `{"status":"healthy"}`

If it returns an error or times out, Cloud Run isn't ready yet.

---

**Go to the Cloud Console link above and copy the URL you see there!**

