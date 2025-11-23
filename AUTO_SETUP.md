# 🚀 Automated Setup - I'll Do It For You!

I'm setting up everything I can automatically. Here's what I'm doing and what you need to do in your browser.

---

## ✅ What I'm Doing Automatically

1. ✅ Verifying project structure
2. ✅ Checking Cloud Run configuration  
3. ⏳ Enabling Cloud Run APIs (needs billing first)
4. ⏳ Deploying to Cloud Run (after billing enabled)
5. ✅ Creating final checklist for you

---

## 🌐 What You Need to Do (2 Browser Tasks)

### TASK 1: Enable Billing for Cloud Run (2 minutes)

**Required:** Billing must be enabled to use Cloud Run (free tier applies, you won't be charged under limits!)

1. **Open:** https://console.cloud.google.com/billing?project=moleculeai-api
2. **Click:** "Link a billing account" (or use existing)
3. **Follow prompts** to link billing
4. **Done!** (You won't be charged - free tier covers your usage)

**Then come back here** - I'll continue deploying!

---

### TASK 2: Fix Vercel Build Settings (3 minutes)

After billing is enabled and Cloud Run deploys, do this:

1. **Go to:** https://vercel.com/dashboard → `antiquant2` → **Settings** → **General**
2. **Verify Root Directory:** Should say `web` (if not, change it)
3. **Expand:** "Build and Output Settings"
4. **Clear these fields** (make them EMPTY):
   - Install Command: **DELETE everything, leave EMPTY**
   - Build Command: **DELETE everything, leave EMPTY**
   - Output Directory: **DELETE everything, leave EMPTY**
5. **Click:** "Save"
6. **Go to:** "Deployments" tab
7. **Click:** "Redeploy"

**That's it for Vercel!**

---

## 📋 After You Enable Billing

Once billing is enabled, I'll:
1. Enable Cloud Run APIs
2. Deploy your Python backend
3. Get your Cloud Run URL
4. Tell you what environment variables to add

**Come back here after enabling billing and I'll continue!**

---

## 🎯 Final Checklist (After Everything)

- [ ] Billing enabled (Task 1)
- [ ] Cloud Run deployed (I'll do this)
- [ ] Vercel build fixed (Task 2)
- [ ] Environment variables set (I'll tell you what to add)
- [ ] Site working! 🎉

---

**Start with Task 1 (Enable Billing), then let me know and I'll continue!**

