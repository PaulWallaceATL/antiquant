# ⚡ Quick Setup Guide

## 🎯 Right Now: Set Root Directory in Vercel

On the "New Project" page you're seeing:

1. **Find "Root Directory"** section
2. **Click "Edit"** button (next to `./`)
3. **Change to:** `web` (type `web` without the slash)
4. **Click "Deploy"** button

That's it! 🎉

---

## After Deployment

### 1. Set Environment Variables

Go to **Vercel Dashboard** → Your Project → **Settings** → **Environment Variables**

Add:
- `OPENAI_API_KEY` = `sk-...` (your OpenAI API key)
- `MOLECULEAI_API_URL` = `https://your-cloud-run-url-uc.a.run.app` (after you deploy to Cloud Run)

### 2. Deploy Python Backend to Cloud Run

Follow `api/cloud-run-deploy.md` to deploy your Python API.

Then come back and add the `MOLECULEAI_API_URL` environment variable.

---

## ✅ That's It!

Once Root Directory is set to `web/`, Vercel will automatically:
- Detect Next.js ✅
- Build correctly ✅
- Deploy successfully ✅

