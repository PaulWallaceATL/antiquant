# 🔴 CRITICAL: Fix Vercel Install Command

## ❌ The Problem:

Your Vercel settings have a **custom Install Command** that says `cd web && npm install`.

When Root Directory is set to `web`, Vercel is ALREADY in that directory, so `cd web` fails!

## ✅ THE FIX (Do This NOW):

Looking at your screenshot, I can see Root Directory is set to `web` ✅, but there must be a **custom Install Command** set somewhere.

### Step-by-Step:

1. **In the same page you're on** (Build and Deployment settings):
   - Scroll down below "Root Directory" section
   - Look for **"Build Command"** or **"Install Command"** sections
   - **You need to find where `cd web && npm install` is set**

2. **If you see these fields, make them EMPTY:**
   - **Install Command:** Delete everything, leave EMPTY
   - **Build Command:** Delete everything, leave EMPTY  
   - **Output Directory:** Delete everything, leave EMPTY

3. **Click "Save"** (should become active after making changes)

4. **Go to Deployments tab**

5. **Click "Redeploy"**

---

## 🔍 Where to Find It:

Since Root Directory is already set, the `cd web` command must be in:
- **Build and Deployment** → Scroll down to find **"Build Command"** section
- OR in **"Override"** settings
- OR it might be in an older deployment setting

**Look for ANY field that has `cd web` or `cd web &&` and DELETE IT!**

---

## ✅ After Fixing:

Once you clear those commands and redeploy, Vercel will:
- Use Root Directory `web` automatically ✅
- Run `npm install` directly (no `cd` needed) ✅
- Build successfully ✅

---

**Do this now, then we'll connect Cloud Run!**

