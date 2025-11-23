# 🔍 Find Vercel Install Command Settings

## Where to Find It:

Since you're on the "Build and Deployment" settings page, here's exactly where to look:

### Option 1: Check "Build Command" Section

1. **Scroll down** on the "Build and Deployment" page you're viewing
2. **Look for these sections** (they might be collapsed/expanded):

   **"Build Command"** section:
   - May say something like "Custom build command"
   - Should be **EMPTY** or just say `npm run build`
   - If it has `cd web &&` in it, **DELETE IT ALL**

   **"Install Command"** section:
   - Should be **EMPTY** or just say `npm install`
   - If it has `cd web && npm install`, **DELETE IT ALL**

   **"Output Directory"** section:
   - Should be **EMPTY**
   - Vercel auto-detects `.next`

### Option 2: Check Project Settings Override

1. On the same "Build and Deployment" page
2. Look for any section that says **"Override"** or **"Custom"**
3. These might have the `cd web` commands

### Option 3: Check Deployment Settings

1. Go to **"Deployments"** tab
2. Click the **latest deployment** (the one that failed)
3. Look for **"Deployment Settings"** section
4. Expand it - there might be overrides there

### Option 4: It Might Be in a Config File

Check if there's a `vercel.json` in your GitHub repo that has build commands. If so, we need to remove or update it.

---

## 🎯 What You're Looking For:

Any field that contains:
- `cd web && npm install`
- `cd web && npm run build`
- Any command that starts with `cd web`

**Delete everything in those fields and leave them EMPTY!**

---

## ✅ Alternative: Screenshot It

If you can't find it:
1. Take a screenshot of the ENTIRE "Build and Deployment" page
2. Scroll down and screenshot all sections
3. I'll tell you exactly which field to clear

---

**Most likely it's in a "Build Command" or "Install Command" section that's below "Root Directory" on the same page.**

