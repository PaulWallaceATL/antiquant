# 🔧 Fix: "cd: web: No such file or directory"

## Problem

Vercel is trying to run `cd web && npm install` even though Root Directory is set to `web/`. When Root Directory is `web/`, Vercel should already be IN that directory.

## Solution: Check Project Settings

### Step 1: Verify Root Directory

1. Go to **Vercel Dashboard** → Your Project (`antiquant2`)
2. Click **Settings** → **General**
3. Scroll to **"Root Directory"**
4. Make sure it says: `web` (not `./` or empty)
5. If it's wrong, click **Edit** and set it to `web`
6. Click **Save**

### Step 2: Remove Custom Build Commands

1. Still in **Settings** → **General**
2. Scroll to **"Build and Output Settings"**
3. Expand it (click the arrow)
4. Check if there are any custom commands:
   - **Install Command** - should be empty or just `npm install`
   - **Build Command** - should be empty or just `npm run build`
   - **Output Directory** - should be empty (Vercel auto-detects `.next`)
5. **Remove any commands that include `cd web`**
6. Click **Save**

### Step 3: Redeploy

1. Go to **Deployments** tab
2. Click **Redeploy** on the latest deployment
3. Or push a new commit to trigger auto-deploy

---

## What Should Happen

When Root Directory is set to `web/`:

✅ Vercel automatically:
- Changes to `web/` directory
- Finds `web/package.json`
- Detects Next.js
- Runs `npm install` (no `cd web` needed)
- Runs `npm run build`
- Serves from `web/.next`

❌ You should NOT need:
- Custom install command: `cd web && npm install`
- Custom build command: `cd web && npm run build`
- Custom output directory

---

## Alternative: If Settings Don't Work

If the Root Directory setting isn't working, you can create a `vercel.json` in the **root** of your repo:

```json
{
  "buildCommand": "npm run build",
  "installCommand": "npm install",
  "outputDirectory": ".next",
  "framework": "nextjs"
}
```

But this should be in the `web/` directory, not the root. Actually, if Root Directory is set correctly, you don't need `vercel.json` at all!

---

## Quick Check

After setting Root Directory to `web/` and removing custom build commands:

1. **Redeploy**
2. **Check build logs** - should see:
   - `Running "install" command: npm install` (NOT `cd web && npm install`)
   - `Running "build" command: npm run build` (NOT `cd web && npm run build`)

---

## Still Not Working?

If it still fails:
1. Double-check Root Directory is `web` (not `./web` or `/web`)
2. Make sure there's no `vercel.json` in the root with build commands
3. Check that `web/package.json` exists and has Next.js
4. Try deleting the project and re-importing with Root Directory set from the start

