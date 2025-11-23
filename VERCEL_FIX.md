# 🔧 Fix Vercel Build Error

## ❌ Current Error

```
Error: No Next.js version detected. Make sure your package.json has "next" in either "dependencies" or "devDependencies". Also check your Root Directory setting matches the directory of your package.json file.
```

## ✅ Solution: Set Root Directory in Vercel Dashboard

**The Root Directory MUST be set to `web/` in Vercel dashboard!**

### Step-by-Step Fix:

1. **Go to Vercel Dashboard**
   - https://vercel.com/dashboard
   - Open your project: **"antiquant"**

2. **Open Settings**
   - Click **"Settings"** in the top navigation
   - Make sure you're on the **"General"** tab

3. **Find "Root Directory"**
   - Scroll down to the **"Root Directory"** section
   - Click **"Edit"** button

4. **Set Root Directory**
   - Enter: `web` (just type `web` without the slash)
   - Click **"Save"**

5. **Redeploy**
   - Go to **"Deployments"** tab
   - Click the three dots (⋯) on the latest deployment
   - Click **"Redeploy"**
   - Or just push a new commit to trigger auto-deploy

### Why This Fixes It:

- ✅ When Root Directory is `web/`, Vercel automatically:
  - Looks for `package.json` in `web/`
  - Detects Next.js from `web/package.json`
  - Runs all commands from `web/` directory
  - Builds from the correct location

- ❌ Without Root Directory set, Vercel:
  - Looks in the repository root
  - Can't find `web/package.json`
  - Can't detect Next.js
  - Build fails

### After Setting Root Directory:

Your `vercel.json` has been simplified to just:
```json
{
  "framework": "nextjs"
}
```

This is all you need when Root Directory is set correctly!

---

## 🧪 Test After Fix:

1. **Set Root Directory** to `web/` in Vercel dashboard
2. **Redeploy** or push a new commit
3. **Check build logs** - should now detect Next.js
4. **Visit your site** - should work!

---

## 📝 Additional Notes:

- The warnings about peer dependencies (React 18 vs 19, etc.) are just warnings and won't break the build
- Make sure **Environment Variables** are set:
  - `OPENAI_API_KEY` (for chat features)
  - `MOLECULEAI_API_URL` (after you deploy to Cloud Run)

---

## ✅ You're Done!

Once Root Directory is set to `web/`, the build should succeed! 🎉

