# Vercel Deployment Setup

## ⚠️ CRITICAL: Root Directory Configuration

Your Next.js app is in the `web/` directory, so Vercel MUST be configured with the correct root directory.

### Steps to Fix the 404 Error:

1. **Go to Vercel Dashboard**: https://vercel.com/dashboard
2. **Open your project**: "antiquant"
3. **Go to Settings**: Click "Settings" in the top navigation
4. **General Tab**: Should be selected by default
5. **Scroll to "Root Directory"**: Find this section
6. **Click "Edit"**: Next to Root Directory
7. **Enter**: `web` (just type `web` without the slash)
8. **Save**: Click "Save" button
9. **Redeploy**: 
   - Go to "Deployments" tab
   - Click the three dots (⋯) on the latest deployment
   - Click "Redeploy"
   - Confirm the redeploy

### Why This is Needed:

- Your Next.js app (`app/`, `components/`, etc.) is in the `web/` directory
- Vercel needs to know where to look for `package.json`, `next.config.ts`, and the app files
- Without this setting, Vercel looks at the repository root and can't find the Next.js app
- This causes the 404 error even though the build succeeds

### After Setting Root Directory:

Once the root directory is set to `web/`, Vercel will:
- ✅ Find your `web/package.json`
- ✅ Detect Next.js automatically
- ✅ Build from the correct directory
- ✅ Serve the app correctly

### Environment Variables Needed:

Make sure these are set in Vercel Dashboard → Settings → Environment Variables:

- **`OPENAI_API_KEY`** (required for chat functionality)
- **`GEMINI_API_KEY`** (optional, only if using Gemini chat)

### Note About Python Backend:

The Python scripts (`scripts/molecular_inference.py`) won't work on Vercel serverless functions. You'll need to:
- Deploy the Python backend separately (Google Cloud Run, AWS Lambda, Railway, etc.)
- Update the API routes to call the external Python API

The build will succeed, but the `/api/analyze`, `/api/analyze-quantum`, and `/api/compare` routes will fail at runtime until the Python backend is deployed separately.

