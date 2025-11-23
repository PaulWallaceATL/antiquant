# Cloud Run Integration Guide

After deploying your Python API to Google Cloud Run, update your Next.js frontend to call it.

---

## 🔧 Update Next.js API Routes

### 1. Add Environment Variable

In **Vercel Dashboard** → Your Project → Settings → Environment Variables:

```
MOLECULEAI_API_URL=https://your-api-url-xxxxx-uc.a.run.app
```

**Or** add to `.env.local` for local development:
```
MOLECULEAI_API_URL=https://your-api-url-xxxxx-uc.a.run.app
```

---

### 2. Update `web/lib/molecular_inference.ts`

Replace the Python spawn implementation with HTTP calls:

```typescript
const API_URL = process.env.MOLECULEAI_API_URL || process.env.NEXT_PUBLIC_MOLECULEAI_API_URL || '';

async function runInference(smiles: string, mode: 'classical' | 'quantum' | 'both'): Promise<MolecularAnalysisResult> {
    if (!API_URL) {
        throw new Error('MOLECULEAI_API_URL environment variable is not set');
    }

    try {
        const response = await fetch(`${API_URL}/analyze`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ smiles, mode }),
        });

        if (!response.ok) {
            const errorData = await response.json().catch(() => ({ error: `HTTP ${response.status}` }));
            throw new Error(errorData.error || `HTTP ${response.status}`);
        }

        const result = await response.json();
        return result;
    } catch (error: any) {
        console.error('API call failed:', error);
        throw new Error(`Failed to analyze molecule: ${error.message}`);
    }
}
```

---

### 3. Update API Routes

The routes in `web/app/api/analyze/route.ts`, `web/app/api/analyze-quantum/route.ts`, and `web/app/api/compare/route.ts` should already work since they use the functions from `molecular_inference.ts`.

Just ensure the environment variable is set!

---

## 🧪 Testing

### Test Locally

1. **Deploy to Cloud Run first** (get the URL)
2. **Set environment variable:**
   ```bash
   export MOLECULEAI_API_URL=https://your-api-url-uc.a.run.app
   ```
3. **Run Next.js locally:**
   ```bash
   cd web
   npm run dev
   ```
4. **Test in browser:** http://localhost:3000

### Test in Production

1. **Set environment variable in Vercel**
2. **Redeploy** (or it will auto-deploy on next push)
3. **Test on your Vercel domain**

---

## 🔒 Security Considerations

### 1. Restrict CORS in Cloud Run API

Update `api/main.py`:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://antiquant.vercel.app",
        "https://antiquant-git-*-pauls-projects-*.vercel.app",  # Preview deployments
    ],
    allow_credentials=True,
    allow_methods=["POST", "GET"],
    allow_headers=["Content-Type"],
)
```

### 2. Add Rate Limiting (Optional)

Use Cloud Run's built-in rate limiting or add middleware in FastAPI.

### 3. Add Authentication (Optional)

For production, you might want to add API keys:

```python
from fastapi import Security, HTTPException
from fastapi.security import APIKeyHeader

api_key_header = APIKeyHeader(name="X-API-Key")

async def verify_api_key(api_key: str = Security(api_key_header)):
    if api_key != os.environ.get("API_KEY"):
        raise HTTPException(status_code=403, detail="Invalid API Key")
    return api_key

@app.post("/analyze")
async def analyze(request: AnalyzeRequest, api_key: str = Security(verify_api_key)):
    ...
```

Then add to Next.js requests:
```typescript
headers: {
    'Content-Type': 'application/json',
    'X-API-Key': process.env.MOLECULEAI_API_KEY,
}
```

---

## 📊 Monitoring

### View Logs

```bash
# Cloud Run logs
gcloud run logs read moleculeai-api --limit 50

# Stream logs
gcloud run logs tail moleculeai-api
```

### Monitor in Console

https://console.cloud.google.com/run

Check:
- Request count
- Latency (p50, p95, p99)
- Error rate
- CPU/Memory usage
- Cost

---

## 💡 Tips

1. **Caching**: Consider caching results for identical SMILES
2. **Batch requests**: If processing multiple molecules, consider batching
3. **Error handling**: Add retry logic for transient failures
4. **Timeout**: Set appropriate timeouts (Cloud Run has 300s max)

---

## 🐛 Troubleshooting

### CORS errors
- Check CORS settings in `api/main.py`
- Verify your Vercel domain is in allowed origins

### Timeout errors
- Increase Cloud Run timeout (max 300s)
- Optimize Python code for faster execution
- Consider pre-warming instances

### Connection errors
- Verify API URL is correct
- Check environment variable is set
- Verify Cloud Run service is deployed and running

### 500 errors
- Check Cloud Run logs
- Verify model files exist
- Check environment variables (OPENAI_API_KEY)

---

## 🎯 Next Steps

1. ✅ Deploy to Cloud Run
2. ✅ Get API URL
3. ✅ Set environment variable in Vercel
4. ✅ Update `molecular_inference.ts`
5. ✅ Test locally
6. ✅ Deploy to Vercel
7. ✅ Test in production
8. ✅ Monitor and optimize

