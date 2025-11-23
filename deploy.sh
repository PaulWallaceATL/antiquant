#!/bin/bash
# Automated deployment script for Cloud Run

set -e  # Exit on error

echo "🚀 Starting Cloud Run deployment..."
echo ""

# Check if project is set
PROJECT=$(gcloud config get-value project 2>/dev/null || echo "")
if [ -z "$PROJECT" ] || [ "$PROJECT" != "moleculeai-api" ]; then
    echo "❌ Project not set correctly. Setting to moleculeai-api..."
    gcloud config set project moleculeai-api
fi

# Check billing
echo "📋 Checking billing status..."
BILLING=$(gcloud beta billing projects describe moleculeai-api --format="value(billingAccountName)" 2>/dev/null || echo "")

if [ -z "$BILLING" ] || [ "$BILLING" == "" ]; then
    echo ""
    echo "⚠️  BILLING NOT ENABLED!"
    echo ""
    echo "Please enable billing first:"
    echo "👉 https://console.cloud.google.com/billing?project=moleculeai-api"
    echo ""
    echo "After enabling billing, run this script again:"
    echo "   bash deploy.sh"
    exit 1
fi

echo "✅ Billing enabled: $BILLING"
echo ""

# Enable APIs
echo "🔧 Enabling Cloud Run APIs..."
gcloud services enable run.googleapis.com cloudbuild.googleapis.com --project=moleculeai-api
echo "✅ APIs enabled"
echo ""

# Check for OpenAI API key
if [ -z "$OPENAI_API_KEY" ]; then
    echo "⚠️  OPENAI_API_KEY not found in environment"
    echo "   You can set it now, or add it after deployment in Cloud Run console"
    echo ""
    read -p "Enter your OpenAI API key (or press Enter to skip): " OPENAI_KEY
    if [ ! -z "$OPENAI_KEY" ]; then
        export OPENAI_API_KEY="$OPENAI_KEY"
        echo "✅ API key set (will be added to Cloud Run)"
    fi
fi

# Deploy
echo "☁️  Deploying to Cloud Run..."
echo "   This will take 5-10 minutes..."
echo ""

DEPLOY_CMD="gcloud run deploy moleculeai-api \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10"

if [ ! -z "$OPENAI_API_KEY" ]; then
    DEPLOY_CMD="$DEPLOY_CMD --set-env-vars OPENAI_API_KEY=$OPENAI_API_KEY"
fi

eval $DEPLOY_CMD

echo ""
echo "✅ Deployment complete!"
echo ""
echo "📋 Next steps:"
echo "1. Copy the Service URL shown above"
echo "2. Go to Vercel Dashboard → antiquant2 → Settings → Environment Variables"
echo "3. Add: MOLECULEAI_API_URL = (your Cloud Run URL)"
echo "4. Add: OPENAI_API_KEY = (if not already set)"
echo "5. Redeploy Vercel"
echo ""

