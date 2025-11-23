# Use Python 3.11 slim image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies for RDKit and image processing
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libgif-dev \
    librsvg2-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better Docker layer caching
COPY api/requirements.txt .

# Upgrade pip first, then install Python dependencies
RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy application code (models, scripts, etc.)
COPY . .

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PORT=8080
ENV PYTHONPATH=/app

# Expose port
EXPOSE 8080

# Run the application
CMD ["python", "-m", "uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8080"]

