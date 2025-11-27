# Use official Python 3.11 image as the base, ensuring the correct architecture
FROM python:3.11-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install system dependencies required for some Python packages
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install the required Python packages
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    dask \
    synapse-counting==0.1.2 \
    seaborn \
    optuna

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["python"]
