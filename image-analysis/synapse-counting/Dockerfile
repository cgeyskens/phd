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
    synapse-counting==0.1.3 \
    dask==2025.2.0 \
    distributed==2025.2.0 \
    optuna==4.3.0 \
    seaborn==0.13.2 \
    numpy==2.2.5 \
    pandas==2.2.3 \
    scipy==1.15.3 \
    scikit-image==0.25.2 \
    matplotlib==3.10.0 \
    pillow==11.1.0 \
    lxml==5.4.0 \
    czifile==2019.7.2.1

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["python"]
