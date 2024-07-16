# Use the official Miniconda3 image as the base
FROM continuumio/miniconda3:latest

# Install Bash (version 3.x or later) and ps
RUN apt-get update && \
    apt-get install -y bash procps && \
    rm -rf /var/lib/apt/lists/*

# Copy the conda environment file
COPY synapse-counting-env.yml /opt

# Initialize Conda
RUN conda init bash

# Create the environment and install required packages
RUN conda env create -f /opt/synapse-counting-env.yml
RUN echo "conda activate synapse-counting" >> ~/.bashrc 

# Ensure the environment is activated before running any command
ENV PATH="/opt/conda/envs/synapse-counting/bin/:/opt/conda/bin:/opt/condabin:${PATH}"

# Set the entrypoint to Bash
# ENTRYPOINT ["/bin/bash", "-l", "-c"] THIS GAVE THE ERROR

# Command to start the container
CMD ["-c", "source ~/.bashrc && exec /bin/bash"]

