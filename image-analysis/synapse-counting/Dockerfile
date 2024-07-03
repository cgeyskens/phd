# Use the official Miniconda3 image as the base
FROM continuumio/miniconda3:24.1.2-0

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the Docker image
COPY synapse-counting-env.yml .

# Install the conda environment specified in the environment.yml file
RUN conda env create -f synapse-counting-env.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "my_environment", "/bin/bash", "-c"]

# Ensure the environment is activated when a container is started
ENTRYPOINT ["conda", "run", "-n", "my_environment", "/bin/bash", "-c"]

# Copy the application code into the Docker image
COPY . .

# Define the default command to run when the container starts
CMD ["your_command_here"]
