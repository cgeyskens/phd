# # Use the official Miniconda3 image as the base
# FROM continuumio/miniconda3:24.5.0-0

# # copy the conda environment file
# COPY synapse-counting-env-2.yml /opt

# # initialize Conda
# RUN conda init bash

# # create the environment and install required packages
# RUN conda env create -f /opt/synapse-counting-env-2.yml
# RUN echo "conda activate synapse-counting-2" >> ~/.bashrc 
# RUN /bin/bash -c "source ~/.bashrc"

# # debugging statement
# RUN echo "Contents of ~/.bashrc:" && cat ~/.bashrc

# ENV PATH="/opt/conda/bin:/opt/condabin:${PATH}"

# # # start bash shell with activated Conda environment
# # CMD ["conda", "run", "-n", "synapse-counting-2", "/bin/bash"]





# # Use the official Miniconda3 image as the base
# FROM continuumio/miniconda3:24.5.0-0

# # copy the conda environment file
# COPY synapse-counting-env-2.yml /opt

# # initialize Conda
# RUN conda init bash

# # create the environment and install required packages
# RUN conda env create -f /opt/synapse-counting-env-2.yml
# RUN echo "conda activate synapse-counting-2" >> ~/.bashrc 
# RUN /bin/bash -c "source ~/.bashrc"

# # Ensure the environment is activated before running any command
# ENV PATH /opt/conda/envs/synapse-counting-2/bin:$PATH

# # Start bash shell
# CMD ["conda", "run", "-n", "synapse-counting-2", "/bin/bash"]




# Use the official Miniconda3 image as the base
FROM continuumio/miniconda3:24.5.0-0

# Install Bash (version 3.x or later) and ps
RUN apt-get update && \
    apt-get install -y bash procps && \
    rm -rf /var/lib/apt/lists/*

# Copy the conda environment file
COPY synapse-counting-env-2.yml /opt

# Initialize Conda
RUN conda init bash

# Create the environment and install required packages
RUN conda env create -f /opt/synapse-counting-env-2.yml
RUN echo "conda activate synapse-counting-2" >> ~/.bashrc 

# Ensure the environment is activated before running any command
ENV PATH="/opt/conda/bin:/opt/condabin:${PATH}"

# Set the entrypoint to Bash
ENTRYPOINT ["/bin/bash", "-l", "-c"]

# Command to start the container
CMD ["-c", "source ~/.bashrc && exec /bin/bash"]

