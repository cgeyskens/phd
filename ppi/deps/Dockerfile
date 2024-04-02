FROM continuumio/miniconda3:24.1.2-0

# Copy the Python file & environment file
COPY ppi-env.yml /opt

# make mounting directory
RUN mkdir /opt/data

# Initialize Conda
RUN conda init bash

# Create the environment and install required packages
RUN conda env create -f /opt/ppi-env.yml && \
    echo "conda activate ppi" >> ~/.bashrc && \
    /bin/bash -c "source ~/.bashrc"

# Debugging statement
RUN echo "Contents of ~/.bashrc:" && cat ~/.bashrc

# Start bash shell with activated Conda environment
CMD ["conda", "run", "-n", "test_env", "/bin/bash"]



# to build this container: docker build -t ppi-con .
# to run this container, see run-docker.sh