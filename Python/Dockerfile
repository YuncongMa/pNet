# Yuncong Ma, 2/13/204
# script to setup a new Conda environment for creating a docker image

# Use the official Miniconda base image
FROM continuumio/miniconda3

# Set environment variables
ENV ENV_NAME pnet

# Copy the conda environment file to the container
COPY environment.yml /tmp/environment.yml

# Set shell compatibility to conda and -c
SHELL ["/bin/bash", "-c", "conda"]

# Create a new Conda environment
RUN conda env create python=3.8 && \
    conda clean -afy

# Activate the Conda environment
RUN echo "source activate $ENV_NAME" > ~/.bashrc
ENV PATH /opt/conda/envs/$ENV_NAME/bin:$PATH

# Install packages using pip or conda
# Do not use requirement.txt which can only use pip install
RUN pip install numpy scipy scikit-learn pandas h5py nilearn
RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -c conda-forge nibabel
RUN conda install pytorch
RUN conda install mesalib vtk
RUN pip install matplotlib surfplot plotnine
RUN python -m pip install -U scikit-image

# Set the working directory in the container
WORKDIR /app

# Copy the local contents into the container
COPY . /app

# Specify the default command to run when the container starts
CMD ["python", "app.py"]