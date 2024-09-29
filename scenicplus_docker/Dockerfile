# Use Ubuntu as the base image
FROM ubuntu:20.04

# Install system dependencies including build tools
RUN apt-get update && apt-get install -y \
    wget \
    git \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    curl \
    git \
    ca-certificates \
    && apt-get clean

# Install Miniconda on x86 or ARM platforms
ENV MINICONDA_PATH="/root/miniconda3"
ENV PATH="${MINICONDA_PATH}/bin:${PATH}"
RUN arch=$(uname -m) && \
    # if [ "$arch" = "x86_64" ]; then \
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    # elif [ "$arch" = "aarch64" ]; then \
    # MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    # else \
    # echo "Unsupported architecture: $arch"; \
    # exit 1; \
    # fi && \
    wget $MINICONDA_URL -O miniconda.sh && \
    bash miniconda.sh -b -p $MINICONDA_PATH && \
    rm -f miniconda.sh

# Check conda version
RUN conda --version

# Create the conda environment
RUN conda create --name scenicplus python=3.11 -y

# Initiate conda
RUN conda init bash

# Clone the SCENIC+ repository
RUN git clone https://github.com/aertslab/scenicplus /app/scenicplus

# Switch to the SCENIC+ directory
WORKDIR /app/scenicplus

# Fetch all branches and output the list of branches
RUN git fetch --all && git branch -a | tee /dev/stderr && git checkout development || echo "Branch 'development' not found."

# Activate the conda environment and install SCENIC+
RUN /bin/bash -c "source activate scenicplus && pip install ."

# Command to run when the container starts
CMD ["/bin/bash"]