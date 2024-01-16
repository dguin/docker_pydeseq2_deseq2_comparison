# Docker inheritance
FROM rocker/r-base

RUN apt-get update && apt-get install -y \
  libcurl4-gnutls-dev \
  libssl-dev \
  libopenblas-dev \
  libz-dev

# Install required Bioconductor package
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("DESeq2")'
# Required for GLM
RUN R -e 'BiocManager::install("apeglm")'

# Install python
RUN apt-get update && \ 
  apt-get install -y --no-install-recommends \
  python3-full \
  python3-pip \
  python3-venv

WORKDIR /opt
RUN python3 -m venv .venv
# Install pydeseq2
ARG PYDESEQ2_VERSION=0.4.4
RUN echo "Install pydeseq2 version=${PYDESEQ2_VERSION}"
RUN /opt/.venv/bin/python3 -m pip install pydeseq2==${PYDESEQ2_VERSION}
# Install dependencies for compare_r_py_deseq2
RUN /opt/.venv/bin/python3 -m pip install pandas

COPY compare_r_py_deseq2.py .
COPY deseq2.R .
COPY __init__.py .
RUN mkdir data
RUN mkdir output

# To just run R DESeq2 and get results from R
# ENTRYPOINT [ "Rscript", "deseq2.R ", "data", "output"]

ENTRYPOINT [ "/opt/.venv/bin/python3", "-m", "compare_r_py_deseq2", "--deseq2_rscript_path=deseq2.R", "--data_directory_path=data", "--output_directory_path=output"]
