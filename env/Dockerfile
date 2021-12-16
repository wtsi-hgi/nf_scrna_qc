FROM rocker/r-ver:4.0.3


############################ Setup: OS #########################################
# Set up locales properly
RUN apt-get update && \
    apt-get install --yes --no-install-recommends locales && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*




RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Use bash as default shell, rather than sh
ENV SHELL /bin/bash

# run critical system updates ... one could also use the artful-curl pack
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
        curl \
        git \
        zip \
        unzip \
        procps && \
    apt-get install --yes \
        libpng-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        zlib1g-dev \
        build-essential && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# set up a user, to avoid running containers as root
ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password \
    --gecos "Default user" \
    ${NB_USER}

# set up the user directory structure
RUN mkdir ${HOME}/tmp && \
    mkdir ${HOME}/bin && \
    mkdir ${HOME}/analysis
ENV PATH ${HOME}/bin:$PATH
################################################################################



############################ Setup: R packages #################################
USER root
WORKDIR ${HOME}/tmp
COPY Rpackages.R ${HOME}/tmp/Rpackages.R
RUN Rscript ${HOME}/tmp/Rpackages.R

############################ Setup: conda ######################################
# install conda
USER root
# ENV CONDA_FILE https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh #- scrublet doesnt work with python 3.9.5 hence we shouldnt use the latest version - mo11 
ENV CONDA_FILE https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh
RUN echo "export PATH=${HOME}/conda/bin:$PATH" > /etc/profile.d/conda.sh && \
    wget --quiet ${CONDA_FILE} -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${HOME}/conda && \
    rm /tmp/miniconda.sh && \
    ${HOME}/conda/bin/conda clean -atipy

ENV PATH ${HOME}/conda/bin:$PATH


# alternative conda install (instead of the manual install from above)
# FROM continuumio/miniconda3:latest
################################################################################


############################ Setup: non-conda tools ############################
# download things to tmp
WORKDIR ${HOME}/tmp

# Example: install GCTA
#RUN wget --quiet \
#    https://cnsgenomics.com/software/gcta/gcta_1.92.1beta5.zip && \
#    unzip gcta_1.92.1beta5.zip && \
#    chmod 755 gcta_1.92.1beta5/gcta64 && \
#    mv gcta_1.92.1beta5/gcta64 /bin/gcta64 && \
#    rm gcta_1.92.1beta5.zip && \
#    rm -r gcta_1.92.1beta5
################################################################################

############################ Setup: conda environment ##########################
# make sure conda files are owned by user
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}

# download the conda environment OR copy from local directory
COPY environment.yml ${HOME}/tmp/environment.yml

# install conda environment
# NOTE: conda clean -atipy
#   removes everything w/o confirmation (leaves all environments functional).

RUN conda config --set always_yes yes --set changeps1 no
# update conda:
RUN conda update -n base -c defaults conda
# use mamba solver to install packages
RUN conda install mamba -n base -c conda-forge
RUN mamba env create --prefix ${HOME}/conda/envs/cenv -f ${HOME}/tmp/environment.yml
RUN rm ${HOME}/tmp/environment.yml
# Set installed Conda env as default:
ENV CONDA_DEFAULT_ENV cenv
ENV PATH ${HOME}/conda/envs/cenv/bin:$PATH
RUN echo $PATH
RUN ls -ltra ${HOME}/conda/envs/cenv/bin

# clean-up  # USER root
RUN conda clean -atipy
RUN rm -rf /tmp/*



RUN python -V
# export the conda path to force container_env since .bashrc necessarily
# sourced when when the image is used
# ENV PATH ${HOME}/conda/envs/container_env/bin:$PATH
# ENV CONDA_DEFAULT_ENV container_env
# ENV CONDA_PREFIX ${HOME}/conda/envs/container_env

############################ Setup: R packages #################################
# R package xml2 does not install - amend configure file according to
# https://www.programmersought.com/article/12136751134/
# COPY xml2-config.patch ${HOME}/tmp/xml2-config.patch
# RUN wget https://cran.r-project.org/src/contrib/xml2_1.3.2.tar.gz \
#  && tar zxvf xml2_1.3.2.tar.gz \
#  && patch xml2/configure ${HOME}/tmp/xml2-config.patch \
#  && R CMD INSTALL ./xml2

# renv via R-packages
# RENV_CONFIG_SANDBOX_ENABLED = FALSE is to avoid automatic re-install of xml2 library
# after installing patched configuration of the xml2 library above
# COPY xml2-config.patch ${HOME}/tmp/xml2-config.patch
# RUN mkdir -p ${HOME}/tmp/renv/local && \
#    wget https://cran.r-project.org/src/contrib/xml2_1.3.2.tar.gz && \
#    tar zxvf xml2_1.3.2.tar.gz && \
#    patch xml2/configure ${HOME}/tmp/xml2-config.patch && \
#    tar zcvf ${HOME}/tmp/renv/local/xml2_1.3.2.tar.gz ./xml2 && \
#    echo 'source("renv/activate.R")' > ${HOME}/tmp/run.R && \
#    echo "Sys.setenv(RENV_PATHS_LOCAL = \"${HOME}/renv/local\")" >> ${HOME}/tmp/run.R && \
#    echo 'renv::install("xml2@1.3.2")' >> ${HOME}/tmp/run.R && \
#    echo 'renv::restore()' >> ${HOME}/tmp/run.R

# copy renv information from local directory
# COPY renv.lock ${HOME}/tmp/renv.lock
# COPY renv/activate.R ${HOME}/tmp/renv/activate.R

# RUN Rscript ${HOME}/tmp/run.R
################################################################################


############################ Setup: cleanup ####################################
# clear tmp if there is anything in there...
USER root
RUN rm -rf ${HOME}/tmp
RUN rm -rf /tmp/*
################################################################################


############################ Deployment ########################################
# set wd to user home
WORKDIR ${HOME}

# set bash as the default entry point
# NOTE: for Google Pipelines APIv2, docker does not execute if entrypoint is set
#ENTRYPOINT [ "/bin/bash", "-c" ]

# define default command as an alternative to ENTRYPOINT
# NOTE: for Google Pipelines APIv2, docker does not execute if CMD is set
#CMD ["/bin/bash"]

# set the user
# NOTE: for Snakemake + Kubernetes we get an error if this != root
# NOTE: for adding a monitoring script with cromwell, the monitoring script is
#       copied into the docker image and needs to be chmod to run. Permission
#       errors arise if not root.
USER ${NB_USER}
################################################################################
