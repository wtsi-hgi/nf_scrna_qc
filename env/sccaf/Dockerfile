FROM rocker/r-ver:4.0.3


############################ Setup: OS #########################################
# Set up locales properly
RUN apt-get update && \
    apt-get install --yes --no-install-recommends locales python3.8-venv && \
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
    apt-get upgrade -y && \
    apt-get --assume-yes install rsync && \
    apt-get install -y git procps
    
# set a python virtual env
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# update pip
RUN /opt/venv/bin/python3 -m pip install --upgrade pip

# install sccaf in virtual env
RUN git clone https://github.com/wtsi-hgi/sccaf.git
RUN pip3 install --no-cache-dir ./sccaf
# current sccaf version fails on umap >= 0.5 due to change in interface
# cf. https://github.com/theislab/scanpy/issues/1579
# also requires to keep numba below 0.53, cf. https://github.com/basilkhuder/Seurat-to-RNA-Velocity/issues/8
RUN pip3 install --no-cache-dir umap-learn==0.4.6 numba==0.52.0 debugpy>=1.2.1

# drop in a python shell by default (on 'docker run -it container:tag')
CMD ["python"]

# HGI manual build and push to Dockerhub:
# nohup sudo docker build . -t wtsihgi/nf_qc_cluster:sccaf_1.5 &
# nohup docker push wtsihgi/nf_qc_cluster:sccaf_1.5 &
