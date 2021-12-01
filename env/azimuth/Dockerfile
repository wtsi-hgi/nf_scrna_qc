FROM satijalab/seurat:4.0.1

RUN apt-get update
RUN apt-get install -y libv8-dev

RUN mkdir lzf
WORKDIR /lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.h
RUN mkdir lzf
WORKDIR /lzf/lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_c.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_d.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzfP.h https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf.h
WORKDIR /lzf
RUN gcc -O2 -fPIC -shared lzf/*.c lzf_filter.c -I /usr/include/hdf5/serial/ -lhdf5_serial -o liblzf_filter.so
WORKDIR /
ENV HDF5_PLUGIN_PATH=/lzf

RUN R --no-echo -e "install.packages(c('ggplot2', 'hdf5r'), repo='https://cloud.r-project.org')"
RUN R --no-echo -e "remotes::install_github(c('mojaveazure/seurat-disk'), dependencies = FALSE)"

RUN mkdir /opt/PBMC_reference
ADD https://zenodo.org/record/4546839/files/idx.annoy https://zenodo.org/record/4546839/files/ref.Rds /opt/PBMC_reference/
RUN chmod 755 -R /opt/PBMC_reference
