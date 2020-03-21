FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-moff-ms1"


RUN apt-get -y install \
	apt-utils \
	apt-transport-https \
	gnupg


RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list


RUN apt-get update


RUN apt-get -y install \
	gcc \
	build-essential \
	libffi-dev \
	git \
	dirmngr \
	ca-certificates \
	mono-devel \
	mono-complete \
	ca-certificates-mono \
	unzip


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a


RUN conda install python=3.6
RUN conda install numpy pandas scikit-learn
RUN conda install -c conda-forge pynumpress
RUN conda install simplejson
RUN conda install -c anaconda cython
RUN conda install -c conda-forge -c bioconda moff


ENV PATH /opt/conda/envs/nf-core-moff-ms1/bin:$PATH


# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "nf-core-moff-ms1", "/bin/bash", "-c"]

RUN conda init bash

#RUN pip install ms-peak-picker pyteomics brain-isotopic-distribution pymzml
#RUN pip install ms_deisotope


COPY configuration_file.ini /
