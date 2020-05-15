FROM ubuntu:19.10

ENV GUPPY_VERSION 3.6.0
ENV NANOPOLISHCOMP_VERSION 0.6.11
ENV PYCOQC_VERSION 2.5.0.21
ENV BEDPARSE_VERSION 0.2.3
ENV MINIMAP2_VERSION 2.17
ENV NANOPOLISH_VERSION 0.12.3
ENV HTSLIB_VERSION 1.10.2
ENV SAMTOOLS_VERSION 1.10
ENV NANOCOMPORE_VERSION 1.0.0rc3-1
ENV NANOCOMPORE_PIPELINE=true
ENV XDG_RUNTIME_DIR=""
ENV PYTHONNOUSERSITE="not_empty_to_avoid_loading_packages_from_user_home"

MAINTAINER Tommaso Leonardi <tommaso.leonardi@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update && \
	apt-get -y install \
		autoconf \
		python3 \
		python3-pip \
		python3-tk \
		python3-numpy \
		python3-scipy \
		git \
		wget \
		lsb-release \
		apt-transport-https \
		vim \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		libcurl4-openssl-dev \
		libidn11

RUN wget -O /guppy-cpu.tgz "https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz" \
	&& tar -zxf /guppy-cpu.tgz \
	&& rm -rf /guppy-cpu.tgz
ENV PATH="/ont-guppy-cpu/bin/:${PATH}"

# Install bedparse
RUN pip3 install bedparse==${BEDPARSE_VERSION}

# Install pycoQC
RUN pip3 install pycoQC==${PYCOQC_VERSION}

# Install NanopolishComp
RUN pip3 install nanopolishcomp=="${NANOPOLISHCOMP_VERSION}"

# Install bedtools
RUN apt-get -y install bedtools=2.27.1+dfsg-4

# Install Minimap2
RUN wget -O /minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" \
	&& tar -jxf /minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
	&& cp /minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin

# Install nanopolish
RUN git clone --recursive https://github.com/jts/nanopolish.git \
	&& cd /nanopolish \
	&& git checkout "v${NANOPOLISH_VERSION}" \
	&& make && cp /nanopolish/nanopolish /usr/local/bin

# HTSlib and Samtools
RUN git clone https://github.com/samtools/htslib.git \
	&& cd /htslib \
	&& git checkout ${HTSLIB_VERSION}

RUN git clone https://github.com/samtools/samtools.git \
	&& cd /samtools \
	&& git checkout ${SAMTOOLS_VERSION} \
	&& autoreconf \
	&& ./configure --without-curses --with-htslib=../htslib \
	&& make \
	&& make install

# Install nanocompore
RUN pip3 install pytest
RUN pip3 install nanocompore==${NANOCOMPORE_VERSION}

