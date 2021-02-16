ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION} as builder

# install needed tools
RUN apt-get update \
    &&  DEBIAN_FRONTEND=noninteractve apt-get install -y --no-install-recommends \
         ca-certificates \
         curl \
         gdebi-core \
         git \
         libboost-all-dev \
         libcurl4-gnutls-dev \
         libxml2-dev \
         lsb-release \
        #  ncurses-dev \
         python3-dev \
         python3-pip \
         python-is-python3

# install R
ARG R_VER=3.6.3
RUN export UBUNTU_VER=$(lsb_release -rs | tr -d '.') \
    && curl -O https://cdn.rstudio.com/r/ubuntu-${UBUNTU_VER}/pkgs/r-${R_VER}_1_amd64.deb \
    && DEBIAN_FRONTEND=noninteractve gdebi --non-interactive r-${R_VER}_1_amd64.deb

ENV PATH=/opt/R/${R_VER}/bin:$PATH

RUN echo 'options(repos = "https://cloud.r-project.org/")' > ~/.Rprofile

# install spp
ARG BIOC_VER="3.10"
RUN git clone https://github.com/kundajelab/phantompeakqualtools \
  && cd phantompeakqualtools/ \
  && R -e 'install.packages( c("snow","snowfall","bitops","caTools", "BiocManager")); \
    BiocManager::install(version = "'$BIOC_VER'"); \
    BiocManager::install("Rsamtools", ask=FALSE); \
    install.packages("spp");'

# install Macs2
# ARG MACS2_VER=2.1.1.20160309
ARG MACS2_VER=2.2.4
RUN pip3 install numpy \
    && pip3 install macs2==${MACS2_VER}

# get KentUtils
RUN curl -fsSLO https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedGraphToBigWig \
  && chmod +x bedGraphToBigWig
RUN curl -fsSLO https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedClip \
  && chmod +x bedClip

# install bedtools
ARG BEDTOOLS_VER=2.26.0
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VER}/bedtools-${BEDTOOLS_VER}.tar.gz | tar xz \
  && cd bedtools2 \
  && make

# install picard
ARG PICARD_VER=2.8.3
RUN curl -fsSLO https://github.com/broadinstitute/picard/releases/download/${PICARD_VER}/picard.jar \
  && chmod +x picard.jar

# install zerone
ARG ZERONE_VERSION=5af03a1
RUN mkdir -p zerone && \
    curl -fsSL https://github.com/nanakiksc/zerone/archive/${ZERONE_VERSION}.tar.gz | tar xz --strip-components 1 -C zerone && \
    make -C zerone

RUN apt-get install -y --no-install-recommends ncurses-dev

# install samtools
ARG SAMTOOLS_VER=1.11
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 \
    | tar xj && cd samtools-${SAMTOOLS_VER} \
    && make

# install GEM
ARG GEM_VER=3.6.1
RUN curl -fsSL https://github.com/smarco/gem3-mapper/archive/v${GEM_VER}.tar.gz \
    | tar xz && cd gem3-mapper-${GEM_VER} \
    && ./configure && make

RUN curl -O http://ppa.launchpad.net/linuxuprising/libpng12/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1.1+1~ppa0~focal_amd64.deb

FROM ubuntu:${UBUNTU_VERSION}

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "ChIP-seq analysis pipeline image"

# install needed tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gawk \
        libblas3 \
        libgomp1 \
        libicu66 \
        liblapack3 \
        libcurl3-gnutls \
        locales \
        python-is-python3 \
        openjdk-8-jre-headless

# set locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# set environment variables
ARG R_VER=3.6.3
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    PATH=/opt/R/${R_VER}/bin:$PATH

ARG SAMTOOLS_VER=1.11
ARG GEM_VER=3.6.1

COPY --from=builder /libpng12-0_1.2.54-1ubuntu1.1+1~ppa0~focal_amd64.deb /
RUN dpkg -i /libpng12-0_1.2.54-1ubuntu1.1+1~ppa0~focal_amd64.deb && \
    rm libpng12-0_1.2.54-1ubuntu1.1+1~ppa0~focal_amd64.deb

# Copy installed libraries from builder
COPY --from=builder /opt/R /opt/R
COPY --from=builder /usr/local/lib/python3.8/dist-packages /usr/local/lib/python3.8/dist-packages
COPY --from=builder /bedtools2/bin/ /usr/local/bin/
COPY --from=builder /bedClip /usr/local/bin
COPY --from=builder /bedGraphToBigWig /usr/local/bin/
COPY --from=builder /zerone/zerone /usr/local/bin/
COPY --from=builder /picard.jar /usr/local/bin/
COPY --from=builder /phantompeakqualtools/run_spp.R /usr/local/bin
COPY --from=builder /gem3-mapper-${GEM_VER}/bin/* /usr/local/bin/
COPY --from=builder /samtools-${SAMTOOLS_VER}/samtools /usr/local/bin/
COPY --from=builder /usr/local/bin/macs2 /usr/local/bin/
