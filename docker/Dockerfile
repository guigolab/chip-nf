FROM grape/mapping:gem-1.7.1

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "ChIP-seq analysis pipeline image"

ENV _ZERONE_VERSION 5af03a1

ENV PATH=/phantompeakqualtools/:$PATH

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
  r-base \
  git \ 
  libboost-dev \
  python-dev \
  bc \
  && apt-get clean \
  && apt-get purge \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install spp  
RUN git clone https://github.com/kundajelab/phantompeakqualtools \
  && cd phantompeakqualtools/ \
  && R -e 'install.packages( c("snow","snowfall","bitops","caTools"), repos="http://cran.us.r-project.org"); \
    source("http://bioconductor.org/biocLite.R"); \
    biocLite("Rsamtools",suppressUpdates=TRUE); \
    install.packages("./spp_1.14.tar.gz");'

# install Macs2
RUN pip install numpy \ 
    && pip install macs2

# get KentUtils
RUN curl -fsSLo /usr/local/bin/bedGraphToBigWig https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedGraphToBigWig \
  && chmod +x /usr/local/bin/bedGraphToBigWig
RUN curl -fsSLo /usr/local/bin/bedClip https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedClip \
  && chmod +x /usr/local/bin/bedClip 

# install bedtools
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz | tar xz \
  && cd bedtools2 \
  && make \
  && make install

# install java
RUN curl --header "Cookie: oraclelicense=accept-securebackup-cookie" -L http://download.oracle.com/otn-pub/java/jdk/8u121-b13/e9e7ea248e2c4826b92b3f075a80e441/jre-8u121-linux-x64.tar.gz | tar xz \
  && update-alternatives --install /usr/bin/java java /jre1.8.0_121/bin/java 100

# install picard
RUN curl -L -o /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.8.3/picard.jar \
  && chmod +x /usr/local/bin/picard.jar

# install zerone
RUN mkdir -p zerone && \
    curl -fsSL https://github.com/nanakiksc/zerone/archive/${_ZERONE_VERSION}.tar.gz | tar xz --strip-components 1 -C zerone && \
    make -C zerone && \
    mv zerone/zerone /usr/local/bin && \
    make clean -C zerone

