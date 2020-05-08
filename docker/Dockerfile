FROM debian:jessie
MAINTAINER Edgar Garriga <edgano@gmail.com>

## Docker image for regressive multiple sequence aligment including all the requirements
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ed \
      less \
      vim-tiny \
      wget \
      git \
      python \
      build-essential \
      cmake \
      curl \
      libargtable2-0 \
      ca-certificates \
      libpng12-0 \
      python-biopython \
      python-numpy \ 
      ruby \ 
      python-setuptools \
      default-jdk \
      valgrind \
      libpng-dev 

# install argtable 2
RUN wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz && \
    tar -zxf argtable2-13.tar.gz && \
    cd argtable2-13 && \
    ./configure && \
    make && \
    make install && \
    rm /argtable2-13.tar.gz

# install CLUSTAL OMEGA
RUN wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz && \
    tar -zxf clustal-omega-1.2.4.tar.gz && \
    cd clustal-omega-1.2.4 && \
    sed -i '1157s/1024/100/' src/clustal/muscle_tree.c && \
    ./configure && \
    make && \
    make install && \
    rm /clustal-omega-1.2.4.tar.gz

#install CLUSTALO w1
RUN wget http://www.clustal.org/download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.83.linux.tar.gz && \
    tar -zxf clustalw1.83.linux.tar.gz && \
    cd clustalw1.83.linux && \
    cp clustalw /bin/. && \
    rm /clustalw1.83.linux.tar.gz

#install CLUSTAL w2
RUN wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz && \
    tar -zxf clustalw-2.1.tar.gz && \
    cd clustalw-2.1 && \
    ./configure && \
    make && \
    make install && \
    rm /clustalw-2.1.tar.gz

#install fasta35 (needed for mafft)
RUN wget https://faculty.virginia.edu/wrpearson/fasta/fasta33-35/fasta-35.4.12.tar.gz && \
    tar -zxf fasta-35.4.12.tar.gz && \
    cd fasta-35.4.12/src && \
    make -f ../make/Makefile.linux_sse2 && \
    cp ../bin/* /bin/. && \
    rm /fasta-35.4.12.tar.gz

# install MAFFT
RUN wget http://mafft.cbrc.jp/alignment/software/mafft-7.397-with-extensions-src.tgz && \
    tar xfvz mafft-7.397-with-extensions-src.tgz && \
    cd mafft-7.397-with-extensions/core/ && \
    sed -i "s/PREFIX = \/usr\/local/PREFIX = \/mafft/g" Makefile && \
    sed -i "s/BINDIR = \$(PREFIX)\/bin/BINDIR = \/mafft\/bin/g" Makefile && \
    make clean && \
    make && \
    make install && \
    wget http://mafft.cbrc.jp/alignment/software/newick2mafft.rb && \
    chmod +x newick2mafft.rb && \
    cp newick2mafft.rb /bin/. && \
    export "PATH=$PATH:/mafft/bin" && \ 
    cp /mafft/bin/* /bin/. && \
    mv /mafft-7.397-with-extensions /mafft && \
    rm /mafft-7.397-with-extensions-src.tgz

# Install PROBCONS
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    probcons

# Install MSAProbs
RUN wget "https://downloads.sourceforge.net/project/msaprobs/MSAProbs-0.9.7.tar.gz" -O msaprobs.tar.gz && \
    tar zxf msaprobs.tar.gz && \
    cd MSAProbs-0.9.7/MSAProbs && \ 
    make && \
    cp msaprobs /usr/bin && \
    rm /msaprobs.tar.gz

# install UPP 
RUN git clone http://github.com/smirarab/sepp.git && \
    cd sepp && \
    mkdir -p /home/lib/python2.7/site-packages/ && \
    export PYTHONPATH=$PYTHONPATH:/home/lib/python2.7/site-packages/ && \
    python setup.py config -c && \
    echo "/home/" > /sepp/home.path && \
    sed -i "s/root/home/g" /sepp/.sepp/main.config && \
    python setup.py install --prefix=/home/ && \
    python setup.py develop && \
    mkdir /pasta-code && \
    cd /pasta-code && \
    git clone https://github.com/smirarab/pasta.git && \
    git clone https://github.com/smirarab/sate-tools-linux.git && \
    cd pasta && \
    python setup.py develop -user && \ 
    export PATH=$PATH:/pasta-code/pasta:/sepp  && \
    cd  /sepp && \
    python setup.py upp -c && \
    sed -i "s/root/home/g" /sepp/.sepp/upp.config && \
    cd /pasta-code/pasta && \ 
    python setup.py develop

# Install MSA
#RUN wget ftp://ftp.ncbi.nih.gov/pub/msa/msa.tar.Z && \
#    tar xfvz msa.tar.Z &&\
#    cd msa && \
#    make clean && \
#    make msa && \
#    rm /msa.tar.Z && \
#    chmod +x /msa/msa && \
#    cp /msa/msa /bin/.

# install FAMSA
RUN git clone https://github.com/refresh-bio/FAMSA.git && \
  cd FAMSA && \ 
  git checkout 1ca9cb2093b07b90a242a54741f5a4f9b39b888a && \ 
  sed -i 's/all: famsa famsa-gpu/all: famsa/' makefile && \
  sed -i 's/NO_AVX = false/NO_AVX = true/' makefile && \
  sed -i 's/NO_AVX2 = false/NO_AVX2 = true/' makefile && \
  sed -i 's/NO_GPU = false/NO_GPU = true/' makefile && \
  make && \
  cp famsa /bin/.
  
# Install TCOFFEE
RUN git clone https://github.com/cbcrg/tcoffee.git tcoffee && \
    cd tcoffee && \
    git checkout 88e58f9c8f4dd0fdaad6b0f32c3186565bfa3cc3 && \
    cd t_coffee/src && \
    make t_coffee && \
    cp t_coffee /bin/.
    #cp t_coffee *.pl /bin/. && \
    #mkdir /tmp/cache && mkdir /tmp/lock && \
    
ENV HOME '/'

ENV FASTA_4_MAFFT '/bin/fasta35'

ENV HOME_4_TCOFFEE '/tmp/.tcoffee/'
ENV CACHE_4_TCOFFEE '/tmp/.tcoffee/cache'
ENV LOCKDIR_4_TCOFFEE '/tmp/.tcoffee/lock'
ENV TMP_4_TCOFFEE '/tmp/.tcoffee/tmp'
ENV TRASH_4_TCOFFEE '/tmp/.tcoffee/trash'

#unset MAFFT_BINARIES
ENV MAFFT_BINARIES ''

#configure perl5 to BLAST templates
ENV PERL5LIB '/usr/share/perl5:/tcoffee/perl/lib/perl5'


# Install BLAST+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.10.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.10.0+/bin/* /usr/bin/.

# Install UNIREF50
#RUN wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz && \
#    gzip -d uniref50.fasta.gz