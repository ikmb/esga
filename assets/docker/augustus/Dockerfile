FROM ubuntu:20.04

# Set timezone in tzdata
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin"

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam 
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y samtools libhts-dev

# Install additional dependencies for homGeneMapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# Install additional dependencies for scripts
RUN apt-get install -y cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl git parallel
RUN apt-get install -y --no-install-recommends python3-biopython

# Clone AUGUSTUS repository

# Build AUGUSTUS
RUN mkdir -p /opt/augustus \
	&& cd /opt/augustus \
	&& git clone https://github.com/Gaius-Augustus/Augustus.git 3.4.0 \
	&& cd 3.4.0 \
	&& git checkout v3.4.0 \
	&& make clean \
	&& make \
	&& make install

ENV PATH="/opt/augustus/3.4.0/bin:/opt/augustus/3.4.0/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH=/opt/augustus/3.4.0/config"

# Test AUGUSTUS
RUN cd /opt/augustus/3.4.0 && make unit_test
