FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/esga-1.1/bin:/opt/bin:$PATH
ENV EVM_HOME /opt/conda/envs/esga-1.1/opt/evidencemodeler-1.1.1
RUN cd /opt/bin && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords && chmod +x faSomeRecords
RUN mkdir -p /opt/blast && cd /opt/blast && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz && tar -xvf ncbi-blast-2.2.31+-x64-linux.tar.gz && \
rm ncbi-blast-2.2.31+-x64-linux.tar.gz && mv ncbi-blast-2.2.31+ 2.2.31
