FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN pip install --upgrade pip \
    && pip install biopython --upgrade \
    && pip install wheel \
    && pip install pandas \
    && pip install openpyxl

RUN apt-get update \
    && apt-get -y install wget

ENV BLAST_VERSION='2.12.0'

RUN cd /opt \
    && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && tar zxvpf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

ENV PATH $PATH:/opt/ncbi-blast-${BLAST_VERSION}+/bin
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
