FROM gcr.io/vanallen-gcp-dfci/scrna/cellranger:3.0.1-bcl2fastq2v2.20

# install lbzip2
RUN curl -o lbzip2.rpm http://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/l/lbzip2-2.5-1.el7.x86_64.rpm \
    && rpm -i lbzip2.rpm
