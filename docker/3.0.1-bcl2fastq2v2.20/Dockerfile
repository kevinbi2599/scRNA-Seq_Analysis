FROM gcr.io/vanallen-gcp-dfci/scrna/cellranger:3.0.1

WORKDIR /tmp

# install zip
RUN curl -o unzip.rpm http://mirror.centos.org/centos/7/os/x86_64/Packages/unzip-6.0-19.el7.x86_64.rpm \
    && rpm -i unzip.rpm

# download & install bcl2fastq
RUN curl -o bcl2fastq.zip https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip \
    && unzip bcl2fastq.zip \
    && rpm -i bcl2fastq*rpm
