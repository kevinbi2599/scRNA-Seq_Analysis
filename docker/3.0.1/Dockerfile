FROM centos:7

WORKDIR /opt

# download cellranger (note that download links expire)
RUN curl -o cellranger.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.1.tar.gz?Expires=1546678354&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4xLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU0NjY3ODM1NH19fV19&Signature=QNkk5l5zI7iXyIzl3cwiRDR2Nz7xRIFpDxjmzAYmFQ-TcUPo8EpRSnt8Iakx9sAy10iLvSQxprA3Lm0jIroVFgxC-MFZ-anlanLCUOgqRh~hdrVzsMgdrmgo~KwObWwFXQMT3lJ8sSTqAv0po~10xzASDV4qjFoFFkg5IR08SrVZZtgpi0gr-0E-6IAemc2WS4nJ6lTkDXhMc6eXI9ENEs347Y4MwhAmZc7FV1nEJF8gUALsF0DLf0FaxKt0O21y4hAB2fAn2h8sqc9sxYPj1sLJQdQ0Jlkqe2xNnLJ3JgmnktL83~7BYy8qlwttAWs6usrx5E9845H9A~Q~a6abaw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

RUN echo "5e1a8b670091bf374655bd7ec2ccafd1 cellranger.tar.gz" > cellranger.md5sum

RUN md5sum -c cellranger.md5sum

RUN tar --absolute-names -xzvf cellranger.tar.gz && mv cellranger-3.0.1 cellranger

RUN yum -y install which

ENV PATH="/opt/cellranger:${PATH}"
