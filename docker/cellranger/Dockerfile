FROM registry.access.redhat.com/rhel7:7.4-120

WORKDIR /opt

# download cellranger (note that download links expire)
RUN curl -o cellranger.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.1.tar.gz?Expires=1546659211&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4xLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU0NjY1OTIxMX19fV19&Signature=Yw-pMnpiKoSkBcayYD37SctxgFCOWn4s~vcB9MBfh14MSVqcHV0~aiNFfSvMYc1fTgsGYh0FSbpIZFIljeunKtgdYKinrPojcO43LvaXNo9nP4TI0xRxRslsy7IawFsIrnlEW6t4L4kfgPrp0DrHEgIJismcuYPEhoGs60PaAcDNw2rJidTfLgOdHqGRt9tmOAT5LT4ErwJXAN4AJKOSJWGmZYF4ARhQazlIQQNcKwBK6gXavGVIuPizPX-h0pUziQ6QlSPKfFXTiaO-rdWRkmz-RxUofzVpzWNgLzsq~x15SuVZ0pb9Kf7cJK9C3-8zfjWkHSFU4POugL4fg4SfTQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

RUN echo "5e1a8b670091bf374655bd7ec2ccafd1 cellranger.tar.gz" > cellranger.md5sum

RUN md5sum -c cellranger.md5sum

RUN tar -xzvf cellranger.tar.gz && mv cellranger-3.0.1 cellranger

ENV PATH="/opt/cellranger:${PATH}"