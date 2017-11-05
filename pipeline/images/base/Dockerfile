# Base image from root desc-pipe image
FROM eiffl/desc-pipe-base:1.0
MAINTAINER francois.lanusse@gmail.com

# Install the wlmassmap package
RUN bash -lc "pip install git+git://github.com/EiffL/gcr-catalogs.git"
COPY wlmassmap-0.0.1.tar.gz tmp
RUN bash -lc "pip install tmp/wlmassmap-0.0.1.tar.gz"

USER vagrant
