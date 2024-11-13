FROM debian:bookworm

ENV QUARTO_VERSION 1.5.53

ARG PUID
ARG PGID

# Install OS Dependencies
RUN apt-get update -y && \
    apt-get install curl wget libssl-dev libcurl4-openssl-dev libxml2-dev libfontconfig1-dev git texlive-full python3 pandoc -y

RUN update-ca-certificates -f

RUN mkdir -p /tmp/quarto && wget -O "/tmp/quarto/quarto.deb" https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb && \
    dpkg -i /tmp/quarto/quarto.deb && \
    rm -rf /tmp/quarto

# Setup jenkins uid/gid
RUN addgroup -gid ${PGID} jenkins
RUN adduser --disabled-password --gecos '' --ingroup jenkins --uid ${PUID} jenkins
RUN usermod -a -G sudo jenkins
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN chown -R jenkins:sudo /usr/local

USER jenkins
WORKDIR /home/jenkins
