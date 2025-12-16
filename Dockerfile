# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install necessary dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-wheel

# Upgrade pip to the latest version
RUN pip3 install --upgrade pip

# Get current version (latest as of July 2023) of Bertini
RUN wget https://bertini.nd.edu/BertiniLinux64_v1.6.tar.gz
RUN tar xzf BertiniLinux64_v1.6.tar.gz

# get CRIU
RUN wget https://mirrorcache-us.opensuse.org/repositories/devel:/tools:/criu/xUbuntu_20.04/amd64/criu_3.17.1-1_amd64.deb
RUN dpkg -i criu_3.17.1-1_amd64.deb; exit 0
RUN apt --fix-broken -y install

# Create a working directory for the application
WORKDIR /home/baccala/

RUN mkdir /home/baccala/helium

# Copy the current.py script into the working directory
# COPY helium-16.6-RQQ-Bertini-input /home/baccala/
# COPY input /home/baccala/

# Run the current.py script as the container command
# CMD ["python3", "current.py", "--world-size=2", "--rank=1"]
CMD ["/bin/bash", "/home/baccala/helium/bertini/bertini.sh"]
