FROM fedora:latest

#ENV HTTP_PROXY "http://10.31.255.65"
#ENV HTTPS_PROXY "https://10.31.255.65"
#ENV FTP_PROXY "ftp://10.31.255.65"

#COPY src /opt/psfex/src
#COPY autogen.sh /opt/psfex/autogen.sh
#COPY Makefile.am /opt/psfex/Makefile.am
#COPY psfex.spec.in /opt/psfex/psfex.spec.in
#COPY configure.ac /opt/psfex/configure.ac

RUN dnf -y update && dnf -y install git make automake autoconf libtool gcc kernel-devel @development-tools
RUN dnf install -y fftw-devel
RUN dnf install -y atlas-devel
RUN dnf install -y plplot-devel

COPY . /opt/psfex/

WORKDIR /opt/psfex
RUN ./autogen.sh
RUN ./configure
RUN make -j
RUN make install
CMD ["./psfex"]
