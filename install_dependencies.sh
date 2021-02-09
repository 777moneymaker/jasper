#!/usr/bin/env bash

function fail {
    printf 'JASPER-FAIL: %s\n' "$1" >&2  ## Send message to stderr. Exclude >&2 if you don't want it that way.
    exit "${2-1}"  ## Return a code specified by $2 or 1 by default.
}


# Check for sudo
if [ "$EUID" -ne 0 ]
  then echo "JASPER-FAIL: Please run as root"
  exit
fi

# autoconf for make
sudo apt install -y cmake autoconf || fail "Unable to download and install cmake or autoconf"

# NCBI-Blast+ package
apt install -y ncbi-blast+ || fail "Unable to install NCBI-Blast+ or tRNAscan-se"
echo "JASPER: NCBI-Blast+ installed"

# wget for dependencies
sudo apt install -y wget

# Infernal download, compile and install
wget http://eddylab.org/infernal/infernal-1.1.4.tar.gz
tar -xvf infernal-1.1.4.tar.gz
rm -rf infernal-1.1.4.tar.gz
cd infernal-1.1.4/
./configure --prefix=/usr/local
make
make install
cd ..
rm -rf infernal-1.1.4/
echo "JASPER: Infernal installed"

# tRNAscan-SE download, compile and install
wget http://trna.ucsc.edu/software/trnascan-se-2.0.7.tar.gz
tar -xvf trnascan-se-2.0.7.tar.gz && rm trnascan-se-2.0.7.tar.gz
cd tRNAscan-SE-2.0/
./configure --prefix=/usr/local
make
make install
cd ..
rm -rf tRNAscan-SE-2.0
echo "JASPER: tRNAscan-SE installed"

# Piler-CR download and install
wget http://www.drive5.com/pilercr/pilercr1.06.tar.gz
tar -xvf pilercr1.06.tar.gz
rm pilercr1.06.tar.gz
cp pilercr1.06/pilercr /usr/local/bin/
rm -rf pilercr1.06/
echo "JASPER: Piler-CR installed"

# WIsH download, compile and install
git clone https://github.com/soedinglab/WIsH.git
cd WIsH
cmake .
make
cd ..
cp WIsH/WIsH /usr/local/bin
rm -rf WIsH
echo "JASPER: WIsH installed"

# Mash download, compile and install
wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar
tar -xvf mash-Linux64-v2.2.tar
rm mash-Linux64-v2.2.tar
cp mash-Linux64-v2.2/mash /usr/local/bin/
rm -rf mash-Linux64-v2.2/
echo "JASPER: Mash installed"


