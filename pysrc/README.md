# Driving DSMACC with Python

Python is a great way to control DSMACC.

There are three files to get things started

 - Makefile can be called (`make install` or `make test`)
            to install the python dsmacc module and optionally
            run box.py, trajectory.py and plot.py
 - setup.py can be called (`python setup.py install`),
            but depends on source code in ../cri/src
            and ../geoschem/src
 - box.py   is the most simple box model driver for dsmacc
            it sets initial conditons and some run options
            to run CRI (or geoschem) to make boxconc.dat
            and boxrate.dat output files
 - trajectory.py conceptually adds emissions and a dynamic
            PBL/TEMP to box.py to make trajectoryconc.dat
            and trajectoryrate.dat output files
 - plot.py  takes an input path (eg. trajectory.dat) and
            and output path to create a figure of the output
            concentrations

# Run Examples including install

## Abstract Steps

1. Install flex library (libfl.a)
2. Clone repository
3. Build and install
4. Execute box.py, trajectory.py

## Ubuntu

```
# Command may vary on non Ubuntu platforms
apt-get install libfl-dev
git clone https://github.com/barronh/dsmacc
cd dsmacc
# Library location may vary on non Ubuntu platforms
LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ CFLAGS=-static ./configure
cd pysrc
make install
python box.py
python trajectory.py
```
