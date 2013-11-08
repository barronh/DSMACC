curl -kLO https://github.com/barronh/DSMACC/archive/master.zip
unzip -u master.zip
cd DSMACC-master
make kpp/bin/kpp
make check
