export LD_LIBRARY_PATH=/home/dhx/anaconda3/lib:$LD_LIBRARY_PATH
../../bin/SPONGE -i infiles/mdin 
python ../../tools/dat2nc.py -n 3 -x mdcrd -xyz 20 20 20 -nc 1.nc
#vmd infiles/water.parm7 1.nc
