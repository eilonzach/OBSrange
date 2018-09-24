#!/bin/bash

user=jrussell1
pass=JoshuaCMEMS2018

# outdir=./dataset-duacs-nrt-global-merged-allsat-phy-l4 

days=($(seq 16 1 29))
year=2018
month=04

for day in ${days[@]}; do
	dayend=$[$day+6]
	wget -nv -m --no-host-directories --cut-dirs=2 --accept gz \
	  ftp://$user:$pass@ftp.sltac.cls.fr/Core/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/dataset-duacs-nrt-global-merged-allsat-phy-l4-v3/nrt_global_allsat_phy_l4_${year}${month}${day}_*.nc.gz
# 	mv ./nrt_global_allsat_phy_l4_${year}${month}${day}_${year}${month}${dayend}.nc $outdir
done
