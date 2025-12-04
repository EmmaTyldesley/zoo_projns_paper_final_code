#!/bin/bash
# remember - no spaces before or after equals sign!
INPATH="/Volumes/T7_EW_IPSL/RECICLE_IPSL-CM5A-MR/bodc/PML230161/1990-2039"
#OUTPATH="/Volumes/T7_EW_IPSL/RECICLE_IPSL-CM5A-MR/test_out"
OUTPATH="/Volumes/T7_EW_IPSL/RECICLE_IPSL-CM5A-MR/bodc/PML230161/1990-2039"

echo directory: $INPATH
echo --------

YEARSTART=2011
YEAREND=2019

#FILE="grid_T" # file type grid_T, bgc_T, grid_U,V,W

#loop over years and months
for (( YEAR=$YEARSTART; YEAR<=$YEAREND; YEAR++ ))
do
	#mkdir -p $OUTPATH/$YEAR
	for MONTH in $(seq -w 1 12)
	do
		for FILE in grid_T bgc_T grid_U grid_V grid_W
		do
			INFILES=$INPATH/$YEAR/AMM7_1d_$YEAR${MONTH}_${FILE}.nc
			echo input: $INFILES
			OUTFILES=$OUTPATH/$YEAR/AMM7_1m_$YEAR${MONTH}_${FILE}_from_daily.nc
			echo output: $OUTFILES
		
			ncra $INFILES $OUTFILES
			echo processed
			echo ---------
		done
	done
done