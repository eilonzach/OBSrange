#!/bin/bash

locator_code_path=../Relocate_postProcess/locate+fileprint_ze.py
survey_file_dir=.

source activate python2

for i in *.txt
do 
	sta=`echo $i | sed s/.txt//g`
	if echo $i | grep -q "orrected" 
		then
			echo " "
		else
			echo "================================================="
			echo "================ RELOCATING $sta ================"
			echo "================================================="
			slat=`awk '$3 ~ /(Latitude)/ {print $4}' < $i`
			slon=`awk '$3 ~ /(Longitude)/ {print $4}' < $i`
			sdep=`awk '$1 ~ /^Depth/ {print $3}' < $i`
			python $locator_code_path Y 1500 $survey_file_dir/$i $slat $slon $sdep "$sta"_SIOcorrected.txt
			#mv _Corrected.txt $i_corrected.txt
	fi
done

