if [ "$#" -ne 2 ]
then 
	echo "arg 1 : prefix out directories and arg 2 : inImg"
	exit 0
fi
sextractor $2 -c sextractor/prepsfex.sex
#mv sextractor/outTest.cat .
/usr/local/bin/psfex outTest.cat -c default.psfex
outDir="$1psfexOut"
echo ${outDir}
if [ -d ${outDir} ] 
then
    echo "Directory ${outDir} already exists, files could replace previous files." 
else
	mkdir ${outDir}
fi
mv *.fits ${outDir}
mv *.psf ${outDir}
mv *.xml ${outDir}
echo "Done"
