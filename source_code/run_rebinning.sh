#!/bin/bash
# this is intended for running DPS jobs; the input directory is where a single file has been pulled because download=TRUE in the algorithm_config.yaml file
mkdir output

#ATL03 File Path
FILENAME03=$(ls -d input_atl03/*)
FILENAME08=$(ls -d input_atl08/*)


#ATL08 File Path

basedir=$( cd "$(dirname "$0")" ; pwd -P )

echo "basedir: ${basedir}"
echo "PWD: ${PWD}"

INPUTFILE03="${PWD}/${FILENAME03[0]}"
FILELIST03=($INPUTFILE03)
INPUTFILE08="${PWD}/${FILENAME08[0]}"
FILELIST08=($INPUTFILE08)
OUTPUTFILE="${PWD}/output/$3"
echo "FILELIST: ${FILELIST}"
echo "OUTPUT FILE: ${OUTPUTFILE}"

pip install pyproj
pip install 

cd ${basedir}/PHOREAL/source_code/

python ${basedir}/icesatBin.py $FILELIST03 $FILELIST08 $OUTPUTFILE





