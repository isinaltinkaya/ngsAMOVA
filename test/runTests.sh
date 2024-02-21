#!/usr/bin/env bash
################################################################################
# runTests.sh
#
set -uo pipefail

SCRIPTPATH=$(realpath "$0")

TESTTYPE=${1:-"regular"}


SCRIPTDIR=$(dirname "$SCRIPTPATH")
DATADIR=$(realpath "$SCRIPTDIR/data")
TESTWD=$(realpath "$SCRIPTDIR/testwd")
EXEC=$(realpath "$SCRIPTDIR/../ngsAMOVA")

GREEN='\033[0;32m'
RED='\033[0;31m'
NOCOLOR='\033[0m'

rm -rfv ${TESTWD}/
mkdir -pv ${TESTWD}/
echo ${TESTWD}

testExec(){
	if ! command -v ${EXEC} &> /dev/null; then
		printf "${RED}\n\n"
		printf "###############################################################################\n"
		printf "# ERROR\n"
		printf "# Executable could not be found at:\n"
		printf "${EXEC}\n"
		printf "###############################################################################\n"
		printf "\n\n"
		printf ${NOCOLOR}
		exit 1;
	fi
}

testValgrind(){
	if ! command -v valgrind &> /dev/null; then
		echo "valgrind could not be found";
		exit 1;
	fi
}


if [ ${TESTTYPE} == "regular" ]; then
	testExec
elif [ ${TESTTYPE} == "vg" ]; then
	testExec
	testValgrind
elif [ ${TESTTYPE} == "all" ]; then
	testExec
	testValgrind
else
	echo "Unknown test type: ${TESTTYPE}"
	exit 1;
fi

printf "\n\n"
printf "###############################################################################\n"
printf "# Starting ${TESTTYPE} tests\n"
printf "\n# Script path:\n"
printf "${SCRIPTPATH}\n"
printf "\n# Script directory:\n"
printf "${SCRIPTDIR}\n"
printf "\n# Data directory:\n"
printf "${DATADIR}\n"
printf "\n# Test working directory:\n"
printf "${TESTWD}\n"
printf "\n# Executable:\n"
printf "${EXEC}\n"
printf "\n# Test data directory:\n"
printf "${DATADIR}\n"
printf "###############################################################################\n"
printf "\n\n"




runTestDiff(){
	if [ ${TESTTYPE} == "regular" ]  || [ ${TESTTYPE} == "all" ]; then
		local id=${1}
		local outFile=${2}
		local refFile=${3}
		local outPref=${TESTWD}/${id}
		local logFile=${outPref}.log
		local diffFile=${outPref}.diff

		local diffcmd="diff -s ${outFile} ${refFile} > ${diffFile} 2>&1"

		printf "# ${id} -> Running diff\n"
		printf "# Command:\n${diffcmd}\n"

		eval ${diffcmd}

		if [ $? -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Diff: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "\n\n"
			printf "${RED}"
			printf "###############################################################################\n"
			printf "# ${id} FAILED\n"

			printf "\n# Command:\n${diffcmd}\n"
			printf "\n# Output file:\n${outFile}\n"
			printf "\n# Reference file:\n${refFile}\n"
			printf "\n# Log file:\n${logFile}\n"
			printf "\n# Diff file:\n${diffFile}\n"
			printf "###############################################################################\n"
			printf "${NOCOLOR}"
			printf "\n\n"
			exit 1
		fi
	fi
}


runTest(){

	local id=${1}
	local inFileName=${2}
	local inOpt=${3}
	local args=${4}
	local outPref=${TESTWD}/${id}
	local outFile=${outPref}.vcf
	local refFile=${SCRIPTDIR}/reference/${id}/${id}.vcf
	local diffFile=${outPref}.diff
	local logFile=${outPref}.log

	local cmd;

	if [ ${TESTTYPE} == "regular" ]; then
		cmd="${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	elif [ ${TESTTYPE} == "vg" ]; then
		cmd="valgrind --leak-check=full -q --error-exitcode=1 --log-fd=9 9>> ${outPref}.vg.log ${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	elif [ ${TESTTYPE} == "all" ]; then
		cmd="valgrind --leak-check=full -q --error-exitcode=1 --log-fd=9 9>> ${outPref}.vg.log ${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	fi

	printf "\n\n"
	printf "###############################################################################\n"
	printf "# ${id} -> Running ${TESTTYPE} test \n"
	printf "\n\n"

	printf "# Command:\n${cmd}\n"

	
	if [ ${TESTTYPE} == "vg" ] || [ ${TESTTYPE} == "all" ]; then

		local vgFile=${outPref}.vg.log
		local vgrun;
		eval ${cmd}
		exitCode=$?

		if [ ${exitCode} -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Valgrind test: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "${RED}\n\n"
			printf "# ${id} -> Valgrind test: FAILED\n"
			printf "# Valgrind output:\n"
			printf "${vgFile}\n"

			printf "\n\n"
			printf "${NOCOLOR}"
			exit 1

		fi

	elif [ ${TESTTYPE} == "regular" ]; then
		eval ${cmd}
		exitCode=$?

		if [ ${exitCode} -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Run: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "${RED}\n\n"
			printf "# ${id} -> Run: FAILED\n"
			printf "# Log file:\n"
			printf "${logFile}\n"

			printf "\n\n"
			printf "${NOCOLOR}"
			exit 1

		fi

	fi

}

################################################################################
# TEST1
ID="test1"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"


doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
metadataFile=${DATADIR}/metadata_Individual_Population.tsv
formula="Individual ~Population "


ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
-m ${metadataFile} \
-f '${formula}'
"


runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
runTestDiff ${ID} ${TESTWD}/${ID}.distance_matrix.csv ${SCRIPTDIR}/reference/${ID}/${ID}.distance_matrix.csv



################################################################################
# TEST 2
ID="test2"
INFILENAME=${TESTWD}/test1.distance_matrix.csv
INOPT="--in-dm"

doEm=0
doAmova=1
minInd=2
doDist=3
metadataFile=${DATADIR}/metadata_Individual_Region_Population_Subpopulation.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 3
ID="test3"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
metadataFile=${DATADIR}/metadata_Individual_Region_Population.tsv
formula="Individual~Region /Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 4
ID="test4"
INFILENAME=${DATADIR}/test_s9_d1_1K_mtd1_maxIter100_tole10.distance_matrix.csv

INOPT="--in-dm"

doEm=0
doAmova=1
minInd=2
doDist=3
nThreads=0
metadataFile=${DATADIR}/metadata_Individual_Region_Population_Subpopulation_groupNotUniq.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
--nThreads ${nThreads} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${TESTWD}/test1.amova.csv

################################################################################
# TEST 5
ID="test5"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"

doEm=0
doAmova=1
minInd=2
doDist=2
metadataFile=${DATADIR}/metadata_Individual_Population.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 6
ID="test6"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"

doEm=0
doAmova=1
minInd=2
doDist=2
metadataFile=${DATADIR}/metadata_Individual_Population.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 7
ID="test7"
INFILENAME=${DATADIR}/sim_demes_v2-model1-1-rep0-d2.bcf
INOPT="--in-vcf"

doEm=1
doAmova=1
doDist=1
maxEmIter=50
emTole="1e-5"
nThreads=4
metadataFile=${DATADIR}/sim_demes_v2-model1_metadata_Individual_Region_Population.tsv 
formula="Individual~Region/Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
-P ${nThreads} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 8
ID="test8"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"

doDist=2
metadataFile=${DATADIR}/metadata_Individual_Population.tsv
formula="Individual~Region/Population"
doPhylo=1

ARGS=" \
--verbose 0 \
-doDist ${doDist} \
-f '${formula}' \
-doPhylo ${doPhylo}
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.newick ${SCRIPTDIR}/reference/${ID}/${ID}.newick

################################################################################
# TEST 9
ID="test9"
INFILENAME=${DATADIR}/test_s9_d1_1K_2contigs.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
metadataFile=${DATADIR}/metadata_Individual_Population.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
# ###############################################################################

################################################################################
# TEST 10
ID="test10"
INFILENAME=${DATADIR}/data0.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
runTestDiff ${ID} ${TESTWD}/${ID}.distance_matrix.csv ${SCRIPTDIR}/reference/${ID}/${ID}.distance_matrix.csv


# ###############################################################################

################################################################################
# TEST 11
ID="test11"
INFILENAME=${DATADIR}/data1.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

# ###############################################################################

################################################################################
# TEST 12
ID="test12"
INFILENAME=${DATADIR}/data2.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/test15/test15.amova.csv

# ###############################################################################


################################################################################
# TEST 13
ID="test13"
INFILENAME=${DATADIR}/data3.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
# ###############################################################################


################################################################################
# TEST 14
ID="test14"
INFILENAME=${DATADIR}/data4.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
# ###############################################################################


################################################################################
# TEST 15
ID="test15"
INFILENAME=${DATADIR}/data5.vcf
INOPT="--in-vcf"

doEm=1
doAmova=1
printDistanceMatrix=1
printJointGenotypeCountMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--printJointGenotypeCountMatrix ${printJointGenotypeCountMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
# ###############################################################################


################################################################################
# TEST 16
ID="test16"
INFILENAME=${DATADIR}/data0.truth.vcf
INOPT="--in-vcf"

doAmova=1
printDistanceMatrix=1
printJointGenotypeCountMatrix=1
minInd=2
doDist=2
rmInvarSites=1
metadataFile=${DATADIR}/data0to5_metadata.txt
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--printDistanceMatrix ${printDistanceMatrix} \
--printJointGenotypeCountMatrix ${printJointGenotypeCountMatrix} \
--minInd ${minInd} \
-doDist ${doDist} \
--rm-invar-sites ${rmInvarSites} \
-m ${metadataFile} \
-f '${formula}'
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv
# ###############################################################################



${EXEC}


echo "__"
# check if logfile looks ok 
wc -l ${TESTWD}/test1.log
echo "__"

wait

printf "\n\n"
echo -e "                                                                           
 -------------------
< ALL TESTS PASSED! >
 -------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/
                ||----w |
                ||     ||
"

printf "\n\n"




