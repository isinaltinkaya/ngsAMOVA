#!/usr/bin/env bash
################################################################################
# runTests.sh
#
set -uo pipefail


SCRIPTPATH=$(realpath "$0")
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



initMainLog(){
	local id=${1}
	printf "\n\n"
	printf "###############################################################################\n"
	printf "# RUNNING TEST: ${id}\n"
	printf "\n\n"
}


printMainLog(){
	local msg=${@}
	printf "\n# ${msg}\n"
}

testSuccess(){
	local id=${1}
	printf "${GREEN}\n\n"
	printf "# FINISHED ${id} -> OK\n"
	printf "${NOCOLOR}"
	printf "###############################################################################\n"
	printf "\n\n"
}

testFail(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local logFile=${4}
	local diffFile=${5}

	printf "\n\n"
	printf "${RED}"
	printf "###############################################################################\n"
	printf "# TEST ${id}: FAILED\n"
	printMainLog "Output file:\n${outFile}"
	printMainLog "Reference file:\n${refFile}"
	printMainLog "Log file:\n${logFile}"
	printMainLog "Diff file:\n${diffFile}"
	printf "###############################################################################\n"
	printf "${NOCOLOR}"
	printf "\n\n"
	exit 1
}

runTestDiff(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local outPref=${TESTWD}/${id}
	local logFile=${outPref}.log
	local diffFile=${outPref}.diff

	diff -s ${outFile} ${refFile} > ${diffFile} 2>&1

	if [ $? -eq 0 ]; then
		testSuccess ${id}
	else
		testFail ${id} ${outFile} ${refFile} ${logFile} ${diffFile}
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
	local cmd="${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args}"
	initMainLog ${id}
	printMainLog "Command:\n${cmd}"

	${cmd} > ${logFile} 2>&1

}


##@@
# 	all:
#         $(RM) -rvf testwd;
#         mkdir -pv testwd/logs;
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 1]: EM, then AMOVA with 1 levels in metadata, print distance matrix\n";
#         ../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --printDistanceMatrix 1 --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1 2> testwd/logs/test_s9_d1_1K_mtd1.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 2]: Read the data from distance matrix file created in the previous step, then AMOVA with 1 levels in metadata\n";
#         ../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm testwd/test_s9_d1_1K_mtd1.distance_matrix.csv --minInd 2 -doDist 3 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_indm 2> testwd/logs/test_s9_d1_1K_mtd1_indm.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd1_indm.amova.csv reference/test_s9_d1_1K_mtd1_indm.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 3]: EM, then AMOVA with 2 levels in metadata\n";
#         ../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_2lvl.tsv -f "Individual~Region/Population" -out testwd/test_s9_d1_1K_mtd2 2> testwd/logs/test_s9_d1_1K_mtd2.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd2.amova.csv reference/test_s9_d1_1K_mtd2.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 4]: Read the data from distance matrix file, then do 1-level AMOVA using formula from metadata with 2 levels\n";
#         ../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm test_s9_d1_1K_mtd1_maxIter100_tole10.distance_matrix.csv --minInd 2 -doDist 3 --nThreads 0 -m metadata_with_header_2lvl.tsv -f "Individual~Population" -out testwd/indm_test_s9_d1_1K_mtd1 2> testwd/logs/indm_test_s9_d1_1K_mtd1.log;
#         bash -c "diff testwd/indm_test_s9_d1_1K_mtd1.amova.csv testwd/test_s9_d1_1K_mtd1.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 5]: Use genotypes in VCF file to perform AMOVA\n";
#         ../ngsAMOVA -doEM 0 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --minInd 2 -doDist 2 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt 2> testwd/logs/test_s9_d1_1K_mtd1_gt.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 6]: Same as test 5, but use bcf input with more data and specify regions to use to get the same results as test 5\n";
#         ../ngsAMOVA -doEM 0 -doAMOVA 1 --in-vcf test_s9_d1_1K_extra.bcf --minInd 2 -doDist 2 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt --region "chr22" 2> testwd/logs/test_s9_d1_1K_mtd1_gt.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "[TEST 7]: 40 individuals, multithreaded EM\n";
#         ../ngsAMOVA -doEM 1 -doAMOVA 1 -doDist 1 --maxEmIter 50 --em-tole 1e-5 -P 4 --in-vcf sim_demes_v2-model1-1-rep0-d2.bcf -m sim_demes_v2-model1_metadata_2lvl_with_header.tsv -f "Individual~Region/Population" -o testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5 2> testwd/logs/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.log;
#         bash -c "diff testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv reference/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 8]: Neighbor Joining with individuals as tips. Use vcf individual names as tip labels\n";
#         ../ngsAMOVA -doDist 2 --in-vcf test_s9_d1_1K.vcf -f "Individual~Region/Population" -doPhylo 1 -o testwd/test_s9_d1_1K_mtd2_nj 2> testwd/logs/test_s9_d1_1K_mtd2_nj.log;
#         bash -c "diff testwd/test_s9_d1_1K_mtd2_nj.newick reference/test_s9_d1_1K_mtd2_nj.newick";
# #
#         @printf "\n\n\n";
#         @printf "===============================================================================\n";
#         @printf "\n\n\n[TEST 9]: Same as TEST1 but with same sites divided into 2 contigs\n";
#         ../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K_2contigs.vcf --printDistanceMatrix 1 --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_2contigs_mtd1 2> testwd/logs/test_s9_d1_1K_2contigs_mtd1.log;
#         bash -c "diff testwd/test_s9_d1_1K_2contigs_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
#         @printf "\n\n\n";


##@@



#TODO
# test input gl reading with unobserved notation <*> (-doUnobserved 1)
# INFILENAME="data1.vcf"

################################################################################
# TEST1
ID="test1"
# INFILENAME="test_s9_d1_1K.vcf"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
# local inFile=${SCRIPTDIR}/data/${inFileName}
INOPT="--in-vcf"



doEm=1
doAmova=1
printDistanceMatrix=1
minInd=2
doDist=1
maxEmIter=100
emTole="1e-10"
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
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
-f ${formula}
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
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f ${formula}
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
metadataFile=${DATADIR}/metadata_with_header_2lvl.tsv
formula="Individual~Region/Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
--maxEmIter ${maxEmIter} \
--em-tole ${emTole} \
-m ${metadataFile} \
-f ${formula}
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
metadataFile=${DATADIR}/metadata_with_header_2lvl.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
--nThreads ${nThreads} \
-m ${metadataFile} \
-f ${formula}
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
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f ${formula}
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
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
formula="Individual~Population"

ARGS=" \
--verbose 0 \
-doEM ${doEm} \
-doAMOVA ${doAmova} \
--minInd ${minInd} \
-doDist ${doDist} \
-m ${metadataFile} \
-f ${formula}
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
metadataFile=${DATADIR}/sim_demes_v2-model1_metadata_2lvl_with_header.tsv
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
-f ${formula}
"

runTest ${ID} ${INFILENAME} ${INOPT} "${ARGS}"
runTestDiff ${ID} ${TESTWD}/${ID}.amova.csv ${SCRIPTDIR}/reference/${ID}/${ID}.amova.csv

################################################################################
# TEST 8
ID="test8"
INFILENAME=${DATADIR}/test_s9_d1_1K.vcf
INOPT="--in-vcf"

doDist=2
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
formula="Individual~Region/Population"
doPhylo=1

ARGS=" \
--verbose 0 \
-doDist ${doDist} \
-f ${formula} \
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
metadataFile=${DATADIR}/metadata_with_header_1lvl.tsv
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
-f ${formula}
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




