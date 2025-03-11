#!usr/bin/env bash
################################################################################
# runTests.sh
#
set -uo pipefail

################################################################################
# TODO add test for --in-majorminor-fn with --a2f
################################################################################

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NOCOLOR='\033[0m'
nl=$'\n' # ANSI C quoting

################################################################################
# FUNCTIONS

aline(){
	printf "%.0s""${@}" {1..80}
	printf "${nl}"
}

print_cmd(){
	printf "${nl}"
	printf "${YELLOW}"
	printf "# Command:${nl}"
	printf "${@}"
	printf "${NOCOLOR}"
}

# END-OF-FUNCTIONS
################################################################################

SCRIPTPATH=$(realpath "$0")
SCRIPTDIR=$(dirname "$SCRIPTPATH")


TEST_TO_RUN=""

# -e EXEC -t TESTTYPE -d DATADIR -r REFDIR -w TESTWD  -x TEST_TO_RUN
# required: -e -t -d -r -w
# optional: -x
while getopts e:t:d:r:w:x: opt; do
	case ${opt} in
		e)
			EXEC=$(realpath ${OPTARG})
			;;
		t)
			TESTTYPE=${OPTARG}
			;;
		d)
			DATADIR=$(realpath ${OPTARG})
			;;
		r)
			REFDIR=$(realpath ${OPTARG})
			;;
		w)
			TESTWD=$(realpath ${OPTARG})
			;;
		x)
			TEST_TO_RUN=${OPTARG}
			;;
		\?)
			echo "Invalid option: $OPTARG" 1>&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." 1>&2
			exit 1
			;;
	esac
done


exit_with_error(){
	printf "${RED}"
	printf "${nl}${nl}"
	aline "#"
	printf "# ERROR${nl}"
	printf "# ""${1}""${nl}"
	aline "#"
	printf "${nl}${nl}"
	printf "${NOCOLOR}"

	exit 1
}

if  [ -z ${REFDIR+x} ]; then
	# REFDIR=${SCRIPTDIR}/reference
	exit_with_error "Reference directory not set (-r REFDIR)"
fi

if [ -z ${DATADIR+x} ]; then
	# DATADIR=$(realpath "$SCRIPTDIR/data")
	exit_with_error "Data directory not set (-d DATADIR)"
fi

if [ -z ${TESTWD+x} ]; then
	# TESTWD=$(realpath "$SCRIPTDIR/testwd")
	exit_with_error "Test working directory not set (-w TESTWD)"
fi



rm -rfv ${TESTWD}/
mkdir -pv ${TESTWD}/
echo ${TESTWD}


check_exit_code(){
	local id=${1}
	local exitCode=${2}
	local logFile;
	local testDescription;

	if [ ${TESTTYPE} == "vg" ];then
		testDescription="Valgrind test"
		logFile=${TESTWD}/${id}.vg.log

	elif [ ${TESTTYPE} == "regular" ]; then
		testDescription="Execution test"
		logFile=${TESTWD}/${id}.log
	fi

	if [ ${exitCode} -eq 0 ]; then
		printf "${GREEN}"
		printf "${nl}"
		printf "# -> Run for ${id}: OK${nl}"
		printf "${NOCOLOR}"
		printf "${nl}${nl}"
	else
		exit_with_error "${nl}-> Run for ${id}: FAILED${nl}# Log file:${nl}${logFile}"
	fi

}

diff_test_out_ref(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local outPref=${TESTWD}/${id}
	local logFile=${outPref}.log
	local diffFile=${outPref}.diff

	local diffcmd="diff -s ${outFile} ${refFile}
	> ${diffFile} 2>&1"

	printf "# -> Running diff for ${id}${nl}"
	printf "${YELLOW}"
	printf "# Command:${nl}${diffcmd}${nl}"
	printf "${NOCOLOR}"

	eval ${diffcmd}

	if [ $? -eq 0 ]; then
		printf "${GREEN}"
		printf "# ${id} -> Diff: OK${nl}"
		printf "${NOCOLOR}"
		printf "${nl}${nl}"
	else
		printf "${nl}${nl}"
		printf "${RED}"
		aline "#"
		printf "# ${id} FAILED${nl}"
		printf "${nl}# Command:${nl}${diffcmd}${nl}"
		printf "${nl}# Output file:${nl}${outFile}${nl}"
		printf "${nl}# Reference file:${nl}${refFile}${nl}"
		printf "${nl}# Log file:${nl}${logFile}${nl}"
		printf "${nl}# Diff file:${nl}${diffFile}${nl}"
		aline "#"
		printf "${NOCOLOR}"
		printf "${nl}${nl}"
		exit 1
	fi
}


################################################################################

if [ ${TESTTYPE} == "regular" ]; then
	:
elif [ ${TESTTYPE} == "vg" ]; then
	if ! command -v valgrind &> /dev/null; then
		echo "valgrind could not be found";
		exit 1;
	fi
else
	echo "Unknown test type: ${TESTTYPE}"
	exit 1;
fi


if ! command -v ${EXEC} &> /dev/null; then
	printf "${RED}${nl}${nl}"
	aline "#"
	printf "# ERROR${nl}"
	printf "# Executable could not be found at:${nl}"
	printf "${EXEC}${nl}"
	aline "#"
	printf "${nl}${nl}"
	printf ${NOCOLOR}
	exit 1;
fi


printf "${nl}${nl}"
aline "#"
printf "# Starting ${TESTTYPE} tests${nl}"
printf "${nl}# Script path:${nl}"
printf "${SCRIPTPATH}${nl}"
printf "${nl}# Script directory:${nl}"
printf "${SCRIPTDIR}${nl}"
printf "${nl}# Data directory:${nl}"
printf "${DATADIR}${nl}"
printf "${nl}# Test working directory:${nl}"
printf "${TESTWD}${nl}"
printf "${nl}# Executable:${nl}"
printf "${EXEC}${nl}"
printf "${nl}# Test data directory:${nl}"
printf "${DATADIR}${nl}"
aline "#"
printf "${nl}${nl}"







################################################################################


declare -a tests
tests=("test0" "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9" "test10" "test11" "test12" "test13" "test14" "test15" "test16" "test17" "test18" "test19" "test20" "test21" "test22" "test23" "test24" "test25" "test26" "test27")


get_cmd(){
	local id=${1}
	local cmd="${@:2}"
	if [ ${TESTTYPE} == "vg" ];then
		cmd="    valgrind --leak-check=full -q --error-exitcode=1 --log-fd=9 9>> ${TESTWD}/${id}.vg.log ${nl}${cmd}"
	fi
	printf "${cmd}"

}

test0(){
	local id=${FUNCNAME}
	local cmd="${EXEC} -doUnitTests 1"
	cmd=$(get_cmd ${id} ${cmd})
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
}



################################################################################

test1(){
	local id=${FUNCNAME}
	local inFile="${DATADIR}/test_s9_d1_1K_2contigs.vcf" #2contig version of test_s9_d1_1K.vcf
	local cmd="${EXEC}"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Population.tsv"
	cmd="${cmd} ""-f 'Individual ~Population '"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
	diff_test_out_ref ${id} ${TESTWD}/${id}.distance_matrix.txt ${REFDIR}/${id}.distance_matrix.txt
}


################################################################################

test2(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test1.distance_matrix.txt"
	cmd="${cmd} ""--in-dm ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Region_Population_Subpopulation.tsv"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""--amova-euclid 1"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test3(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 2"
	cmd="${cmd} ""--in-majorminor ${DATADIR}/test_s9_d1_1K_majorminor.tsv"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Region_Population.tsv"
	cmd="${cmd} ""-f 'Individual~Region /Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test4(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test1.distance_matrix.txt"
	cmd="${cmd} ""--in-dm ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 2"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Region_Population_Subpopulation_groupNotUniq.tsv"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/test2.amova.csv
}

################################################################################

test5(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 0"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Population.tsv"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test6(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 2"
	cmd="${cmd} ""--in-majorminor ${DATADIR}/test_s9_d1_1K_majorminor.tsv"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 0"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Population.tsv"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test7(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/sim_demes_v2-model1-1-rep0-d2_filled.bcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 5"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""--em-tole 1e-5"
	cmd="${cmd} ""-P 2"
	cmd="${cmd} ""-m ${DATADIR}/sim_demes_v2-model1_metadata_Individual_Region_Population.tsv"
	cmd="${cmd} ""-f 'Individual~Region/Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test8(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""--print-jgtm 1"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-f 'Individual~Region/Population'"
	cmd="${cmd} ""-doPhylo 1"
	cmd="${cmd} ""--print-tree 1"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.newick ${REFDIR}/${id}.newick
}

################################################################################

test9(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K_2contigs.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--amova-euclid 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Population.tsv"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test10(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data0.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
	diff_test_out_ref ${id} ${TESTWD}/${id}.distance_matrix.txt ${REFDIR}/${id}.distance_matrix.txt
}

################################################################################

test11(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data1.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test12(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data2.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--amova-euclid 1"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test13(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data3.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test14(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data4.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test15(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data5.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--print-jgtm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test16(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data0_filled.truth.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--print-jgtm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 0"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f Individual~Population"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}

################################################################################

test17(){
# same data as data4 but with different alleles read with alleles file
# so result is expected to be the same as test14
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data6.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 2"
	cmd="${cmd} ""--in-majorminor ${DATADIR}/data6_majorminor.tsv"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--amova-euclid 1"
	cmd="${cmd} ""-m ${DATADIR}/data0to5_metadata.txt"
	cmd="${cmd} ""-f 'Individual~Population'"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
}


################################################################################

test18(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/data8_1.bcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doBlockBootstrap 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--seed 2"
	cmd="${cmd} ""--block-size 50000"
	cmd="${cmd} ""--print-blocks 1"
	cmd="${cmd} ""--nbootstraps 1"
	cmd="${cmd} ""--allow-mispairs 1"
	cmd="${cmd} ""--min-npairs 0"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.distance_matrix.txt ${REFDIR}/${id}.distance_matrix.txt
}

################################################################################

test19(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test1_ibd_100_withAC.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doIbd 2"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--min-a2c 2"
	cmd="${cmd} ""--ibdlod 3"
	cmd="${cmd} ""--ibdtrim 0"
	cmd="${cmd} ""--errormax 0.001"
	cmd="${cmd} ""--errorprop 0.25"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/${id}.ibd_segments.tsv
}



################################################################################

test20(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test2_ibd_1000_err0_filled.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 2"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doIbd 2"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--min-a2c 2"
	cmd="${cmd} ""--ibdlod 3"
	cmd="${cmd} ""--ibdtrim 0"
	cmd="${cmd} ""--errormax 0.001"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/${id}.ibd_segments.tsv
}


################################################################################

test21(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test2_ibd_1000_err0_filled.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doIbd 1"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--errormax 0"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/${id}.ibd_segments.tsv
}

################################################################################

test22(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test2_ibd_1000_err0_filled.vcf"
	local a2fFile="${DATADIR}/test2_ibd_1000_err0_filled.maf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doIbd 1"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--errormax 0"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""--a2f ${a2fFile}"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/test21.ibd_segments.tsv
}


################################################################################

test23(){

	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test2_ibd_1000_err0_filled_alleles_flipped.vcf"
	local a2fFile="${DATADIR}/test2_ibd_1000_err0.maf"
	local majorMinorFile="${DATADIR}/test2_ibd_1000_err0_majorminor.tsv"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 2"
	cmd="${cmd} ""--in-majorminor ${majorMinorFile}"
	cmd="${cmd} ""-doIbd 1"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--errormax 0"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""--a2f ${a2fFile}"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/test22.ibd_segments.tsv
}


################################################################################

test24(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test3_small_multiallelic_filled_trimmed.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 3"
	cmd="${cmd} ""-doIbd 1"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--errormax 0"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--rm-multiallelic-sites 1"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/${id}.ibd_segments.tsv
}

################################################################################
# should give the same result as test24

test25(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test3_small_biallelic_filled_trimmed.vcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 3"
	cmd="${cmd} ""-doIbd 1"
	cmd="${cmd} ""--print-ibd 1"
	cmd="${cmd} ""--errormax 0"
	cmd="${cmd} ""--errorprop 0"
	cmd="${cmd} ""--rm-invar-sites 1"
	cmd="${cmd} ""--rm-multiallelic-sites 0"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.ibd_segments.tsv  ${REFDIR}/${id}.ibd_segments.tsv
}


################################################################################

test26(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K_2contigs_emptyBlock.bcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""--print-dm 3"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""-doBlockBootstrap 1"
	cmd="${cmd} ""--print-blocks 1"
	cmd="${cmd} ""--verbose 3"
	cmd="${cmd} ""--minInd 9"
	cmd="${cmd} ""--block-size 100000"
	cmd="${cmd} ""--print-bootstrap 1"
	cmd="${cmd} ""--nBootstraps 2"
	cmd="${cmd} ""--seed 3"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.bootstrap_samples.tsv  ${REFDIR}/${id}.bootstrap_samples.tsv
	diff_test_out_ref ${id} ${TESTWD}/${id}.blocks.tsv  ${REFDIR}/${id}.blocks.tsv
	diff_test_out_ref ${id} ${TESTWD}/${id}.distance_matrix.txt  ${REFDIR}/${id}.distance_matrix.txt
}

################################################################################

test27(){
	local id=${FUNCNAME}
	local cmd="${EXEC}"
	local inFile="${DATADIR}/test_s9_d1_1K_2contigs_emptyBlock.bcf"
	cmd="${cmd} ""--in-vcf ${inFile}"
	cmd="${cmd} ""-o ${TESTWD}/${id}"
	cmd="${cmd} ""-doMajorMinor 1"
	cmd="${cmd} ""-doDist 1"
	cmd="${cmd} ""-doJGTM 1"
	cmd="${cmd} ""-doEM 1"
	cmd="${cmd} ""-doAMOVA 1"
	cmd="${cmd} ""-doBlockBootstrap 1"
	cmd="${cmd} ""--print-blocks 1"
	cmd="${cmd} ""--print-bootstrap 1"
	cmd="${cmd} ""--nBootstraps 2"
	cmd="${cmd} ""--seed 3"
	cmd="${cmd} ""--block-size 100000"
	cmd="${cmd} ""--print-amova 1"
	cmd="${cmd} ""--print-dm 1"
	cmd="${cmd} ""--minInd 2"
	cmd="${cmd} ""--maxEmIter 10"
	cmd="${cmd} ""--bcf-src 1"
	cmd="${cmd} ""--amova-euclid 1"
	cmd="${cmd} ""-m ${DATADIR}/metadata_Individual_Population.tsv"
	cmd="${cmd} ""-f 'Individual ~Population '"
	cmd="${cmd} ""${nl}    2>${TESTWD}/${id}.log"
	cmd=$(get_cmd ${id} "${cmd}")
	print_cmd "${cmd}"
	eval ${cmd}
	local exitCode=$?
	check_exit_code ${id} ${exitCode}
	diff_test_out_ref ${id} ${TESTWD}/${id}.amova.csv ${REFDIR}/${id}.amova.csv
	diff_test_out_ref ${id} ${TESTWD}/${id}.distance_matrix.txt ${REFDIR}/${id}.distance_matrix.txt
}

################################################################################

# if TEST_TO_RUN is defined, run only the test with the corresponding index
if [ -n "${TEST_TO_RUN}" ]; then
	if [ ${TEST_TO_RUN} -ge 0 ] && [ ${TEST_TO_RUN} -lt ${#tests[@]} ]; then
		aline "_"
		# printf "${nl}"
		# printf "# [${tests[${TEST_TO_RUN]}]${nl}${nl}"
		# printf "# -> Running test for ${tests[${TEST_TO_RUN}]}"
		${tests[${TEST_TO_RUN}]}
	else
		exit_with_error "Invalid test index: ${TEST_TO_RUN}"
	fi
	exit 0
fi
# if []

for testNo in "${!tests[@]}"; do
	aline "_"
	printf "${nl}"
	printf "# [${tests[${testNo}]}]${nl}${nl}"
	printf "# -> Running test for ${tests[${testNo}]}"
	${tests[${testNo}]}
done


################################################################################
aline "_"

${EXEC}


################################################################################
aline "_"

# check if logfile looks ok 
wc -l ${TESTWD}/test1.log

################################################################################
aline "_"


printf "${nl}${nl}"
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

printf "${nl}${nl}"




