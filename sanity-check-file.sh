#!/bin/bash

HELP_MSG="usage: ${0} <file-with-auts> [optional parameters]"

# Check the number of command-line arguments
if [ \( "$#" -lt 1 \) ] ; then
	echo ${HELP_MSG}
	exit 1
fi

INPUT=$1
shift
params=$*

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

TMP=$(mktemp)
j=1
cat ${INPUT} | while read line ; do
	echo -n "$j: checking ${line}..."
	./sanity-check.sh ${line} ${params} > ${TMP} 2>&1
	retval=$?
	if [ \( ${retval} == 0 \) ] ; then
		echo -e "${GREEN}OK${NC}"
	elif $(grep -q 'Assertion failed' ${TMP}) ; then
		echo -e "${RED}INTERNAL ERROR${NC}"
	else
		echo -e "${RED}FAIL${NC}"
		cat ${TMP}
	fi
	j=$((j+1))
done
