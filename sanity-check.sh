#!/bin/bash

HELP_MSG="usage: ${0} <input-ba> [optional parameters]"
TIMEOUT=60
AUTCROSS_CMD="timeout ${TIMEOUT} autcross -T ${TIMEOUT}"
BINDIR=$(dirname $(readlink -f $0))
KOFOLA="${BINDIR}/build/src/kofola"

# Check the number of command-line arguments
if [ \( "$#" -lt 1 \) ] ; then
	echo ${HELP_MSG}
	exit 1
fi

INPUT=$1
shift

cat ${INPUT} | ${AUTCROSS_CMD} \
	'autfilt --complement %H >%O' \
	"${KOFOLA} $* %H >%O"
