#!/bin/bash

HELP_MSG="usage: ${0} [--all] <input-ba> [optional parameters]"
TIMEOUT=60
AUTCROSS_CMD="timeout ${TIMEOUT} autcross -T ${TIMEOUT}"

# Check the number of command-line arguments
if [ \( "$#" -lt 1 \) ] ; then
	echo ${HELP_MSG}
	exit 1
fi

ALL=0
if [ \( "$1" == "--all" \) ] ; then
	if [ \( "$#" -lt 2 \) ] ; then
		echo ${HELP_MSG}
		exit 1
	fi

	ALL=1
	shift
fi

INPUT=$1
shift
if [ \( ${ALL} == 1 \) ] ; then
	cat ${INPUT} | ${AUTCROSS_CMD} \
	'./bin/autfilt --complement %H >%O' \
	"./bin/ranker $* %H >%O" \
	'./bin/seminator --complement %H >%O' \
	'./bin/ltl2dstar --complement-input=yes --input=nba --output=nba -H %H %O' \
	'./bin/roll-autcross-wrap.sh %H %O' \
	'./bin/goal-autcross-wrap.sh %H %O -m safra' \
	'./bin/goal-autcross-wrap.sh %H %O -m piterman' \
	'./bin/goal-autcross-wrap.sh %H %O -m rank -tr -ro' \
	'./bin/goal-autcross-wrap.sh %H %O -m fribourg' \
	'./bin/ranker-tight %H >%O' \
	'./bin/ranker-composition %H >%O'
else
	echo  running "./kofola --algo=comp $* %H >%O"
	cat ${INPUT} | ${AUTCROSS_CMD} \
		'autfilt --complement %H >%O' \
		"./kofola --algo=comp $* %H >%O"
		# './kofola --algo=comp %H >%O'
fi
