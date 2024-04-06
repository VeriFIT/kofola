#!/bin/bash
to_duration="600s"
exec_dir_prefix="../build/src/"
stats="/mnt/c/Users/ondro/OneDrive/Desktop/Výskum/Projektova_praxe/perf_dir/WSL2-Linux-Kernel/tools/perf/perf"

for FILE in ../inclusion_test_benchmark/*; do
    # Get the last 5 characters of the file's name
    last5="${FILE: -5}"

    # Check if last5 is equal to "A.hoa"
    if [ "$last5" == "A.hoa" ]; then
        fileB="${FILE%?????}B.hoa"
        echo "$FILE $fileB"
        a=$("$stats" stat -r 5 timeout "$to_duration" "$exec_dir_prefix"kofola $FILE $fileB --inclusion 2>&1 | awk '/msec/ {print $1}')
        #b=$("$stats" stat -r 5 timeout "$to_duration" "$exec_dir_prefix"kofola_spot $FILE $fileB --inclusion 2>&1 | awk '/msec/ {print $1}')
        #c=$(timeout "1000s" "$exec_dir_prefix"kofola_spot $FILE $fileB --inclusion 2>&1)
	#$exec_dir_prefix"kofola $FILE --complement
        b=$(timeout "$to_duration" "$exec_dir_prefix"kofola $FILE $fileB --inclusion)
	   # echo $a " vs " $b
        echo $a
	echo $b
        echo "----------------------------"
    fi
done
