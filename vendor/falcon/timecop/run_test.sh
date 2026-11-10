#!/bin/bash

export ROOT="$(pwd)"

usage() {
    echo "Various possible options :"
    echo -e "\t-p: only perform the parsing of ctgrind results"
    echo -e "\t-n N_CALLERS: set the number of caller to consider with ctgrind (default=20)"
    echo -e "\t-s FILENAME: add a suppression file to ctgrind (default=None)"
    echo -e "\t-d DEPTH: call depth to display when parsing the error tree (default=6)"
    echo -e "\t-b FILENAME: name of the binary to test (default=a.out)"
    exit -1
}

bin=build/test_timecop
VALGRIND_OPTS="--leak-check=full --error-limit=no"

depth=100
parse_only=0
num_callers=20
while getopts "pd:b:s:n:" opt; do
    case "$opt" in
    p)
		parse_only=1
        ;;
    n)
        num_callers=$OPTARGS
        ;;
    s)
        VALGRIND_OPTS="$VALGRIND_OPTS -v --suppressions=$OPTARG"
        ;;
    d)
        depth=$OPTARG
        ;;
    b)
        bin=$OPTARG
        ;;
    *)
        usage
    esac
done
VALGRIND_OPTS="$VALGRIND_OPTS --num-callers=$num_callers"

LOG_PATH="${ROOT}/logs/"
mkdir -p "${LOG_PATH}"

if [[ $parse_only == 0 || -z $bin ]]
then
    valgrind ${VALGRIND_OPTS} ./$bin 2>&1 | grep "==" | sed 's/==.*== //g' > $LOG_PATH/.timecop_raw.log
fi

python "$ROOT"/timecop/parse.py -d $depth -o $LOG_PATH/timecop.log $LOG_PATH/.timecop_raw.log
