#!/bin/bash

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
file=
verbose=0

while :; do
    case $1 in
        -h|-\?|--help)
            echo ""
            # echo "    -m, --model: run Torsten example model tests"
            echo "     -u, --unit: run Torsten unit tests"
            # echo "-s, --signature: test Torsten signature"
            exit
            ;;
        # -s|--signature)
        #     echo 'Torsten function signature tests'
        #     pushd cmdstan/stan/
        #     ./runTests.py -j3 src/test/unit/lang/parser/torsten_functions_test.cpp
        #     popd
        #     ;;
        -u|--unit)
            echo ""
            echo 'Torsten library tests'
            pushd cmdstan/stan/lib/stan_math/
            ./runTests.py -j3 stan/math/torsten/test/unit
            popd
            ;;
        # -m|--model)
        #     echo ""
        #     echo 'Torsten model tests'
        #     pushd tests/
        #     rm *.rds
        #     popd
        #     for f in $( ls tests/test_*.R ); do
        #         rcode='testthat::test_file('"\"${f}\""')'
        #         Rscript -e ${rcode}
        #     done
        #     ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            echo ""
            # echo "    -m, --model: run Torsten example model tests"
            echo "     -u, --unit: run Torsten unit tests"
            # echo "-s, --signature: test Torsten signature"
            exit
            ;;
        *)               # Default case
            break
    esac
    shift
done


