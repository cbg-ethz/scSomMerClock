#!/bin/sh

SNAKE_CMD="snakemake --keep-going --restart-times=0"
while [ "$1" != "" ]; do
    case $1 in
        -c | --config ) shift        
                        SNAKE_CMD+=" --configfile ${1}"
                        ;;                      
        *)              SNAKE_CMD+=" ${1}"
    esac
    shift
done

# Run workflow
echo ''
echo 'Running Snakemake:'
echo "${SNAKE_CMD}"
module purge
module load snakemake && ${SNAKE_CMD}
