#!/bin/sh

# module purge
# module load snakemake/5.4.5-python-3.6.8


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
${SNAKE_CMD}
