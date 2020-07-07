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

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
    JOB_ID = ${SLURM_JOB_ID}
else
    JOBID=$(date +%Y-%m-%d.%H-%M-%S)  
fi
LOGFILE="logs/snakelog.${JOBID}.out"

# Run workflow
echo ""
echo "=========================================="
echo "Running Snakemake:"
echo "${SNAKE_CMD}"
echo ""
echo "Logging output to: ${LOGFILE}"
echo "=========================================="
echo ""

module --quiet purge
module --quiet load snakemake && ${SNAKE_CMD} &> ${LOGFILE}
