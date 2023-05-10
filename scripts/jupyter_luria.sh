#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH -J jupyter
export JUPYTER_LURIA_PREFIX=${JUPYTER_LURIA_PREFIX:-"99"}
if hostname -s | grep -qE '^b[0-9]+$'; then
    ### we are on a compute node
    conda_env=${1:-${CONDA_DEFAULT_ENV}}
    PORT=$(hostname -s | tr -d 'b' | xargs printf "${JUPYTER_LURIA_PREFIX}%02d\n")
    echo "Activating Jupyter notebook on (${conda_env}) for port ${PORT} in dir $PWD"
    source ~/src/conda_init.sh
    conda activate ${conda_env}
    jupyter notebook --no-browser --port=${PORT}
else
    PORT=$(echo $1 | tr -d 'b' | xargs printf "${JUPYTER_LURIA_PREFIX}%02d\n")
    ssh ${ARGS} -N -L ${PORT}:localhost:${PORT} "$@"
fi
