#!/usr/bin/env bash

batchdir=""
RMSK="/net/bmc-lab5/data/kellis/group/Benjamin/ref/RepeatMasker_GRCh38.gtf"
ZGTF=""

OPTIONS=b:r:g:
LONGOPTS=batch:,repeatmasker:,gtf:

# Use getopt and store the output into $PARSED
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    exit 2
fi

# Read getopt’s output this way to handle the quoting right
eval set -- "$PARSED"

# Now go through all the options with a case and use shift to analyze one option at a time.
while true; do
    case "$1" in
        -b|--batch)
            batchdir="`readlink -f $2`"
            shift 2
            ;;
	-r|--repeatmasker)
	    RMSK="`readlink -f $2`"
	    shift 2
	    ;;
	-g|--gtf)
	    ZGTF="`readlink -f $2`"
	    shift 2
	    ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: Invalid option: $1" >&2
            exit 3
            ;;
    esac
done

# Check if batch was set
if [[ -z "${batchdir}" ]]; then
    echo "Error: --batch or -b is required." >&2
    exit 1
fi

if ! [[ -f "${batchdir}/aggr.csv" ]]; then
    echo "Error: aggr.csv in batch dir is required." >&2
    exit 1
fi

if ! [[ -f "${RMSK}" ]]; then
    echo "Error: --repeatmasker or -r is required." >&2
    exit 1
fi

nsample=$(wc -l < ${batchdir}/aggr.csv)
echo "#!/usr/bin/env bash"
echo "#SBATCH --array=1-$((nsample-1))%5"
echo "#SBATCH -J velocyto"
echo "cd ${batchdir}"
BENJ_CONDA_DIR=$(dirname $(dirname $CONDA_EXE))
cat <<EOF
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="\$('$CONDA_EXE' 'shell.bash' 'hook' 2> /dev/null)"
if [ \$? -eq 0 ]; then
    eval "\$__conda_setup"
else
    if [ -f "$BENJ_CONDA_DIR/etc/profile.d/conda.sh" ]; then
        . "$BENJ_CONDA_DIR/etc/profile.d/conda.sh"
    else
        export PATH="`dirname $CONDA_EXE`/:\$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
EOF

echo "conda activate ${CONDA_DEFAULT_ENV}"
echo ""
echo "SAMPLEDIR=\$(< aggr.csv awk -F, -v TASK=\$((SLURM_ARRAY_TASK_ID + 1)) 'NR == TASK { print \$1 }')"
echo "OUTDIR=\"\$SAMPLEDIR/outs/\""
echo "echo Sampledir: \"\${SAMPLEDIR}\""
echo "echo Outdir: \"\${OUTDIR}\""
echo "GTF=\$(mktemp --suffix .gtf)"
if [[ -f "$ZGTF" ]]; then
    echo "zcat \"$ZGTF\" > \"\$GTF\""
else
    echo "REFDATA=\$(awk -F'= ' '/reference_path/{gsub(/[\" ,]/, \"\", \$2); print \$2}' \"\${SAMPLEDIR}/_invocation\")"
    echo "zcat \"\${REFDATA}/genes/genes.gtf.gz\" > \"\$GTF\""
fi
echo "if [[ -f \"\${OUTDIR}/gex_possorted_bam.bam\" || -L \"\${OUTDIR}/gex_possorted_bam.bam\" ]]; then"
echo "  \`type -P time\` velocyto run -m \"${RMSK}\" \"\${OUTDIR}/gex_possorted_bam.bam\" \"\${GTF}\" -o \"\${OUTDIR}/velocyto\""
echo "elif [[ -f  \"\${OUTDIR}/gex_possorted.bam\" || -L  \"\${OUTDIR}/gex_possorted.bam\" ]]; then"
echo "  \`type -P time\` velocyto run -m \"${RMSK}\" \"\${OUTDIR}/gex_possorted.bam\" \"\${GTF}\" -o \"\${OUTDIR}/velocyto\""
echo "else"
echo "  \`type -P time\` velocyto run -m \"${RMSK}\" \"\${OUTDIR}/possorted_genome_bam.bam\" \"\${GTF}\" -o \"\${OUTDIR}/velocyto\""
echo "fi"
echo "unlink \"\$GTF\""
