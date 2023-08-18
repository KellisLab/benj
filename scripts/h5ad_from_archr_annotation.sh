#!/usr/bin/env bash
# usage: export PEAKS=/path/to/peaks; export BLACKLIST=/path/to/blacklist; < sample_list.txt xargs -P8 -I{} bash -c 'h5ad_from_archr_annotation.sh {}'

outdir="$1"
shift;

if [[ -n "$1" && "$1" != -* ]]; then
    ### if we want to rename, use the new name
    sample="$1"
    shift;
else
    sample=$(basename "${outdir}")
fi

# if [[ -z "${REFDATA}" ]]; then
#     refdata=$(awk -F'= ' '/reference_path/{gsub(/[" ,]/, "", $2); print $2}' "${outdir}/_invocation")
# elif [[ -d "${REFDATA}" ]]; then
#     refdata="${REFDATA}"
# else
#     echo "refdata directory must be set, by either having a 10x style _invocation file, or providing the REFDATA env variable."
#     exit 1
# fi

md=$(dirname "${outdir}")
md="${md}/ArrowFiles/archr_metadata.tsv.gz"
if [[ ! -f "${md}" ]]; then
    echo "ArchR metadata at ${md} does not exist!"
    exit 1
fi


if [[ -z "${PEAKS}" ]]; then
    echo "PEAKS must be set and be a file!"
    exit 1
elif [[ ! -f "${PEAKS}" ]]; then
    echo "PEAKS must be set and be a file!"
    exit 1
fi

extraopts=""
if [[ -f "${outdir}/outs/atac_fragments.tsv.gz" ]]; then
    ## multi-ome
    fragments="${outdir}/outs/atac_fragments.tsv.gz"
elif [[ -f "${outdir}/outs/fragments.tsv.gz" ]]; then
    ## single-ome
    fragments="${outdir}/outs/fragments.tsv.gz"
else
    echo "Fragments not found at ${outdir}/outs/*fragments.tsv.gz"
    exit 1
fi

mkdir -p "H5AD/ATAC"
output="H5AD/ATAC/${sample}.h5ad"
if [[ -f "${BLACKLIST}" ]]; then
    h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${PEAKS}" -o "${output}" -b "${BLACKLIST}"
else
    h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${PEAKS}" -o "${output}"
fi
