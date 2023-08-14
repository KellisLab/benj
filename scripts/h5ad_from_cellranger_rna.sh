#!/usr/bin/env bash
## usage: < sample_list.txt xargs -P8 -n1 h5ad_from_cellranger_rna.sh
outdir="$1"
shift;
if [[ -n "$1" && "$1" != -* ]]; then
    ### if we want to rename, use the new name
    sample="$1"
    shift;
else
    sample=$(basename "${outdir}")
fi
MIN_N_GENES=${MIN_N_GENES:-3}
### find refdata from the MRO file
if [[ -z "${REFDATA}" ]]; then
    refdata=$(awk -F'= ' '/reference_path/{gsub(/[" ,]/, "", $2); print $2}' "${outdir}/_invocation")
elif [[ -d "${REFDATA}" ]]; then
    refdata="${REFDATA}"
else
    echo "refdata directory must be set, by either having a 10x style _invocation file, or providing the REFDATA env variable."
    exit 1
fi
tss="${refdata}/regions/tss.bed"
gene_info="${refdata}/star/geneInfo.tab"
gtf="${refdata}/genes/genes.gtf.gz"

if [[ -z "${REFDATA}" ]]; then
    refdata=$(awk -F'= ' '/reference_path/{gsub(/[" ,]/, "", $2); print $2}' "${outdir}/_invocation")
else
    refdata="${REFDATA}"
fi
if [[ ! -d "${refdata}" ]]; then
    echo "Refdata ${refdata} is not a directory!"
    exit 1
fi

if [[ -f "${outdir}/outs/filtered_feature_bc_matrix.h5" ]]; then
    mkdir -p H5AD/filtered
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/filtered_feature_bc_matrix.h5 --sample "${sample}" --output H5AD/filtered/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" --gtf "${gtf}" --min-n-genes="${MIN_N_GENES}" $@
fi

if [[ -f "${outdir}/outs/raw_feature_bc_matrix.h5" ]]; then
    mkdir -p H5AD/raw
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/raw_feature_bc_matrix.h5 --sample "${sample}" --output H5AD/raw/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" --gtf "${gtf}" --no-use-velocyto --no-use-scrublet $@
fi

if [[ -f "${outdir}/outs/cellbender_filtered.h5" ]]; then
    mkdir -p H5AD/cellbender
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/cellbender_filtered.h5 --sample "${sample}" --output H5AD/cellbender/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}"  --gtf "${gtf}" --min-n-genes="${MIN_N_GENES}" $@
fi

if [[ -f "${outdir}/outs/cellbender.h5" ]]; then
    mkdir -p H5AD/cb_raw
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/cellbender.h5 --sample "${sample}" --output H5AD/cb_raw/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" --gtf "${gtf}" --min-n-genes="${MIN_N_GENES}" $@
fi
