#!/usr/bin/env bash
## usage: < sample_list.txt xargs -P8 -n1 h5ad_from_cellranger_rna.sh
mkdir -p H5AD/{filtered,raw,cellbender}
outdir="$1"
shift;
if [[ -n "$1" && "$1" != -* ]]; then
    sample="$1"
    shift;
else
    sample=$(basename "${outdir}")
fi

### find refdata from the MRO file
refdata=$(awk -F'= ' '/reference_path/{gsub(/[" ,]/, "", $2); print $2}' "${outdir}/_invocation")
tss="${refdata}/regions/tss.bed"
gene_info="${refdata}/star/geneInfo.tab" 

if [[ -f "${outdir}/outs/filtered_feature_bc_matrix.h5" ]]; then
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/filtered_feature_bc_matrix.h5 --sample "${sample}" --output H5AD/filtered/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" $@
fi

if [[ -f "${outdir}/outs/raw_feature_bc_matrix.h5" ]]; then
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/raw_feature_bc_matrix.h5 --sample "${sample}" --output H5AD/raw/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" $@
fi

if [[ -f "${outdir}/outs/cellbender_filtered.h5" ]]; then
    h5ad_from_cellranger_rna.py --h5 "${outdir}"/outs/cellbender_filtered.h5 --sample "${sample}" --output H5AD/cellbender/"${sample}.h5ad" --tss "${tss}" --gene-info "${gene_info}" $@
fi
