#!/usr/bin/env bash
# usage: < sample_list.txt xargs -P8 -I{} bash -c 'h5ad_from_archr_annotation.sh {}'

outdir="$1"
shift;

if [[ -n "$1" && "$1" != -* ]]; then
    ### if we want to rename, use the new name
    sample="$1"
    shift;
else
    sample=$(basename "${outdir}")
fi

if [[ -z "${REFDATA}" ]]; then
    refdata=$(awk -F'= ' '/reference_path/{gsub(/[" ,]/, "", $2); print $2}' "${outdir}/_invocation")
elif [[ -d "${REFDATA}" ]]; then
    refdata="${REFDATA}"
else
    echo "refdata directory must be set, by either having a 10x style _invocation file, or providing the REFDATA env variable."
    exit 1
fi
gtf="${refdata}/genes/genes.gtf.gz"
tss="${refdata}/regions/tss.bed"
md=$(dirname "${outdir}")
atac_dir="${md}/ArrowFiles/"
md="${atac_dir}/archr_metadata.tsv.gz"

if [[ ! -f "${md}" ]]; then
    echo "ArchR metadata at ${md} does not exist!"
    exit 1
fi

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


for tilebed in "${atac_dir}/"tile_*.tsv.gz; do
    # TODO; suffix
    suffix=$(basename "${tilebed}" | sed 's/^tile_//g;s/.tsv.gz$//g')
    hdir="H5AD/TileMatrix${suffix}"
    mkdir -p "${hdir}"
    h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${tilebed}" -o "${hdir}/${sample}.h5ad"
done

## now get the largest tile matrix (for gene score, more accurate)
gs_in=$(ls -1S "${atac_dir}/tile_*.tsv.gz" | head -n1)
gs_in=$(basename "${gs_in}" | sed 's/^tile_/TileMatrix/g;s/.tsv.gz$//g')

Tile4GeneH5AD="H5AD/${gs_in}/${sample}.h5ad"
if [[ -f "${Tile4GeneH5AD}" ]]; then
    mkdir -p "H5AD/GeneScoreMatrix"
    estimate_gene_accessibility.py -i "${Tile4GeneH5AD}" -o "H5AD/GeneScoreMatrix/${sample}.h5ad" --gtf "${gtf}" --tss "${tss}"
fi
