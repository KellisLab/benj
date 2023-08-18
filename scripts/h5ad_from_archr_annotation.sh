#!/usr/bin/env bash
# usage: export BLACKLIST=/path/to/blacklist; < sample_list.txt xargs -P8 -I{} bash -c 'h5ad_from_archr_annotation.sh {}'

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
md="${md}/ArrowFiles/archr_metadata.tsv.gz"
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

mytempbed=$(mktemp)
bed5k=$(mktemp)
bed500=$(mktemp)
bedtools makewindows -b <(awk '{ print $1"\t0\t"$2 }' < "${refdata}/star/chrNameLength.txt") -w 5000 | grep -Fw -f <(tabix -l "${fragments}") > "${mytempbed}"
Rscript -e "benj::bed.dump(benj::peak_annotation(benj::bed.slurp(\"${mytempbed}\"), \"${refdata}/genes/genes.gtf.gz\"), \"${bed5k}\")"
bedtools makewindows -b <(awk '{ print $1"\t0\t"$2 }' < "${refdata}/star/chrNameLength.txt") -w 500 | grep -Fw -f <(tabix -l "${fragments}") > "${mytempbed}"
Rscript -e "benj::bed.dump(benj::peak_annotation(benj::bed.slurp(\"${mytempbed}\"), \"${refdata}/genes/genes.gtf.gz\"), \"${bed500}\")"
rm "${mytempbed}"


mkdir -p "H5AD/TileMatrix5k"

if [[ -f "${bed5k}" ]]; then
    hdir="H5AD/TileMatrix5k/"
    mkdir -p "${hdir}"
    peaks="${bed5k}"
    if [[ -f "${BLACKLIST}" ]]; then
	h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${peaks}" -o "${hdir}/${sample}.h5ad" -b "${BLACKLIST}"
    else
	h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${peaks}" -o "${hdir}/${sample}.h5ad"
    fi
    gs_in="${hdir}/${sample}.h5ad"
    rm "${bed5k}"
fi

if [[ -f "${bed500}" ]]; then
    hdir="H5AD/TileMatrix500/"
    peaks="${bed500}"
    mkdir -p "${hdir}"
    if [[ -f "${BLACKLIST}" ]]; then
	h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${peaks}" -o "${hdir}/${sample}.h5ad" -b "${BLACKLIST}"
    else
	h5ad_from_archr_annotation.py -f "${fragments}" -s "${sample}" --cell-metadata "${md}" --peaks "${peaks}" -o "${hdir}/${sample}.h5ad"
    fi
    gs_in="${hdir}/${sample}.h5ad"
    rm "${bed500}"
fi

mkdir -p "H5AD/GeneScoreMatrix"
estimate_gene_accessibility -i "${gs_in}" -o "H5AD/GeneScoreMatrix/${sample}.h5ad" --gtf "${gtf}" --tss "${tss}"
