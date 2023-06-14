#!/usr/bin/env bash

batchdir=""
atac=true
splice=true
rna=true



# The options in getopt format
OPTIONS=b:aAsSrR
LONGOPTS=batch:,use-atac,no-use-atac,use-splice,no-use-splice,use-rna,no-use-rna

! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    exit 2
fi

eval set -- "$PARSED"

# Now go through all the options with a case and use shift to analyze one option at a time.
while true; do
    case "$1" in
        -b|--batch)
            batchdir="$2"
            shift 2
            ;;
        -a|--use-atac)
            atac=true
            shift
            ;;
        -A|--no-use-atac)
            atac=false
            shift
            ;;
        -s|--use-splice)
            splice=true
            shift
            ;;
        -S|--no-use-splice)
            splice=false
            shift
            ;;
	-r|--use-rna)
	    rna=true
	    shift
	    ;;
	-R|--no-use-rna)
	    rna=false
	    shift
	    ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: Invalid option: $1"
            exit 3
            ;;
    esac
done

# Check if batch was set
if [[ -z "${batchdir}" ]]; then
    echo "Error: --batch or -b is required."
    exit 1
fi

# if [[ -d "$HOME/data-group/${batch}" ]]; then
#     batchdir="$HOME/data-group/${batch}/"
# elif [[ -d "$HOME/data6/${batch}" ]]; then
#     batchdir="$HOME/data6/${batch}/"
# fi

if [[ ! -d "${batchdir}" ]]; then
    echo "Batch dir ${batchdir} does not exist"
fi
function die() {
    if [ "$?" -ne 0 ]; then
        echo "${batchdir}: $1"
        exit 1
    fi
}

function check_h5() {
    h5="$1"
    if [[ "$splice" == "true" ]]; then
	datasets=("layers/ambiguous" "layers/matrix" "layers/spliced" "layers/unspliced" "var/gene_type" "var/pc" "var/interval" "obs/Sample")
    else
	datasets=("var/gene_type" "var/pc" "var/interval" "obs/Sample")
    fi
    for ds in "${datasets[@]}"; do
	h5dump -n "$h5" | grep -q "$ds"
	if [ "$?" -ne 0 ]; then
	    echo "Dataset $ds does not exist in $h5."
	    exit 1
	fi
    done
}
test -f "${batchdir}/sample_list.txt" || die "Sample list does not exist"
test -f "${batchdir}/metrics_summary.tsv" || die "Metrics summary does not exist"
test -f "${batchdir}/aggr.csv" || die "Aggr.csv does not exist"

if [[ "$rna" == "true" ]]; then
    test -d "${batchdir}/H5AD" || die "H5AD dir does not exist"
    test -d "${batchdir}/H5AD/raw" || die "H5AD raw dir does not exist"
    test -d "${batchdir}/H5AD/filtered" || die "H5AD filtered dir does not exist"
    test -d "${batchdir}/H5AD/cellbender" || die "H5AD cellbender dir does not exist"
fi

for sample in $(cut -d , -f 1 < "${batchdir}/aggr.csv" | tail -n+2); do
    if [[ "$rna" == "true" ]]; then
	check_h5 "${batchdir}/H5AD/raw/${sample}.h5ad"
	check_h5 "${batchdir}/H5AD/filtered/${sample}.h5ad"
	check_h5 "${batchdir}/H5AD/cellbender/${sample}.h5ad"
    fi
    if [[ "$atac" == "true" ]]; then
	test -f "${batchdir}/ArrowFiles/${sample}.arrow" || die "Arrow file ${batchdir}/ArrowFiles/${sample}.arrow does not exist"
    fi
done


if [[ "$atac" == "true" ]]; then
    test -f "${batchdir}/ArrowFiles/archr_metadata.tsv.gz" || die "ArchR metadata does not exist"
    zgrep -q '\t' "${batchdir}/ArrowFiles/archr_metadata.tsv.gz" || die "ArchR metadata is not tabbed"
    zcat "${batchdir}/ArrowFiles/archr_metadata.tsv.gz" | head -n1 | grep -Fwq "TSSEnrichment" || die "TSSEnrichment is not in a column in archr metadata"
    tsscol=$(zcat "${batchdir}/ArrowFiles/archr_metadata.tsv.gz" | head -n1 | tr '\t' '\n' | awk '$0 == "TSSEnrichment" { print NR }')
    mintss=$(zcat "${batchdir}/ArrowFiles/archr_metadata.tsv.gz" | cut -f "${tsscol}" | tail -n+2 | sort -g | head -n1)
    tss_lt_11=$(echo "$mintss < 1.1" | bc)
    if [[ $tss_lt_11 -eq 0 ]]; then
	die "Minimum TSS enrichment is ${mintss}"
    fi
fi
exit 0
