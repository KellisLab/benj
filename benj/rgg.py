
def extract_rank_genes_groups(adata, key="rank_genes_groups_filtered", group_name="group"):
    names = pd.DataFrame(adata.uns[key]["names"])
    R, C = np.where((~names.isna()) & (pvals < pval_cutoff))
    df = pd.DataFrame({"names": names.values[R, C], group_name: names.columns.values[C]})
    for subkey in adata.uns[key].keys():
        if type(adata.uns[key][subkey]) == np.recarray:
            xf = pd.DataFrame(adata.uns[key][subkey])
            df[subkey] = xf.values[R, C]
    if "pts" in adata.uns[key].keys():
        idx = adata.uns[key]["pts"].index.get_indexer(df["names"].values)
        df["pts"] = adata.uns[key]["pts"].values[idx, C]
    return df
