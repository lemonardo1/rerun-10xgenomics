import scanpy as sc
import pandas as pd
import rerun as rr

# 1) AnnData 파일 불러오기
adata = sc.read_10x_h5("yourdata.h5")

# 2) UMAP 임베딩이 없다면 계산 (이미 있으면 이 부분 생략)
if "X_umap" not in adata.obsm:
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

# 3) obs와 obsm을 DataFrame으로 변환
metadata = adata.obs.copy().reset_index().rename(columns={"index": "cell_id"})
umap = (
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["umap1", "umap2"],
    )
    .reset_index()
    .rename(columns={"index": "cell_id"})
)
df = metadata.merge(umap, on="cell_id", how="left")

# 4) Rerun SDK 초기화
rr.init("sc_anndata_visualization", spawn=True)

# 5) send_columns로 데이터 전송
columns = rr.AnyValues.columns(
    **{col: df[col].values for col in df.columns}
)
indexes = [rr.AnyValues.index_column("cell_id", df["cell_id"].values)]

rr.send_columns(
    "cells",
    indexes=indexes,
    columns=columns,
) 