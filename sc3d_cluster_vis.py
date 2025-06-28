import scanpy as sc
import pandas as pd
import numpy as np
import rerun as rr

# 1. Load AnnData
adata = sc.read_h5ad("yourdata.h5ad")

# 2. Preprocessing (if needed)
if 'X_pca' not in adata.obsm:
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)

# 3. Compute 3D UMAP if not present
if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[1] != 3:
    sc.tl.umap(adata, n_components=3)

# 4. Compute Leiden clustering if not present
if 'leiden' not in adata.obs.columns:
    sc.tl.leiden(adata)

# 5. Prepare DataFrame for visualization
umap3d = pd.DataFrame(
    adata.obsm['X_umap'],
    index=adata.obs_names,
    columns=['umap1', 'umap2', 'umap3']
).reset_index().rename(columns={'index': 'cell_id'})

metadata = adata.obs.copy().reset_index().rename(columns={'index': 'cell_id'})
df = metadata.merge(umap3d, on='cell_id', how='left')

# 6. Assign colors to clusters
def get_cluster_colors(clusters):
    import matplotlib.pyplot as plt
    unique_clusters = sorted(clusters.unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    cmap = plt.get_cmap('tab20', len(unique_clusters))
    color_map = {cl: tuple((np.array(cmap(i)[:3]) * 255).astype(np.uint8)) for i, cl in enumerate(unique_clusters)}
    # Convert to a numpy array of shape (n_cells, 3)
    return np.array([color_map[cl] for cl in clusters])

colors = get_cluster_colors(df['leiden'])

# 7. Rerun visualization
rr.init("sc3d_cluster_vis", spawn=True)
rr.log(
    "cells/umap3d",
    rr.Points3D(
        positions=df[['umap1', 'umap2', 'umap3']].values,
        colors=colors,
        labels=df['leiden'].astype(str).values
    )
)
print("3D UMAP visualization with cluster coloring sent to Rerun.") 