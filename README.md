# Single Cell 3D UMAP Visualization with Cluster Coloring

This repository provides a script to visualize single-cell RNA-seq data in 3D UMAP space, colored by cluster, using [Scanpy](https://scanpy.readthedocs.io/) and [Rerun](https://www.rerun.io/).

## Features
- Loads a `.h5ad` single-cell dataset
- Computes PCA, neighbors, 3D UMAP, and Leiden clustering if not already present
- Visualizes the 3D UMAP embedding with points colored by cluster in the Rerun viewer
- Easy configuration of data path via `.env` file

## Requirements
- Python 3.8+
- [Scanpy](https://scanpy.readthedocs.io/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [rerun-sdk](https://www.rerun.io/)
- [python-dotenv](https://pypi.org/project/python-dotenv/)
- [matplotlib](https://matplotlib.org/)

Install requirements (if needed):
```bash
pip install scanpy pandas rerun-sdk python-dotenv matplotlib
```

## Setup
1. **Configure your data path:**
   - Copy `.env.example` to `.env` and edit the `H5AD_PATH` variable to point to your `.h5ad` file.
   ```bash
   cp .env.example .env
   # Edit .env to set your H5AD_PATH
   ```

2. **Run the script:**
   ```bash
   python sc3d_cluster_vis.py
   ```
   - The script will launch the Rerun viewer and display the 3D UMAP embedding colored by cluster.

## Notes
- If your `.h5ad` file does not have PCA, neighbors, UMAP, or clustering, the script will compute them automatically.
- The script uses Leiden clustering by default. You can modify the script to use other clustering methods if desired.
- The `.env` file is ignored by git for security and portability. Use `.env.example` as a template.

## Customization
- To color by a specific gene or other metadata, or to add more advanced visualizations, see the comments in `sc3d_cluster_vis.py` or ask for help!

## License
MIT OR Apache-2.0 