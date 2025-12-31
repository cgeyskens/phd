### Reading up on the literature

Packages: SpatialData, AnnData, Scanpy, Sparrow, Harpy.

Cell segmentation: Cellpose, Segger, Baysor.

# Todo:
V- Check the papers from the hackathon
V- Check papers from spatial transcripomics in hippocampus in pubmed

V-Get the data inside a SpatialData object, with Harpy Help
V-Do image preprocessing 
V- Cell segmentation: Try CellPose with Dask (2000 x 2000 pixel takes 11min without dask) - try to implement it with Dask/MPS
V- Write it to zarr and try to load it again
V- Allocate transcripts
- Do single-cell analysis

Further refinement:
- stitching of tiles 
- cell segmentation fine tuning with cellpose params
- cell segmentation expanding

Nextup:
- sample integration
