## Reading up on the literature

Preprocessing: SpatialData, AnnData, Scanpy, Sparrow, Harpy.

Cell segmentation: Cellpose.

Sample integration: harmony, scVI, Scanorama, scanpy internal (BBKNN)

Differential expression analysis between clusters: scVI, PyDESeq2

## Animals
Animal 1: A1, A2
Animal 2: B2, C2
Animal 3: D1

## Todo:
V- Check the papers from the hackathon
V- Check papers from spatial transcripomics in hippocampus in pubmed

V-Get the data inside a SpatialData object, with Harpy Help
V-Do image preprocessing 
V- Cell segmentation: Try CellPose with Dask (2000 x 2000 pixel takes 11min without dask) - try to implement it with Dask/MPS
V- Write it to zarr and try to load it again
V- Allocate transcripts
V- Do single-cell analysis
V- do the preprocessing for each sample of ST, and make the AnnData object
    V- transfer data, code and dependencies to HPC
    V- do the actual cell segmentation using cellpose on HPC GPUs interactively, not through the harpy but with cellpose directly. Check the docs.
- Single-cell analysis:

    V- filtering of cells

    V- leiden clustering

    - Integration
        
        V- Test if integration is necessary

        V- Try out Scanorama: https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html

        V- Try ou scVI: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html

        V- Try out BBKNN

        V- Try out ResolVI: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/spatial/resolVI_tutorial.html

        V- Benchmark with scib: https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html
    
    V- tranfer labeling from integrated to the different samples and map back into spatial context

    V- Cell type annotation

        V- manual (based on expression markers and spatial location)

         X- automated: CellTypist: Mouse_Isocortex_Hippocampus or Developing_Mouse_Hippocampus
    
    V- diferential expression: Vcam1, Gpr37l1, Hepacam, Lsamp


V- readup on scVI a bit better, the model training and how do you save the model for example for reproducibility.


Qs:
V- How do I make it aware of the animals/samples? sample per sample clusters, sample_id.

Input from Dani: 
    - Transfer labels with scVI, then do scVI DE with 3 samples/animals.do 
    - Transfer labels with scVI, then do pseudobulk with PyDESeq2 on original raw counts/log counts.
    - Transfer labels with scVI, then do CA2/CA3 astro vs CA2/CA3 neurons.
Input from Joris/Dan:
    - Keep it simple: just astrocytes vs CA1, CA2/3 and DG Subgranular Neurons

Nextup
V- write the zarr with the cell type annotation onto disk and try the DE analysis using scvi. 
V- all DEG analysis in scVI
V- Try out PyDESeq that can account for the animal/sample structure.

Idea:
V- For visualization: tryout a heatmap. Check other papers. With cell_type x genes (with normalized expression)


Next:
V- Option 1: Only show Volcano plot from scVI from Astrocytes vs All Neurons.
V- Option 2: 3 different heatmaps of Astro vs CA1, Astro vs CA2/CA2 and Astro vs DG. Then crossreference them. Further followup on candidates.


V- Idea: segment the hippocampus in CA1, CA3 and DG. Then do DEG analysis of CA1 astro vs CA1 neurons, CA2/3 astro vs CA2/3 neuron, DG astr vs DG Granule Cell/ Subgranular layer
V- Joris meeting: keep it very simple

Further finetuning:
X- stitching of tiles 
V- cell segmentation fine tuning with cellpose params
X- cell segmentation expanding

Notes:
V- Shouldn't include p values from DEG analysis between different groups (because there were selected already based on differences, so you get inflated p-values). So a heatmap of the 
V- Workflow:
    - filter cells/genes
    - concat adatas with sample_id
    - integrate with scVI (raw counts)
    - then pca, neighbors and umap

https://docs.scvi-tools.org/en/stable/tutorials/index_spatial.html
