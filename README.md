# Single-Cell-RNA-seq-Benchmark

cellranger count --id=run_count_3kpbmc --fastqs=~/pmbc/pbmc3k_fastqs --transcriptome=./cellranger/refdata-cellranger-GRCh38-3.0.0 --expect-cells=3000

cellranger mat2csv run_c

cellranger mat2csv filtered_feature_bc_matrix pbmc3k.csv

## For seqtk
1. Install
>
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make

