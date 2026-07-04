# R environment notes

The 125 server did not expose `R` or `Rscript` in the default PATH when this
project was initialized.

## Recommended options

### Option A: use an existing R module/environment

If the server has a module system or an existing conda environment, activate it
before running the benchmark:

```bash
source ~/.bashrc
conda activate <env-with-R-Seurat-CellChat>
cd /home/dzf/cellchat_acceleration
bash code/00_environment_check.sh
```

### Option B: create a conda/mamba environment

If mamba/conda is available:

```bash
mamba env create -f /home/dzf/cellchat_acceleration/environment.yml
conda activate cellchat-acceleration
Rscript -e 'remotes::install_github("sqjin/CellChat")'
Rscript -e 'remotes::install_github("Blake-Deng/CellChatAccelRcpp")'
```

If the acceleration repository is private, install it from the local clone or
configure GitHub credentials first.

### Option C: system R

If you have sudo access:

```bash
sudo apt-get update
sudo apt-get install -y r-base r-base-dev
```

Then install R packages from R/CRAN/Bioconductor/GitHub.

## Required final check

Before any benchmark run:

```bash
bash /home/dzf/cellchat_acceleration/code/00_environment_check.sh
```

