# bbi-genetic-demux

### *work in progress...* 
This is a genetic demultimplexing pipeline intended for BBI's sci-RNA-seq-3 libraries. This pipeline uses the workflow management system, Nextflow, and incorporates the [Souporcell](https://github.com/wheaton5/souporcell) genetic demutiplexing method developed by Haynes Heaton and scripts adapted from the [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Installation.html) platfrom.

## Prerequisites 

1. This pipeline requires Nextflow DSL1 versions >=20.07.1 and <= 22.10.4. 
2. Singularity container
3. Nextflow runs interactively so a terminal multiplexer such as tmux or screen is recommended for the pipeline to continue if you get disconnected from your terminal. 

**Modules** 
1. python 3.12.1
2. samtools/1.19


## Installation

1. Install Nextflow DSL1 through the following command: 
```
export NXF_VER=22.10.0
curl -s https://get.nextflow.io | bash
```

- Check that Nextflow installed correctly and that you have the correct version: 
```
nextflow -v
```

2. Clone the Github repo to a designated directory, for example: 
```
mkdir ~/git
cd ~/git
git clone https://github.com/bbi-lab/bbi-genetic-demux.git
```


## Running the pipeline 

### Input Data
To run the pipeline, you will need a `run_info.csv` sheet and an `experiment.config` sheet. The `run_info.csv` contains the run name and paths for each sample. The `experiment.config` sheet contains various parameters and specifications such as paths to the reference genome fasta/vcf files, and the range of k-values to run genetic demultiplexing on. Refer to the example sheets in the GitHub repository on how to set up these two sheets. 

The sample required inputs are `.bam` and `umis_per_cell_barcode.txt` files generated from the `bbi-sci` pipeline. The pipeline assumes that these files are organized within directories named after each sample and this sample path should be specified in your `run_info.csv` sheet.

Here is an example tree set up: 

```
/path/to/data/
├── Sample1 
│    ├── Sample1.bam
│    └── umis_per_cell_barcode.txt
└──  Sample2
     ├── Sample2.bam
     └── umis_per_cell_barcode.txt
```

A test data set and reference files are provided in the `test_data` directory in the Github repository. If you are running on the test dataset, set the k-value to 2 in the `experiment.config` sheet (eg. `params.k_range= ['2']`) since there are two donors in the test database. We have also provided the SoupOrCell Singularity Image in the `souporcell_image` directory. Set the `params.singularity_path` with the path to `souporcell.sif`. 


### Run Nextflow 

Within your terminal multiplexers and an interactive qlogin session, cd to your current work directory and run the Nextflow pipeline using the following command:

```
nextflow run ~/git/bbi-genetic-demux/main.nf -c experiment.config
```

You can also add the option for Nextflow to output a log file that contains some informative processing information. You can specify the file name. 
```
nextflow run ~/git/bbi-genetic-demux/main.nf -c experiment.config -with-trace genetic_demux_trace.tsv
```

### Output files 

If the pipeline was successful, you should see the following files in your output directory:

```
/path/to/output/
├── alt.mtx
├── ambient_rna.txt
├── cluster_genotypes.vcf
├── clustering.done
├── clusters.err
├── clusters_tmp.tsv
├── clusters.tsv
├── common_variants_covered_tmp.vcf
├── common_variants_covered.vcf
├── consensus.done
├── depth_merged.bed
├── doublets.err
├── ref.mtx
├── remapping.done
├── souporcell_summary.tsv
├── troublet.done
├── variants.done
└── vartrix.done
```
