# Nanocompore pipeline
This repo contains the code for running the [Nanocompore](https://github.com/tleonardi/nancompore) analysis workflow using [Nextflow](https://www.nextflow.io/).

## Pipeline steps
The pipeline runs the following steps:

* Basecalls the raw fast5 files using Albacore
* Run [pycoQC](https://github.com/a-slide/pycoQC) on Albacore's output
* Prepare a fasta file for the reference transcriptome
* Map the basecalled data to the reference using [minimap2](https://github.com/lh3/minimap2)
* Realign the raw signal-level data to the kmers of the reference with [Nanopolish](https://github.com/jts/nanopolish)
* Collapse Nanopolish output by kmer using [NanopolishComp](https://github.com/a-slide/NanopolishComp)
* Run Nanocompore (to be implemented)

All steps are executed in a Conda environment that includes all dependencies.

## How to run

### Sample annotation
Prepare a tab-separated file that describes the samples:
```
SampleName  Condition DataPath
Ctrl1       Ctrl      /path/to/fast5dir
Ctrl2       Ctrl      /path/to/fast5dir
Treated1    Trt       /path/to/fast5dir
Treated2    Trt       /path/to/fast5dir
```

### Configure the pipeline
Configure the pipeline by editing the _nextflow.config_ file.

All parameters are described in the comments and should be self explanatory.

The only two options that deserve an explanation are `target_trancripts` and `input_is_basecalled`.
The former allows you to provide the path to a text file that lists Ensembl transcript IDs of interest. 
Any transcript not present in this list will be discarded from the reference.
`input_is_basecalled` allows you to start the pipeline __after__ the Albacore step. In order for this to work,
the paths in the sample annotation file must point to basecalled fast5 files.

### Run the pipeline
To run the pipeline just execute:
`nextflow run pipeline.nf`

### Profiles
The pipeline is shipped with 2 profiles: local (default) and yoda, which respectively run the processing on the local
machine or on the Yoda HPC cluster at the EBI-EMBL.
The local profile limits the number of cuncurrent tasks to 10, which is reduced to 2 and 3 for Albacore and Nanopolish respectively.
The Yoda profile if configured to use the LSF executor with dynamic memory management.
