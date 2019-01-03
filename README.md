# ATACseqQCnextflow
pipeline for ATACseqQC

## install

### requirements

[nextflow](https://www.nextflow.io/)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[bwa](http://bio-bwa.sourceforge.net/)

[samtools](http://samtools.sourceforge.net/)

[macs2](https://github.com/taoliu/MACS)

[R](https://www.r-project.org/)

[UCSC genome tools](http://hgdownload.soe.ucsc.edu/admin/exe/)

### install

```
git clone https://github.com/jianhong/ATACseqQCnextflow.git
```


## configuration

edit the nextflow.config file depend on your local setting.

## prepare fastq files

put the original fastq files in the folder of ${dataDir}


## run pipeline

```
nextflow run main.nf --dataDir fastq
```

### run it in cluster

change 

`docker.enabled` from `true` to `false`
`singularity.enabled` from `false` to `true`

## output

The output files will be in the folder of results be default. You can change it in the nextflow.config file.

