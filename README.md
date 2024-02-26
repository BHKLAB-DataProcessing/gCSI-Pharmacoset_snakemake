#OWNER: Jermiah Joseph
#DATE: 2024-02-20
#DESC: A pipeline to process the gCSI dataset from Genentech. At the moment it only does rnaseq data processing.

# Note: 
RNA-seq processing was done on HPC4Health cluster and uploaded to the cloud.


# GR-based drug response dataset from Genentech
[Link to Data Page](http://research-pub.gene.com/gCSI_GRvalues2019/)

- Drugs and Cell Lines per drug
  - http://research-pub.gene.com/gCSI_GRvalues2019/Cell_lines_per_drug.txt
- http://research-pub.gene.com/gCSI_GRvalues2019/gCSI_GRdata_v1.3.rds.tar.gz
- http://research-pub.gene.com/gCSI_GRvalues2019/gCSI_GRdata_v1.3.tsv.tar.gz



conda environment starting off with 
``` bash
mamba create -n gcsi -c conda-forge -c bioconda r-base=4.3.1  r-essentials=4.3.0 r-data.table=1.14.8 bioconductor-CoreGx bioconductor-PharmacoGx bioconductor-annotationdbi -y 

```

Download pset
``` R
PharmacoGx::downloadPSet(name = "gCSI_2019", saveDir = ".", verbose = TRUE, timeout=3600) 
```
