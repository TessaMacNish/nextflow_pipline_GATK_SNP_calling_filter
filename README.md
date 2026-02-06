# Nextflow Pipeline for Soft Filtering GATK Called SNPs

## What does this pipeline do?

This pipline filters GATK called SNPs (in VCF format) using ApplyVQSR. The tools used in this pipeline are gatk4 version 4.6.1.0 and bcftools version 1.15. This pipeline submits all individual steps in the pipeline as SLURM scripts.  

## Running this pipeline using only default parameters

1) Set up your directories. In your chosen work directory (the directory you are going to run nextflow from) set up the following directories as shown by the tree command below.
``` bash
tree
.
├── data
├── Genome
|    └── Reference
├── VCF
```

2) In the data/ directory download the txt files with the coordinates for the good and medium quality SNPs. By defualt they are named VQSR_training_good.txt and VQSR_training_medium.txt. An example of what one of these coordinate txt files should look like is below:
``` bash
A01     10001078
A01     10001098
A01     10001332
A01     10001642
A01     10003584
A01     10004280
A01     10004284
A01     10005942
```

3) Download main.nf and nextflow.config into your working directory. Do not change the names of these files.

4) In the data/ directory download your vcfs with the truth sites (known sites) and the corresponding index file.
 ``` bash
Truth_Set.DV10.SNPs.vcf.gz
Truth_Set.DV10.SNPs.vcf.gz.idx
```

5) In the directory Reference download all reference files. You should have the following files:
``` bash
DarmorV10_Chromosomes_Only.dict
DarmorV10_Chromosomes_Only.fa
DarmorV10_Chromosomes_Only.fa.fai
```

6) In the directory VCF/ directory download the raw vcf and corresponding index file from a GATK SNP calling pipeline such as the one available at https://github.com/TessaMacNish/nextflow_pipeline_for_GATK_SNP_calling. You should have the following files:
``` bash
Samples_minicore_1-135_whole_genome.DV10.raw.vcf.gz
Samples_minicore_1-135_whole_genome.DV10.raw.vcf.gz.tbi
```

7) To run this nextflow pipeline you will need sigularity installed. In the nextflow.config file the following line loads singularity.
``` bash
beforeScript = 'module load singularity/4.1.0-nompi'
```
You will need to change this line to the version of singularity available on your CPU. The above singularity module is available on Pawsey super computer's Setonix server. If a singularity module is not available on your CPU you can install singlarity using conda or another package manager. If using a version of singularity isntalled inside a conda environment you can change the module load line to the example below, ensuring you replace myenv with the name of your conda environment.
``` bash
beforeScript = 'conda activate myenv'
```

8) You can now run nextflow from your working directory.
``` bash
nextflow run main.nf
```

You can run nextflow with the additional parameters if you want better logs and traceback, which are useful if something goes wrong.
``` bash
nextflow run main.nf -with-trace -with-report -with-timeline
```

9) nextflow will copy all relevant output into the VCF/ directory.

nextflow will also make a directory called work, which will have all of these results as well as the logs. 

## Running this pipeline defining your own parameters
The purpose of the parameters is to make the pipline more flexible, so that the user can change the input and output file and directory names. All parameters are described below.

`Ref_Abbr` - This controls how multiple output files are named. The default is DV10, which represents an abbreviation of the refernce genome originally used for this pipeline. All output files would then have DV10 in their name. You can change this parameter to be an abbreviation of the reference you use or you can change it to any other identifying string you wish to use for your file names.

`rawVCF` - This is the path to the raw VCF from a GATK SNP calling pipeline. If defining a path to your own raw VCF please use the absolute path to avoid errors.

`referencePrefix` - This is represents how your reference files are named. The default is DarmorV10_Chromosomes_Only but if you had ference files called Reference_1.fa, Reference_1.fa.fai, and Reference_1.dict, you could change this parameter to "Reference_1".

`VCF` - This is the output directory.

`GoodCoord` - The path to the txt file with the coordinates for the good SNPs.

`MediumCoord` - The path to the txt file with the coordinates for the medium SNPs.

`ref` - This is the directory that holds the reference files.

`KnownSites` - The path to your Truth dataset. The defualt is data/Truth_Set.DV10.SNPs.vcf.gz

`outputVCFprefix` - This parameter determines how the soft-filtered final output VCFs are named. 

`tranche` - The tranche level you want to filter the VCF file by. The defualt is 99.9.

An example of how to use the parameters when running your nextflow pipeline is below. All other parameters can be used in the same way.

``` bash
nextflow run main.nf --rawVCF "/path/to/your/rawVCF.vcf.gz" --tranche 99.5
```

Alternatively, you can edit the parameters directly in the main.nf file. 
