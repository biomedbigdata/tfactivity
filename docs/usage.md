# nf-core/tfactivity: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/tfactivity/usage](https://nf-co.re/tfactivity/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The following inputs can be processed by the pipeline:

- Chromatin accessibility data (at least one of the following)
  - Peaks in BED3-compatible format (e.g. broadPeak), parameter `--input`
  - BAM files, parameter `--input_bam`
- Gene expression data (all of the following)
  - Raw count matrix (e.g. from nf-core/rnaseq), parameter `--counts`
  - Design matrix assigning conditions and optionally batch information to the samples provided in the count matrix, parameter `--counts_design`

The conditions in the peak/BAM samplesheets need to match the conditions in the design matrix.

### Peaks samplesheet

The samplesheet for peak files can look as follows:

```csv title="samplesheet.csv" caption="Required columns for peak files"
sample,condition,assay,peak_file
condition1_H3K27ac_1,condition1,H3K27ac,condition1_H3K27ac_1.broadPeak
condition1_H3K27ac_2,condition1,H3K27ac,condition1_H3K27ac_2.broadPeak
condition1_H3K4me3,condition1,H3K4me3,condition1_H3K4me3.broadPeak
condition2_H3K27ac,condition2,H3K27ac,condition2_H3K27ac.broadPeak
condition3_H3K27ac,condition3,H3K27ac,condition3_H3K27ac.broadPeak
condition3_H3K4me3,condition3,H3K4me3,condition3_H3K4me3.broadPeak
```

:::note
Only the first three columns (chromosome, start, end) of the `bed` format are used.
:::

There are some optional columns which can be added to the samplesheet to configure the footprinting:

```csv title="samplesheet.csv" caption="Optional columns for footprinting"
sample,condition,assay,peak_file,footprinting,include_original,max_peak_gap
condition1_H3K27ac_1,condition1,H3K27ac,condition1_H3K27ac_1.broadPeak,true,true,500
condition1_H3K27ac_2,condition1,H3K27ac,condition1_H3K27ac_2.broadPeak,true,true,500
condition1_H3K4me3,condition1,H3K4me3,condition1_H3K4me3.broadPeak,true,true,500
condition2_H3K27ac,condition2,H3K27ac,condition2_H3K27ac.broadPeak,true,true,500
condition3_H3K27ac,condition3,H3K27ac,condition3_H3K27ac.broadPeak,true,true,500
condition3_H3K4me3,condition3,H3K4me3,condition3_H3K4me3.broadPeak,true,true,500
condition1_ATAC-seq,condition1,ATAC-seq,condition1_ATAC-seq.broadPeak,false,,
```

- `footprinting`: Whether to perform footprinting analysis on the peaks. If enabled, the regions between close peaks will be scanned for transcription factor affinity. This is recommended for Histone modification ChIP-seq data, but not for ATAC-Seq and DNase-Seq data. Default: `true`
- `include_original`: Whether to include the original peaks in the footprinting analysis. Default: `true`
- `max_peak_gap`: Maximum number of base pairs between two peaks to be considered as a single region for footprinting analysis. Default: `500`

### BAM samplesheet

The samplesheet for BAM files can look as follows:

```csv title="samplesheet_bam.csv" caption="Required columns for BAM files"
sample,condition,assay,signal,control
condition1_H3K27ac_1,condition1,H3K27ac,condition1_H3K27ac_1.bam,condition1_control.bam
condition1_H3K27ac_2,condition1,H3K27ac,condition1_H3K27ac_2.bam,condition1_control.bam
condition1_H3K4me3,condition1,H3K4me3,condition1_H3K4me3.bam,condition1_control.bam
condition2_H3K27ac,condition2,H3K27ac,condition2_H3K27ac.bam,condition2_control.bam
condition3_H3K27ac,condition3,H3K27ac,condition3_H3K27ac.bam,condition3_control.bam
condition3_H3K4me3,condition3,H3K4me3,condition3_H3K4me3.bam,condition3_control.bam
```

The first three columns are the same as in the peak file samplesheet. The `signal` column should contain the path to the signal BAM file. The `control` column should contain the path to the control BAM file.

These files are used to predict enhancer regions in the following way:

- Train chromHMM on the signal and control BAM files
- Identify states that are enriched for either `H3K27ac` or `H3K4me3`
- Extract the regions of these states
- Merge close regions to enhancer regions using the ROSE algorithm

The resulting enhancer regions are then used as if they were peak files provided in the peak samplesheet. However, footprinting analysis is not performed on these regions.

### Gene expression data

Gene expression data can be provided in two ways. In both ways, it should be raw counts per gene ID across samples.

1. A single count matrix with gene IDs as rows and samples as columns. This matrix should be provided with the `--counts` parameter. The `--counts_design` parameter is used to match samples in the count matrix to conditions (and optionally batches).
2. A gene list file and one count file per sample. In this case, provide the gene list file with the `--counts` parameter and use the `counts_file` column in `--counts_design` to specify the count files. The files will be merged into a single count matrix (as in the first option) before further processing.

#### Single count matrix

The count matrix (`--counts`) should look like this:

```csv
gene_id,sample1,sample2,sample3
ENSG00000000001,10,20,30
ENSG00000000002,5,10,15
ENSG00000000003,2,4,6
```

The design matrix (`--counts_design`) should look like this:

```csv
sample,condition
sample1,condition1
sample2,condition1
sample3,condition2
```

#### Gene list and multiple count files

The gene list file (`--counts`) should look like this:

```csv
ENSG00000000001
ENSG00000000002
ENSG00000000003
```

The design matrix (`--counts_design`) should look like this:

```csv
sample,condition,counts_file
sample1,condition1,sample1_counts.txt
sample2,condition1,sample2_counts.txt
sample3,condition2,sample3_counts.txt
```

In this case, the count files should look like this:

```
10
20
30
```

:::warning
The number of rows in each count file needs to match the number of rows in the gene list file.
:::

#### Batch effect correction

Optionally, you can specify a column `batch` in the design matrix to correct for batch effects. The batch effect correction is performed using DESeq2. This is possible for both the single count matrix and the gene list with multiple count files.

```csv
sample,condition,batch
sample1,condition1,batch1
sample2,condition1,batch2
sample3,condition2,batch2
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/tfactivity --input ./samplesheet.csv --counts ./count_matrix.csv --counts_design ./counts_design.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/tfactivity -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/tfactivity
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/tfactivity releases page](https://github.com/nf-core/tfactivity/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
