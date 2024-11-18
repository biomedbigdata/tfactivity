import os
import sys
import argparse
import string
import pandas as pd
from typing import List

parser = argparse.ArgumentParser()

parser.add_argument('--fasta', type=str, required=False, help='Path to the FASTA file')
parser.add_argument('--gtf', type=str, required=False, help='Path to the GTF file')
parser.add_argument('--genome', type=str, required=False, help='genome name/version')
parser.add_argument('--motifs', type=str, required=False, help='Path to JASPAR motif file')
parser.add_argument('--taxon_id', type=str, required=False, help='JASPAR taxon ID for motif download')
parser.add_argument('--outdir', type=str, required=False, help='Path to output directory')

parser.add_argument('--rna_seq', type=str, required=False, help='RNA-Seq directory')
parser.add_argument('--chip_seq', type=str, required=False, help='ChIP-Seq directory')

parser.add_argument('--profile', type=str, required=False, help='nextflow profile', default="apptainer")
parser.add_argument('--apptainer_cache', type=str, required=False, help='Path to the apptainer cache directory')
parser.add_argument('--process_executor', type=str, required=False, help='Executor for the nextflow pipelines', default='local')
parser.add_argument('--process_queue', type=str, required=False, help='scheduler queue')

args = parser.parse_args()

if args.taxon_id is None:
    if args.genome is None or args.motifs is None:
        print('Error: If taxon_id is not specified, both fasta and motifs must be provided.')
        sys.exit(1)

fasta = args.fasta
gtf = args.gtf
genome = args.genome
motifs = args.motifs
taxon_id = args.taxon_id
outdir = args.outdir

rna_seq = args.rna_seq
chip_seq = args.chip_seq

profile = args.profile
apptainer_cache = args.apptainer_cache
process_executor = args.process_executor
process_queue = args.process_queue


os.environ["NXF_APPTAINER_CACHEDIR"] = apptainer_cache

# Create output directory if not exisiting
if not os.path.exists(outdir):
    os.makedirs(outdir)


def list_files(directory, prefix=None, suffix=None, recursive=False) -> List[str]:
    """
    List files in the specified directory, optionally including subdirectories.

    Parameters:
        directory (str): The root directory to start the search.
        prefix (str, optional): Filter files that start with this prefix.
        suffix (str, optional): Filter files that end with this suffix.
        recursive (bool, optional): If True, search subdirectories as well.
                                    If False, search only the specified directory.
                                    Default is False.

    Returns:
        List[str]: List of file paths matching the criteria.
    """

    file_list = []

    if recursive:
        for root, _, files in os.walk(directory):
            for file in files:
                if (prefix is None or file.startswith(prefix)) and (suffix is None or file.endswith(suffix)):
                    file_list.append(os.path.join(root, file))
    else:
        for file in os.listdir(directory):
            full_path = os.path.join(directory, file)
            if os.path.isfile(full_path):
                if (prefix is None or file.startswith(prefix)) and (suffix is None or file.endswith(suffix)):
                    file_list.append(full_path)

    return file_list


# Run nf-core/rnaseq pipeline
path_nfcore_rnaseq = os.path.join(outdir, 'nfcore-rnaseq')
path_outdir_rnaseq = os.path.join(path_nfcore_rnaseq, 'output')
path_samplesheet_rnaseq = os.path.join(path_nfcore_rnaseq, 'samplesheet_rnaseq.csv')

if not os.path.exists(path_nfcore_rnaseq):
    os.makedirs(path_nfcore_rnaseq)

file_paths = list_files(rna_seq, recursive=True)

samples = {}
for file_path in file_paths:
    condition = file_path.split('/')[-2]
    basename = os.path.basename(file_path)

    # Remove .fq.gz and split file name
    sample, rep, read = basename.split('.')[0].split('_')

    sample_name = "_".join([condition, sample, rep])

    if sample_name not in samples:
        samples[sample_name] = {'sample': sample_name,
                                'fastq_1': '',
                                'fastq_2': '',
                                'strandedness': 'auto'}

    if read == 'R1':
        samples[sample_name]['fastq_1'] = file_path
    elif read == 'R2':
        samples[sample_name]['fastq_2'] = file_path

df = pd.DataFrame.from_dict(samples, orient='index').to_csv(path_samplesheet_rnaseq, index=False)

rnaseq_run = f"""
nextflow run \
    nf-core/rnaseq \
    --input {path_samplesheet_rnaseq} \
    --outdir {path_outdir_rnaseq} \
    --gtf {gtf} \
    --fasta {fasta} \
    --igenomes_ignore \
    --genome null \
    -profile {profile} \
    -process.executor {process_executor} \
    -process.queue {process_queue} \
    -resume \
    --skip_deseq2_qc true
"""
os.chdir(path_nfcore_rnaseq)
#os.system(rnaseq_run)


# Run nfcore/chipseq pipeline
path_nfcore_chipseq = os.path.join(outdir, 'nfcore-chipseq')
path_outdir_chipseq = os.path.join(path_nfcore_chipseq, 'output')
path_samplesheet_chipseq = os.path.join(path_nfcore_chipseq, 'samplesheet_chipseq.csv')

if not os.path.exists(path_nfcore_chipseq):
    os.makedirs(path_nfcore_chipseq)

file_paths = list_files(chip_seq, recursive=True)

samples = {}
for file_path in file_paths:
    condition = file_path.split('/')[-3]
    antibody = file_path.split('/')[-2]
    basename = os.path.basename(file_path)

    # Remove .fq.gz and split file name
    sample, rep, read = basename.split('.')[0].split('_')

    condition_antibody_sample = "_".join([condition, antibody, sample])

    sample_rep = '_'.join([condition_antibody_sample, rep])

    # Strip 'REP' from replicate
    if sample_rep not in samples:
        samples[sample_rep] = {'sample': condition_antibody_sample,
                               'fastq_1': '',
                               'fastq_2': '',
                               'replicate': rep.strip(string.ascii_letters),
                               'antibody': antibody,
                               'control': f'{condition}_CONTROL_{sample}',
                               'control_replicate': rep.strip(string.ascii_letters)}

    if read == 'R1':
        samples[sample_rep]['fastq_1'] = file_path
    elif read == 'R2':
        samples[sample_rep]['fastq_2'] = file_path

chipseq_samplesheet = pd.DataFrame.from_dict(samples, orient='index')

# Remove control files that do not exist
chipseq_samplesheet['control'] = chipseq_samplesheet['control'].apply(lambda x: x if x in chipseq_samplesheet['sample'].values else '')

# Remove antibody, control and control_replicate column from dataframe
chipseq_samplesheet.loc[chipseq_samplesheet['sample'].str.contains('CONTROL'), ['antibody', 'control', 'control_replicate']] = ''
chipseq_samplesheet.to_csv(path_samplesheet_chipseq, index=False)

# Runs on revision dev as the pipeline has not had an release for some time
# TODO change parameter read_length
# TODO change parameter skip_preseq
# TODO disable processes not needed for downstream pipelines
chipseq_run = f"""
nextflow run \
    nf-core/chipseq \
    -r dev \
    --input {path_samplesheet_chipseq} \
    --outdir {path_outdir_chipseq} \
    --gtf {gtf} \
    --fasta {fasta} \
    --read_length 50 \
    --save_align_intermeds \
    -profile {profile} \
    -process.executor {process_executor} \
    -process.queue {process_queue} \
    -resume \
    --skip_preseq true
"""

os.chdir(path_nfcore_chipseq)
#os.system(chipseq_run)


# Run nfcore/tfactivity pipeline
path_nfcore_tfactivity = os.path.join(outdir, 'nfcore-tfactivity')
path_outdir_tfactivity = os.path.join(path_nfcore_tfactivity, 'output')
path_samplesheet_tfactivity_peaks = os.path.join(path_nfcore_tfactivity, 'samplesheet_tfactivity_peaks.csv')
path_samplesheet_tfactivity_bams = os.path.join(path_nfcore_tfactivity, 'samplesheet_tfactivity_bams.csv')

path_design_tfactivity_rna = os.path.join(path_nfcore_tfactivity, 'design_tfactivity_rna.csv')
path_counts_tfactivity_rna = os.path.join(path_nfcore_tfactivity, 'counts_tfactivity_rna.csv')

if not os.path.exists(path_nfcore_tfactivity):
    os.makedirs(path_nfcore_tfactivity)

# Create samplesheet peaks
file_paths = list_files(os.path.join(path_outdir_chipseq, 'bwa', 'merged_library', 'macs3', 'broad_peak'), suffix='.broadPeak')

samples = {}
for file_path in file_paths:
    basename = os.path.basename(file_path)
    condition, antibody, sample, rep, _ = basename.split('_')

    sample_id = '_'.join([condition, antibody, sample, rep])

    if sample_id not in samples:
        # TODO: Set 'footprinting' for ChIP-Seq but not for ATAC-Seq and DNase-Seq
        samples[sample_id] = {
            'sample': sample_id,
            'condition': condition,
            'assay': antibody,
            'peak_file': file_path,
        }
    else:
        raise ValueError('Duplicate sample_id detected!')

pd.DataFrame.from_dict(samples, orient='index').sort_values(by='sample').to_csv(path_samplesheet_tfactivity_peaks, index=False)

# Create samplesheet bams
file_paths = list_files(os.path.join(path_outdir_chipseq, 'bwa', 'merged_library'), suffix='.mLb.clN.sorted.bam')

samples = {}
controls = {}
for file_path in file_paths:
    basename = os.path.basename(file_path)
    condition, antibody, sample, rep = basename.split('.')[0].split('_')

    sample_id = '_'.join([condition, antibody, sample, rep])

    if antibody != 'CONTROL' and sample_id not in samples:
        samples[sample_id] = {
            'sample': sample_id,
            'condition': condition,
            'assay': antibody,
            'signal': file_path,
            'merge_key': f'{condition}_{sample}_{rep}',
        }
    elif antibody == 'CONTROL' and sample_id not in controls:
        controls[sample_id] = {
            'merge_key': f'{condition}_{sample}_{rep}',
            'control': file_path,
        }
    else:
        raise ValueError('Duplicated sample_id detected!')

samplesheet_bam = pd.DataFrame.from_dict(samples, orient='index')
controls = pd.DataFrame.from_dict(controls, orient='index')

samplesheet_bam = samplesheet_bam.merge(controls, on='merge_key', how='inner').drop(columns='merge_key').sort_values(by='sample')
samplesheet_bam.to_csv(path_samplesheet_tfactivity_bams, index=False)

# Convert rna counts to csv
# TODO: Remove hardcoded aligner 'star_salmon' when adding aligner choice
counts_rna_tsv = pd.read_csv(os.path.join(path_outdir_rnaseq, 'star_salmon', 'salmon.merged.gene_counts.tsv'), sep='\t')
counts_rna_tsv['gene_id'].to_csv(path_counts_tfactivity_rna, header=False, index=False)

counts_design = {}
for sample in counts_rna_tsv.drop(columns=['gene_id', 'gene_name']).columns:
    condition = sample.split('_')[0]
    path_sample_counts = os.path.join(path_nfcore_tfactivity, f'{sample}_counts.txt')

    counts_rna_tsv[sample].to_csv(path_sample_counts, index=False, header=False)
    counts_design[sample] = {
        'sample': sample,
        'condition': condition,
        'counts_file': path_sample_counts,
    }

pd.DataFrame.from_dict(counts_design, orient='index').to_csv(path_design_tfactivity_rna, index=False)

tfactivity_run = f"""
nextflow run \
    /nfs/home/students/l.hafner/inspect/tfactivity/dev/main.nf \
    nf-core/tfactivity \
    --input {path_samplesheet_tfactivity_peaks} \
    --input_bam {path_samplesheet_tfactivity_bams} \
    --counts {path_counts_tfactivity_rna} \
    --counts_design {path_design_tfactivity_rna} \
    --motifs {motifs} \
    --genome {genome} \
    --fasta {fasta} \
    --gtf {gtf} \
    --outdir {path_outdir_tfactivity} \
    -profile {profile} \
    -process.executor {process_executor} \
    -process.queue {process_queue} \
    -resume
"""

os.chdir(path_nfcore_tfactivity)
os.system(tfactivity_run)
