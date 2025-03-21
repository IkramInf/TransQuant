# TransQuant

## Overview
TransQuant is a user-friendly RNA-Seq analysis pipeline for transcript quantification. This tool provides a streamlined workflow for processing RNA-Seq data, from quality control through alignment to gene/transcript counting. TransQuant is designed to be accessible to researchers with limited bioinformatics experience while providing the flexibility needed for more advanced users.

## Features
- Quality control analysis using FastQC
- Read trimming with Trimmomatic
- Genome indexing and read alignment with HISAT2
- BAM file processing with Samtools
- Transcript quantification with HTSeq-count
- Support for both single-end and paired-end sequencing data
- Detailed logging and error reporting
- Ability to resume analysis from any step
- Configurable multithreading support

## Requirements
TransQuant requires the following software to be installed and available in your PATH:

- Python 3.6 or higher
- FastQC
- Trimmomatic
- HISAT2
- Samtools
- HTSeq

Additionally, you'll need:
- Reference genome in FASTA format
- Gene annotations in GTF format
- Illumina adapter files (included in Trimmomatic installation)

## Installation

### Option 1: Clone the repository
```bash
git clone https://github.com/IkramInf/TransQuant.git
cd TransQuant
```

### Option 2: Install dependencies with conda
```bash
conda create -n transquant python=3.8 fastqc trimmomatic hisat2 samtools htseq
conda activate transquant
```

## Usage

### Basic usage:

For single-end sequencing:
```bash
python transquant.py -m SE -r1 sample_R1.fastq -r reference.fasta -g annotation.gtf
```

For paired-end sequencing:
```bash
python transquant.py -m PE -r1 sample_R1.fastq -r2 sample_R2.fastq -r reference.fasta -g annotation.gtf
```

### All available options:

```
usage: transquant.py [-h] -m {SE,PE} -r1 R1_FASTQ -r REFERENCE_GENOME -g GTF
                     [-r2 R2_FASTQ] [-o OUTPUT_DIR] [-t THREADS]
                     [--skip_trim] [--skip_fastqc] [--skip_index]
                     [--resume {index,fastqc,trim,align,bam,count}] [-v]

RNA-Seq Analysis Pipeline

Required arguments:
  -m {SE,PE}, --mode {SE,PE}
                        Mode of sequencing: SE (single-end) or PE (paired-end)
  -r1 R1_FASTQ, --r1_fastq R1_FASTQ
                        Path to the R1 FASTQ file (forward reads)
  -r REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        Path to the reference genome FASTA file
  -g GTF, --gtf GTF     Path to the GTF annotation file

Optional arguments:
  -r2 R2_FASTQ, --r2_fastq R2_FASTQ
                        Path to the R2 FASTQ file (reverse reads, required for PE mode)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to the output directory (default: rna_seq_output)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 4)
  --skip_trim           Skip the trimming step
  --skip_fastqc         Skip the FastQC step
  --skip_index          Skip the genome indexing step if index exists
  --resume {index,fastqc,trim,align,bam,count}
                        Resume the pipeline from a specific step
  -v, --verbose         Enable verbose logging
  -h, --help            Show this help message and exit
```

## Pipeline Steps

TransQuant performs the following steps:

1. **Genome Indexing** - Creates a HISAT2 index from the reference genome
2. **Quality Control** - Runs FastQC on raw reads to assess quality
3. **Read Trimming** - Trims adapter sequences and low-quality bases using Trimmomatic
4. **Alignment** - Maps reads to the reference genome using HISAT2
5. **BAM Processing** - Converts, sorts, and indexes the alignment files using Samtools
6. **Read Counting** - Quantifies gene/transcript expression using HTSeq-count

## Output Structure

TransQuant generates the following directory structure:

```
output_dir/
├── fastqc/              # FastQC reports for raw and trimmed reads
├── trimmed/             # Trimmed FASTQ files
├── genome_index/        # HISAT2 genome index files
├── aligned_reads.sam    # Raw alignment file
├── aligned_reads.bam    # Binary alignment file
├── aligned_reads_sorted.bam  # Sorted alignment file
├── aligned_reads_sorted.bam.bai  # BAM index
├── counts.txt           # Raw HTSeq count output
├── final_counts.txt     # Filtered count data (without __no_feature, etc.)
├── alignment_summary.txt  # HISAT2 alignment statistics
├── flagstat.txt         # Samtools alignment statistics
└── pipeline.log         # Detailed log of the entire analysis
```

## Example Workflow

### 1. Prepare your reference genome and annotation

Download your reference genome and annotation files from a source like Ensembl or UCSC.

### 2. Run the complete pipeline

```bash
python transquant.py -m PE \
    -r1 sample_1_R1.fastq.gz \
    -r2 sample_1_R2.fastq.gz \
    -r reference/genome.fa \
    -g reference/genes.gtf \
    -o results/sample_1 \
    -t 8
```

### 3. Run multiple samples

```bash
for sample in sample_1 sample_2 sample_3; do
    python transquant.py -m PE \
        -r1 ${sample}_R1.fastq.gz \
        -r2 ${sample}_R2.fastq.gz \
        -r reference/genome.fa \
        -g reference/genes.gtf \
        -o results/${sample} \
        -t 8 \
        --skip_index
done
```

## Troubleshooting

**Issue**: Pipeline fails at trimming step
- **Solution**: Ensure that adapter files (TruSeq3-PE-2.fa for paired-end or TruSeq3-SE.fa for single-end) are available in your working directory or specify their full path.

**Issue**: "Command not found" errors
- **Solution**: Make sure all required tools are installed and in your PATH.

**Issue**: High memory usage during indexing
- **Solution**: For large genomes, run the indexing step on a high-memory machine or increase your system's swap space.

## Citation

If you use TransQuant in your research, please cite it as:

```
TransQuant: A user-friendly RNA-Seq analysis pipeline. https://github.com/yourusername/TransQuant
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Contact

Email - ikraminf.dev@gmail.com
Project Link: [https://github.com/IkramInf/TransQuant](https://github.com/IkramInf/TransQuant)
