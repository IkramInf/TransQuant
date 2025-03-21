#!/usr/bin/env python3
"""
RNA-Seq Analysis Pipeline

A user-friendly pipeline for RNA-Seq data analysis that includes:
- Quality control with FastQC
- Read trimming with Trimmomatic
- Read alignment with HISAT2
- BAM file processing with Samtools
- Read counting with HTSeq
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
import time


class RNASeqPipeline:
    """RNA-Seq analysis pipeline manager"""
    
    def __init__(self, args):
        """Initialize the pipeline with command line arguments"""
        self.mode = args.mode
        self.r1_fastq = args.r1_fastq
        self.r2_fastq = args.r2_fastq
        self.ref_genome = args.reference_genome
        self.gtf_file = args.gtf
        self.output_dir = args.output_dir
        self.threads = str(args.threads)
        self.skip_trim = args.skip_trim
        self.skip_fastqc = args.skip_fastqc
        self.skip_index = args.skip_index
        self.resume = args.resume
        self.verbose = args.verbose
        
        # Set up logging
        log_level = logging.DEBUG if self.verbose else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(os.path.join(self.output_dir, "pipeline.log")),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger("RNA-Seq Pipeline")
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Set up file paths
        self.setup_file_paths()
        
    def setup_file_paths(self):
        """Set up file paths for pipeline outputs"""
        # Index files
        self.index_dir = os.path.join(self.output_dir, "genome_index")
        self.index_prefix = self.index_dir + "/genome"
        
        # Alignment files
        self.output_sam = os.path.join(self.output_dir, "aligned_reads.sam")
        self.output_bam = os.path.join(self.output_dir, "aligned_reads.bam")
        self.sorted_bam = os.path.join(self.output_dir, "aligned_reads_sorted.bam")
        
        # Count files
        self.count_file = os.path.join(self.output_dir, "counts.txt")
        self.final_count_file = os.path.join(self.output_dir, "final_counts.txt")
        
        # Trimmed files paths will be set during trimming
        self.r1_trim_path = None
        self.r2_trim_path = None
        
        # FastQC output directory
        self.fastqc_dir = os.path.join(self.output_dir, "fastqc")
        os.makedirs(self.fastqc_dir, exist_ok=True)
        
        # Trimmed output directory
        self.trim_dir = os.path.join(self.output_dir, "trimmed")
        os.makedirs(self.trim_dir, exist_ok=True)
        
        # Index directory
        os.makedirs(self.index_dir, exist_ok=True)
        
    def run_command(self, cmd, description, output_file=None):
        """Run a shell command with logging and error handling"""
        self.logger.info(f"Starting: {description}")
        self.logger.debug(f"Command: {' '.join(cmd)}")
        
        start_time = time.time()
        
        try:
            if output_file:
                with open(output_file, 'w') as out_f:
                    result = subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE)
            else:
                result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
            elapsed_time = time.time() - start_time
            self.logger.info(f"Completed: {description} in {elapsed_time:.2f} seconds")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error in {description}: {e}")
            self.logger.error(f"Error output: {e.stderr.decode()}")
            return False
            
    def check_file_exists(self, filepath, description):
        """Check if a file exists and log the result"""
        if os.path.exists(filepath):
            self.logger.info(f"{description} exists: {filepath}")
            return True
        else:
            self.logger.warning(f"{description} does not exist: {filepath}")
            return False
            
    def build_hisat2_index(self):
        """Build HISAT2 index for the reference genome"""
        # Check if index already exists
        index_files = [f"{self.index_prefix}.{i+1}.ht2" for i in range(8)]
        if all(os.path.exists(file) for file in index_files) and self.skip_index:
            self.logger.info("HISAT2 index already exists. Skipping index building.")
            return True
        
        cmd = ["hisat2-build", "-p", self.threads, self.ref_genome, self.index_prefix]
        return self.run_command(cmd, "Building HISAT2 index")
        
    def run_fastqc(self):
        """Run FastQC on input FASTQ files"""
        if self.skip_fastqc:
            self.logger.info("Skipping FastQC as requested")
            return True
            
        fastq_files = [self.r1_fastq]
        if self.mode == "PE":
            fastq_files.append(self.r2_fastq)
            
        cmd = ["fastqc"] + fastq_files + ["-o", self.fastqc_dir, "-t", self.threads]
        return self.run_command(cmd, "Running FastQC on input reads")
        
    def run_trimmomatic(self):
        """Run Trimmomatic for read trimming"""
        if self.skip_trim:
            self.logger.info("Skipping trimming as requested")
            self.r1_trim_path = self.r1_fastq
            if self.mode == "PE":
                self.r2_trim_path = self.r2_fastq
            return True
            
        if self.mode == "PE":
            # For paired-end data
            fname1 = os.path.basename(self.r1_fastq).split(".")[0]
            fname2 = os.path.basename(self.r2_fastq).split(".")[0]
            
            self.r1_trim_path = os.path.join(self.trim_dir, f"{fname1}_trimmed.fastq")
            r1_unpaired_path = os.path.join(self.trim_dir, f"{fname1}_unpaired.fastq")
            self.r2_trim_path = os.path.join(self.trim_dir, f"{fname2}_trimmed.fastq")
            r2_unpaired_path = os.path.join(self.trim_dir, f"{fname2}_unpaired.fastq")
            
            adapter_file = "TruSeq3-PE-2.fa"
            
            cmd = [
                "trimmomatic", "PE", "-threads", self.threads, "-phred33",
                self.r1_fastq, self.r2_fastq,
                self.r1_trim_path, r1_unpaired_path,
                self.r2_trim_path, r2_unpaired_path,
                f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3",
                "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
            
        else:
            # For single-end data
            fname = os.path.basename(self.r1_fastq).split(".")[0]
            self.r1_trim_path = os.path.join(self.trim_dir, f"{fname}_trimmed.fastq")
            
            adapter_file = "TruSeq3-SE.fa"
            
            cmd = [
                "trimmomatic", "SE", "-threads", self.threads, "-phred33",
                self.r1_fastq, self.r1_trim_path,
                f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3",
                "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
            
        success = self.run_command(cmd, "Trimming reads with Trimmomatic")
        
        # Run FastQC on trimmed reads if trimming was successful
        if success and not self.skip_fastqc:
            fastq_files = [self.r1_trim_path]
            if self.mode == "PE":
                fastq_files.append(self.r2_trim_path)
                
            cmd = ["fastqc"] + fastq_files + ["-o", self.fastqc_dir, "-t", self.threads]
            self.run_command(cmd, "Running FastQC on trimmed reads")
            
        return success
        
    def align_reads(self):
        """Align reads to the reference genome using HISAT2"""
        if self.mode == "PE":
            cmd = [
                "hisat2", "-p", self.threads, "-x", self.index_prefix,
                "-1", self.r1_trim_path, "-2", self.r2_trim_path,
                "-S", self.output_sam, "--summary-file", 
                os.path.join(self.output_dir, "alignment_summary.txt")
            ]
        else:
            cmd = [
                "hisat2", "-p", self.threads, "-x", self.index_prefix,
                "-U", self.r1_trim_path, "-S", self.output_sam,
                "--summary-file", os.path.join(self.output_dir, "alignment_summary.txt")
            ]
            
        return self.run_command(cmd, "Aligning reads with HISAT2")
        
    def process_bam(self):
        """Convert SAM to BAM, sort and index"""
        # Convert SAM to BAM
        cmd1 = ["samtools", "view", "-@ ", self.threads, "-bS", self.output_sam, "-o", self.output_bam]
        success1 = self.run_command(cmd1, "Converting SAM to BAM")
        
        if not success1:
            return False
            
        # Sort BAM file
        cmd2 = ["samtools", "sort", "-@ ", self.threads, self.output_bam, "-o", self.sorted_bam]
        success2 = self.run_command(cmd2, "Sorting BAM file")
        
        if not success2:
            return False
            
        # Index sorted BAM file
        cmd3 = ["samtools", "index", self.sorted_bam]
        success3 = self.run_command(cmd3, "Indexing BAM file")
        
        # Generate alignment stats
        cmd4 = ["samtools", "flagstat", self.sorted_bam]
        flagstat_file = os.path.join(self.output_dir, "flagstat.txt")
        success4 = self.run_command(cmd4, "Generating alignment statistics", flagstat_file)
        
        return success3 and success4
        
    def count_reads(self):
        """Count reads using HTSeq-count"""
        cmd = [
            'htseq-count', '-m', 'union', '-f', 'bam', 
            '--additional-attr=transcript_id', '-s', 'yes', 
            self.sorted_bam, self.gtf_file
        ]
        
        success = self.run_command(cmd, "Counting reads with HTSeq", self.count_file)
        
        if success:
            # Remove __no_feature, __ambiguous, etc. lines
            with open(self.count_file, 'r') as input_file:
                with open(self.final_count_file, 'w') as output_file:
                    for line in input_file:
                        if not line.startswith('__'):
                            output_file.write(line)
            
            self.logger.info(f"Final counts written to {self.final_count_file}")
            
        return success
        
    def run_pipeline(self):
        """Run the complete RNA-Seq analysis pipeline"""
        self.logger.info("=" * 50)
        self.logger.info("Starting RNA-Seq Analysis Pipeline")
        self.logger.info("=" * 50)
        
        # Check if we're resuming from a specific step
        steps = ["index", "fastqc", "trim", "align", "bam", "count"]
        start_idx = 0
        
        if self.resume and self.resume in steps:
            start_idx = steps.index(self.resume)
            self.logger.info(f"Resuming pipeline from {self.resume} step")
        
        # Track overall success
        all_success = True
        
        # Step 1: Build HISAT2 index
        if start_idx <= steps.index("index"):
            success = self.build_hisat2_index()
            if not success and not self.skip_index:
                self.logger.error("Failed to build HISAT2 index. Pipeline halted.")
                return False
        
        # Step 2: Run FastQC
        if start_idx <= steps.index("fastqc"):
            success = self.run_fastqc()
            all_success = all_success and (success or self.skip_fastqc)
        
        # Step 3: Run Trimmomatic
        if start_idx <= steps.index("trim"):
            success = self.run_trimmomatic()
            if not success and not self.skip_trim:
                self.logger.error("Failed to trim reads. Pipeline halted.")
                return False
        
        # Step 4: Align reads
        if start_idx <= steps.index("align"):
            success = self.align_reads()
            if not success:
                self.logger.error("Failed to align reads. Pipeline halted.")
                return False
        
        # Step 5: Process BAM files
        if start_idx <= steps.index("bam"):
            success = self.process_bam()
            if not success:
                self.logger.error("Failed to process BAM files. Pipeline halted.")
                return False
        
        # Step 6: Count reads
        if start_idx <= steps.index("count"):
            success = self.count_reads()
            all_success = all_success and success
        
        if all_success:
            self.logger.info("=" * 50)
            self.logger.info("RNA-Seq Analysis Pipeline completed successfully!")
            self.logger.info("=" * 50)
            self.logger.info(f"Final counts are available in: {self.final_count_file}")
            self.logger.info(f"Full logs available in: {os.path.join(self.output_dir, 'pipeline.log')}")
            return True
        else:
            self.logger.error("RNA-Seq Analysis Pipeline completed with errors.")
            return False


def main():
    """Parse command line arguments and run the pipeline"""
    parser = argparse.ArgumentParser(
        description='RNA-Seq Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-m', '--mode', type=str, choices=['SE', 'PE'], required=True,
                        help='Mode of sequencing: SE (single-end) or PE (paired-end)')
    required.add_argument('-r1', '--r1_fastq', type=str, required=True,
                        help='Path to the R1 FASTQ file (forward reads)')
    required.add_argument('-r', '--reference_genome', type=str, required=True,
                        help='Path to the reference genome FASTA file')
    required.add_argument('-g', '--gtf', type=str, required=True,
                        help='Path to the GTF annotation file')
    
    # Optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-r2', '--r2_fastq', type=str,
                        help='Path to the R2 FASTQ file (reverse reads, required for PE mode)')
    optional.add_argument('-o', '--output_dir', type=str, default='rna_seq_output',
                        help='Path to the output directory')
    optional.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use')
    optional.add_argument('--skip_trim', action='store_true',
                        help='Skip the trimming step')
    optional.add_argument('--skip_fastqc', action='store_true',
                        help='Skip the FastQC step')
    optional.add_argument('--skip_index', action='store_true',
                        help='Skip the genome indexing step if index exists')
    optional.add_argument('--resume', type=str, choices=['index', 'fastqc', 'trim', 'align', 'bam', 'count'],
                        help='Resume the pipeline from a specific step')
    optional.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.mode == 'PE' and not args.r2_fastq:
        parser.error("Paired-end (PE) mode requires --r2_fastq")
        
    for file_path in [args.r1_fastq, args.reference_genome, args.gtf]:
        if not os.path.exists(file_path):
            parser.error(f"File does not exist: {file_path}")
            
    if args.mode == 'PE' and not os.path.exists(args.r2_fastq):
        parser.error(f"R2 FASTQ file does not exist: {args.r2_fastq}")
    
    # Run the pipeline
    pipeline = RNASeqPipeline(args)
    success = pipeline.run_pipeline()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
