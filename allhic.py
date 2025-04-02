#!/usr/bin/env python
import argparse  
import subprocess  
import os  

def run_command(command, step_name):   
    print(f"Running command (Step: {step_name}): {command}")  
    result = subprocess.run(command, shell=True, text=True, capture_output=True)  
    if result.returncode != 0:  
        print(f"Error in {step_name}: {result.stderr.strip()}")  
        raise RuntimeError(f"Command failed: {command}\nError: {result.stderr.strip()}")  
    else:  
        print(f"{step_name} completed successfully.")  

def main(genome, fq1, fq2, threads, m):  
    try:  
        # Step 1: Index genome and create .fai file  
        run_command(f"bwa index {genome}", "Index Genome")  
        run_command(f"samtools faidx {genome}", "Create .fai file")  

        # Step 2: Perform alignment  
        run_command(f"bwa mem -t {threads} {genome} {fq1} {fq2} > bwa.aln.sam", "Alignment")  

        # Step 3: Filter SAM file  
        run_command(f"PreprocessSAMs.pl bwa.aln.sam {genome} MBOI", "Filter SAM File")  

        # Step 4: Convert to BAM and clean  
        run_command(f"samtools view -@ {threads} -bt {genome}.fai bwa.aln.REduced.paired_only.bam > allhic.clean.bam", "Convert SAM to BAM")  

        # Step 5: Partition into groups  
        run_command(f"ALLHiC_partition -b allhic.clean.bam -r {genome} -e GATC -k 11", "Partition")  

        # Step 6: Extract CLM and frequency information  
        run_command(f"allhic extract allhic.clean.bam {genome} --RE GATC", "Extract CLM")  

        # Step 7: Optimize counts  
        for i in range(1, 12):  
            run_command(f"allhic optimize allhic.clean.counts_GATC.11g{i}.txt allhic.clean.clm", f"Optimize Group {i}")  

        # Step 8: Generate final results  
        run_command(f"ALLHiC_build {genome}", "Build Final Results")  

        # Step 9: Generate assembly file from groups.agp  
        run_command(f"python juicebox_scripts/agp2assembly.py groups.agp allhic.assembly", "Generate Assembly File")  
        
        # Step 10: Generate Hi-C file  
        run_command(f"3d-dna/visualize/run-assembly-visualizer.sh -p false allhic.assembly {m}", "Generate Hi-C File")
        
        print("ALLHiC pipeline completed successfully.")  

    except RuntimeError as error:  
        print(f"Pipeline terminated due to an error: {error}")  

if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="Run ALLHiC pipeline.")  
    parser.add_argument('-g', '--genome', required=True, help='Path to the genome file.')  
    parser.add_argument('-fq1', '--fq1', required=True, help='double-ended hic reads 1.')  
    parser.add_argument('-fq2', '--fq2', required=True, help='double-ended hic reads 2.')  
    parser.add_argument('-t', '--threads', type=int, default=64, help='Number of additional threads to use [default: 64]')  
    parser.add_argument('-m', '--merged_nodups', required=True, help='Results of juicer.')
    args = parser.parse_args()  

    # Run the pipeline  
    main(args.genome, args.fq1, args.fq2, args.threads, args.merged_nodups)
