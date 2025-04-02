import os
import subprocess
import argparse

def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def main(genome, threads, curatedlib, deepte_dir, plants_model_dir):
    # Step 1: RepeatModeler
    run_command(f"BuildDatabase -name {genome} -engine rmblast {genome}")
    run_command(f"RepeatModeler -database {genome} -threads {threads}")
    
    # Step 2: EDTA
    run_command(f"EDTA.pl --genome {genome} --species others --sensitive 1 --step all --anno 1 -t {threads} --force 1 --curatedlib {curatedlib}")

    # Step 3: LTR classification
    run_command("grep 'LTR/unknown' *.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq *.mod.EDTA.TElib.fa - > LTR_unknown.fa")
    run_command("grep -v 'LTR/unknown' *.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq *.mod.EDTA.TElib.fa - > LTR_known.fa")
    run_command("conda activate DeepTE && "
                f"python {deepte_dir}/DeepTE.py -i LTR_unknown.fa -sp P -m_dir {plants_model_dir} -fam LTR")
    
    # Step 4: Post-processing DeepTE output
    run_command("sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta | "
                "sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' | "
                "sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' > LTR_unknown_DeepTE.fa")
    
    # Step 5: Merge repeat libraries
    run_command("cat LTR_unknown_DeepTE.fa LTR_known.fa > EDTA.TElib.fa")
    run_command("cat EDTA.TElib.fa *-families.fa > repeat.lib.fa")
    
    # Step 6: RepeatMasker
    run_command(f"RepeatMasker {genome} -lib repeat.lib.fa -poly -html -gff -pa {threads} -nolow -no_is -norna")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for repeat annotation using RepeatModeler, EDTA, DeepTE and RepeatMasker")
    parser.add_argument("-g", "--genome", required=True, help="Path to the genome file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("-l", "--curatedlib", required=True, help="Path to curatedlib. EDTA provided rice, maize and Arabidopsis in https://github.com/oushujun/EDTA/tree/master/database")
    parser.add_argument("-d", "--deepte_dir", required=True, help="Path to DeepTE installation directory")
    parser.add_argument("-p", "--plants_model_dir", required=True, help="Path to Plants model directory for DeepTE")

    args = parser.parse_args()
    
    main(args.genome, args.threads, args.curatedlib, args.deepte_dir, args.plants_model_dir)
