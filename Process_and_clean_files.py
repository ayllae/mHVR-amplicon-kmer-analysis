# This script: Searches a folder for all .fastq files; Reads sequences from each FASTQ file
#Filters reads by sequence length (keeps only reads between 200 and 500 bp);
#Renames all reads into a standardized format: sampleName_Read1, sampleName_Read2, etc.;
#Merges all filtered reads from all samples into a single FASTA output file
#Prints how many reads were kept and filtered for each sample.

from Bio import SeqIO  
import os 
import glob  

# folder containing FASTQ files
input_dir = "Raw_fastq_renamed_dogs"

# full path to merged FASTA output
output_fasta = "all_sequences.fasta"

# minimum read size allowed
min_length = 200

# maximum read size allowed
max_length = 500


def seq_dir_to_fasta():

    # find all .fastq files in input folder
    seq_files = glob.glob(
        os.path.join(input_dir, "*.fastq")
    )

    # print files detected
    print("Sequence files found:")

    # loop through files and print names
    for f in sorted(seq_files):
        print(" ", f)

    # open final merged fasta output file
    with open(output_fasta, "w") as out_f:

        # process one FASTQ file at a time
        for filepath in sorted(seq_files):

            # use filename as sample name
            sample_name = os.path.basename(filepath).split(".")[0]

            # counter to rename reads sequentially
            read_counter = 1

            # count reads that passed filter
            kept_reads = 0

            # count reads removed by filter
            filtered_reads = 0

            # open current FASTQ file
            with open(filepath, "r") as handle:

                # read sequences from FASTQ
                for record in SeqIO.parse(handle, "fastq"):

                    # get read sequence length
                    seq_len = len(record.seq)

                    # skip reads outside size range
                    if seq_len < min_length or seq_len > max_length:
                        filtered_reads += 1
                        continue

                    # rename read ID into standardized format
                    record.id = (
                        f"{sample_name}_Read{read_counter}"
                    )

                    # remove extra FASTQ description text
                    record.description = ""

                    # write read into merged FASTA file
                    SeqIO.write(record, out_f, "fasta")

                    # increase read counter
                    read_counter += 1

                    # increase passed reads counter
                    kept_reads += 1

            # print summary for current sample
            print(
                f"{sample_name}: "
                f"{kept_reads} sequences written | "
                f"{filtered_reads} filtered by length"
            )

    # print final output file name
    print(f"\nFASTA written to: {output_fasta}")


# run function
seq_dir_to_fasta()
