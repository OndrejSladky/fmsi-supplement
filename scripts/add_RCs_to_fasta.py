import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement_fasta(input_file, output_file):
    # Open the output file
    with open(output_file, "w") as output_handle:
        # Parse the input FASTA file
        for record in SeqIO.parse(input_file, "fasta"):
            # Write the original sequence
            SeqIO.write(record, output_handle, "fasta")

            # Calculate the reverse complement of the sequence
            rev_comp_seq = record.seq.reverse_complement()

            # Create a new record for the reverse complement
            rev_comp_record = record
            rev_comp_record.seq = rev_comp_seq
            rev_comp_record.id = record.id + "_RC"
            rev_comp_record.description = "RC of " + record.description

            # Write the reverse complement sequence
            SeqIO.write(rev_comp_record, output_handle, "fasta")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Add reverse complement of each sequence in a FASTA file.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_fasta", help="Path to the output FASTA file")

    # Parse the arguments
    args = parser.parse_args()

    # Run the reverse complement function
    reverse_complement_fasta(args.input_fasta, args.output_fasta)