# MutFinder by Joshua Bugis
# Developed for the Barkai Lab at the Weizmann Institute of Science
# Created as part of Gábor Szabó's python course
# 03/02/2025

# Import libraries
import pandas as pd
import argparse
import os
import re
import requests
import gget
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


# Parse arguments
def parse_arguments():
    """
    Parse command-line arguments for the MutFinder script.
    """
    parser = argparse.ArgumentParser(description="MutFinder: Detect mutations in genomic sequences.")
    parser.add_argument("-seqfile", required=True, help="Path to the FASTA file containing the sequence to analyze.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-gene", help="Gene symbol (e.g., FOXP2).")
    group.add_argument("-ensemblid", help="Ensembl gene or transcript ID (e.g., ENSG00000128573).")
    parser.add_argument("-species", help="Species name (e.g., Homo sapiens). Required if -gene is provided.")
    parser.add_argument("--noncanonical", action="store_true", help="Expand analysis to the noncanonical isoform.")
    return parser.parse_args()


# Read FASTA sequence
def read_fasta_sequence(file_path):
    """
    Read the first nucleotide sequence from a FASTA file.
    """
    for seq_record in SeqIO.parse(file_path, "fasta"):
        return str(seq_record.seq)  # Return the first sequence as a string
    raise ValueError("The FASTA file is empty.")  # Raise an error if no sequences are found


# Extract the ensembl gene ID for your reference sequence
def get_ensembl_gene_id(gene=None, species=None, ensemblid=None):
    """
    Retrieve the Ensembl gene ID from either a gene symbol (with species) or directly via an ensemblid.
    """
    if gene:  # If the gene name is given:
        if not species:
            raise ValueError("Species must be provided with gene symbol.")
        response = requests.get(
            f"https://rest.ensembl.org/lookup/symbol/{species}/{gene}?content-type=application/json"
        )
        if response.ok:
            data = response.json()
            return data['id']
        else:
            raise ValueError("Error: Could not find relevant Ensembl ID associated with the gene and species provided.")
    elif ensemblid:  # If the ensembl ID is directly given:
        return ensemblid
    else:
        raise ValueError("Either gene or ensemblid must be provided.")


# Parse the reference data from gget into a workable format
def parse_ref_data(ref_data):
    """
    Convert raw FASTA output from gget into a dictionary mapping transcript IDs to sequences.
    """
    if isinstance(ref_data, list):
        sequences = {}
        header = None
        for line in ref_data:
            if line.startswith(">"):
                header = line.strip()[1:].split()[0]
                sequences[header] = ""
            elif header:
                sequences[header] += line.strip()
        return sequences
    elif isinstance(ref_data, dict):
        return ref_data
    else:
        return {"canonical": ref_data}


# Get the reference transcript sequence based on ensembl id
def retrieve_reference_sequence_by_ensembl(ensembl_gene_id, noncanonical=False):
    if noncanonical:
        ref_data = gget.seq(ensembl_gene_id, isoforms=True)
        if not ref_data:
            raise ValueError("No valid reference sequence was retrieved for noncanonical analysis.")
        return parse_ref_data(ref_data)
    else:
        gene_info = gget.info(ensembl_gene_id)
        canonical_transcript = gene_info.get("canonical_transcript")
        if isinstance(canonical_transcript, (pd.Series,)):
            canonical_transcript = canonical_transcript.iloc[0]
        if not canonical_transcript:
            raise ValueError("Canonical transcript not found in gene info.")
        ref_seq = gget.seq(canonical_transcript, isoforms=False)
        if not ref_seq:
            raise ValueError("No valid canonical reference sequence was retrieved.")
        if isinstance(ref_seq, (list, dict)):
            sequences = parse_ref_data(ref_seq)
            return sequences.get(canonical_transcript, max(sequences.values(), key=len))
        return ref_seq


# Find the open reading frame of the input and reference sequences
def find_longest_orf(sequence):
    """
    Find the longest open reading frame (ORF) in a nucleotide sequence.
    """
    longest_orf = ""
    seq_len = len(sequence)
    stop_codons = {"TAA", "TAG", "TGA"}
    
    # Loop over each possible start codon.
    for i in range(seq_len - 2):
        if sequence[i:i+3] == "ATG":
            found_stop = False
            # Iterate codon by codon starting after this ATG.
            for j in range(i + 3, seq_len, 3):
                codon = sequence[j:j+3]
                if len(codon) < 3:  # Reached the end with a partial codon.
                    break
                if codon in stop_codons:
                    candidate_orf = sequence[i:j+3]
                    found_stop = True
                    if len(candidate_orf) > len(longest_orf):
                        longest_orf = candidate_orf
                    break  # Stop at the first in-frame stop codon.
            if not found_stop:
                # If no stop codon is found, assume one immediately after the sequence.
                candidate_orf = sequence[i:]
                if len(candidate_orf) > len(longest_orf):
                    longest_orf = candidate_orf
    return longest_orf


# Align input and reference sequences
def align_sequences(seq1, seq2, display = True):
    """
    Align two sequences using global alignment and return the best alignment.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(seq1, seq2)
    if display:
        print(f"\n\n\n{str(alignments[0])}\n\n\n")
    return alignments[0]  # Return the best alignment (first alignment has highest score)


# Extract alignment sequences with dashes
def extract_sequences_from_alignment(alignment):
    """
    Extract and stitch together the target and query sequences from a formatted alignment string.
    """
    alignment_string = str(alignment)
    lines = alignment_string.strip().split('\n')
    input_sequence = []
    ref_sequence = []
    
    # Loop through the lines and extract sequences based on line number
    for i, line in enumerate(lines):
        # Process each line and classify into input or reference sequence:
        if not line.strip() or ' ' not in line:
            continue
        sequence_part = line.split(' ', 1)[1].strip()
        sequence_part = re.sub(r'\d+|\s+', '', sequence_part)
        if i % 4 == 0:
            input_sequence.append(sequence_part)
        elif i % 4 == 2:
            ref_sequence.append(sequence_part)
    input_sequence_str = ''.join(input_sequence)
    ref_sequence_str = ''.join(ref_sequence)

    return input_sequence_str, ref_sequence_str


# Find the mutations in the input compared to the reference
def get_mutations(s1, s2):
    """
    Compare two aligned sequences to detect mutations
    """
    differences = []
    i = j = 0
    while i < len(s1) or j < len(s2):
        # If both indices are valid and the characters match, simply advance.
        if i < len(s1) and j < len(s2) and s1[i] == s2[j]:
            i, j = i + 1, j + 1
            continue
        pos = j
        # Get the current characters (ignoring gaps)
        char_input = s1[i] if i < len(s1) and s1[i] != '-' else ''
        char_ref   = s2[j] if j < len(s2) and s2[j] != '-' else ''
        # Check for substitutions
        if char_input and not char_ref and i + 1 < len(s1) and j + 1 < len(s2):
            next_char_input = s1[i + 1] if s1[i + 1] != '-' else ''
            next_char_ref   = s2[j + 1] if s2[j + 1] != '-' else ''
            if not next_char_input and next_char_ref:
                differences.append((pos, 'sub', f"{char_input}->{next_char_ref}"))
                i, j = i + 2, j + 2
                continue
        # Determine the mutation type.
        if char_input and not char_ref:
            differences.append((pos, 'ins', char_input))
        elif char_ref and not char_input:
            differences.append((pos, 'del', char_ref))
        else:
            differences.append((pos, 'sub', f"{char_input}->{char_ref}"))
        i, j = i + 1, j + 1
    return differences


# Report adjacent indels as the same mutation
def merge_adjacent_mutations(mutations):
    """
    Merge adjacent mutations of the same type.
    """
    if not mutations:
        return []
    merged = [mutations[0]]
    for pos, mtype, seq in mutations[1:]:
        last_pos, last_mtype, last_seq = merged[-1]
        # Check if the current mutation is adjacent to the last one
        if mtype == last_mtype and pos == last_pos + len(last_seq):
            # Merge the mutation by appending the new sequence
            merged[-1] = (last_pos, mtype, last_seq + seq)
        else:
            merged.append((pos, mtype, seq))
    return merged


# Create output csv
def save_results_to_csv(output_data, output_file="mutfinder_results.csv"):
    """
    Save mutation analysis results to a CSV file.
    """
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(output_file, mode="a", header=not os.path.exists(output_file), index=False)
    print(f"Results written to {output_file}")


# Combine together the functions above to create a mutational analysis pipeline
def mutational_analysis(original_nuc_seq, input_longest_rf, input_aa_seq, transcript_id, ref_seq, gene_or_id):
    """
    Mutational analysis between the input sequence and a reference transcript.
    """
    # Determine the reference reading frame and translate to amino acids
    reference_longest_rf = find_longest_orf(ref_seq)
    reference_aa_seq = str(Seq(reference_longest_rf).translate())
    
    # Nucleotide mutational analysis
    nt_alignment = align_sequences(input_longest_rf, reference_longest_rf, display=True)
    nt_aligned_input, nt_aligned_ref = extract_sequences_from_alignment(nt_alignment)
    nt_mutations = get_mutations(nt_aligned_input, nt_aligned_ref)
    nt_mutations = merge_adjacent_mutations(nt_mutations)
    has_mutations = bool(nt_mutations)
    alignment_warning = ""
    if nt_aligned_input.count('-') > 25 or nt_aligned_ref.count('-') > 25:
        alignment_warning = ("Caution! Percentage of alignment is lower than expected for the input and reference sequences. Ensure that you have chosen the correct reference sequence, or try another isoform.")
    
    # Amino acid mutational analysis
    aa_alignment = align_sequences(input_aa_seq, reference_aa_seq, display=False)
    aa_aligned_input, aa_aligned_ref = extract_sequences_from_alignment(aa_alignment)
    aa_mutations = get_mutations(aa_aligned_input, aa_aligned_ref)
    aa_mutations = merge_adjacent_mutations(aa_mutations)
    
    # Prepare data for the output csv
    output_data = {
        "Search Date and Time": pd.Timestamp.now(),
        "Input Gene/ID": gene_or_id,
        "Input Nucleotide Sequence": original_nuc_seq,
        "Translated Input Amino Acid Sequence": input_aa_seq,
        "Reference Transcript ID": transcript_id,
        "Reference Nucleotide Sequence": ref_seq,
        "Translated Reference Amino Acid Sequence": reference_aa_seq,
        "Input Contains Mutations": has_mutations,
        "List of Nucleotide Mutations": ", ".join([f"{pos}-{change_type}-{seq}" for pos, change_type, seq in nt_mutations]),
        "List of Amino Acid Mutations": ", ".join([f"{pos}-{change_type}-{seq}" for pos, change_type, seq in aa_mutations]),
        "Alignment Warning": alignment_warning
    }
    save_results_to_csv(output_data)


# MAIN
def main():
    args = parse_arguments()

    input_nucleotide_seq = read_fasta_sequence(args.seqfile)

    ensembl_gene_id = get_ensembl_gene_id(gene=args.gene, species=args.species, ensemblid=args.ensemblid)
    if args.noncanonical:
        reference_data = retrieve_reference_sequence_by_ensembl(ensembl_gene_id, noncanonical=True)
    else:
        reference_data = retrieve_reference_sequence_by_ensembl(ensembl_gene_id, noncanonical=False)
    
    input_longest_rf = find_longest_orf(input_nucleotide_seq)
    input_aa_seq = str(Seq(input_longest_rf).translate())
    
    if args.noncanonical:
        transcripts = reference_data 
    else:
        transcripts = {"canonical": reference_data}
    
    for transcript_id, ref_seq in transcripts.items():
        mutational_analysis(
            original_nuc_seq=input_nucleotide_seq,
            input_longest_rf=input_longest_rf,
            input_aa_seq=input_aa_seq,
            transcript_id=transcript_id,
            ref_seq=ref_seq,
            gene_or_id=args.gene or args.ensemblid
        )


if __name__ == "__main__":
    main()