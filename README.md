# mutfinder

**mutfinder** is a Python-based tool designed for genomic research to identify and classify mutations in DNA sequences. By comparing user-provided sequences against reference sequences retrieved from genomic databases, **mutfinder** provides detailed reports on mutations, including their impact on amino acid sequences.

---

## How does **mutfinder** work?

- Accepts input sequences in FASTA format.
- Retrieves reference sequences based on:
  - Gene symbol and species name.
  - Ensembl gene or transcript ID.
- Automatically identifies the longest reading frame for accurate translation.
- Translates nucleotide sequences into amino acids using BioPython's `translate` method.
- Generates a CSV file with detailed information on mutations, including:
  - Mutation type classification (silent/non-silent).
  - Full sequence comparison results.
  - Session-wide cumulative outputs.

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/b-joshua/mutfinder.git
   cd mutfinder
   ```

2. Install dependencies using `requirements.txt`:
   ```bash
   pip install -r requirements.txt
   ```

### Dependencies:

- [BioPython](https://biopython.org/docs/dev/api/Bio.Seq.html)
- [gget](https://gget.readthedocs.io/)
- [Pandas](https://pandas.pydata.org/)

---

## How to use **mutfinder**

Run **mutfinder** from the command line using the following syntax:

```bash
python mutfinder.py -seqfile <sequence_file.fasta> [(-gene <gene_symbol> -species <species>) | -ensemblid <ensembl_id>] [--canonical]
```

### Parameters:
- `-seqfile` (required): Path to the FASTA file containing the sequence to analyze.
- Either:
  - `-gene` and `-species`: Specify the gene symbol (e.g., FOXP2) and species name (e.g., Homo sapiens).
  - `-ensemblid`: Provide either an Ensembl gene ID (e.g., ENSG00000128573) or transcript ID (e.g., ENST00000408937.3).
- `--canonical` (optional): Restricts analysis to the canonical isoform of the provided gene.

**Note:** You must provide either `-gene` and `-species` together or `-ensemblid`. Both cannot be provided simultaneously.

---

## Output

**mutfinder** generates a CSV file with the following columns:

1. **Search date and time**: Timestamp of the analysis.
2. **Input gene name/ID**: Identifier for the input gene.
3. **Input nucleotide sequence**: Sequence from the input FASTA file.
4. **Translated input amino acid sequence**: Protein translation of the input sequence.
5. **Reference sequence ID**: Identifier for the retrieved reference sequence.
6. **Reference nucleotide sequence**: Sequence of the reference.
7. **Translated reference amino acid sequence**: Protein translation of the reference sequence.
8. **Input contains mutations**: Boolean indicating if mutations exist.
9. **Change in amino acid sequence**: Boolean indicating if the amino acid sequence differs between input and reference.
10. **List of all mutations**: Mutations in the format: `M30T(pos.204C>T), T45A(pos.265A>G)`
11. **List of silent mutations**: Silent mutations in the format: `MutationType(pos.PositionBaseChange)`
12. **List of non-silent mutations**: Non-silent mutations in the format: `MutationType(pos.PositionBaseChange)`

---

## Example

### Command:

```bash
python mutfinder.py -seqfile sample.fasta -gene FOXP2 -species Homo_sapiens
```

### Sample Output CSV:

| Date & Time       | Gene/ID   | Input Seq. | Trans. Input Seq. | Ref. ID | Ref. Seq. | Trans. Ref. Seq. | Mutations | AA Changes | All Mutations                    | Silent Mutations | Non-Silent Mutations             |
|--------------------|-----------|------------|--------------------|---------|-----------|------------------|-----------|------------|----------------------------------|------------------|----------------------------------|
| 2025-01-10 12:34  | FOXP2     | ATG...TAA  | M...*              | ENST... | ATG...TAA | M...*            | True      | True       | M30T(pos.204C>T), T45A(pos.265A>G) | None             | M30T(pos.204C>T), T45A(pos.265A>G) |

---

## Citation

If you use **mutfinder**, please cite the following works for the tools and databases it leverages:

- **gget**:  
  Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. *Bioinformatics.* https://doi.org/10.1093/bioinformatics/btac836  
- **Ensembl**:  
  Martin FJ et al., Ensembl 2023. *Nucleic Acids Res.* 2023 Jan 6;51(D1):D933-D941. https://doi.org/10.1093/nar/gkac958  
- **UniProt**:  
  The UniProt Consortium, UniProt: the Universal Protein Knowledgebase in 2023, *Nucleic Acids Res.* 2023 Jan 6;51(D1):D523–D531. https://doi.org/10.1093/nar/gkac1052

