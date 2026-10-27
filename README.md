# Muscle-Contractility-MYBPC3

## Molecular Mechanisms of Muscle Contractility — Protein Structure and Mutation Analysis

**Student:** Willhelm Wien Galvez  
**Assignment:** #39  
**Gene:** MYBPC3  
**Transcript:** NM_000256.3  
**ClinVar Variation ID:** 42540  
**Variant:** NM_000256.3:c.1504C>T (p.Arg502Trp)  
**Variant type:** Missense  

## Research Question

How does the MYBPC3 gene and the p.Arg502Trp variant relate to cardiac muscle contractility and hypertrophic cardiomyopathy?

## Methods

### 1. Project setup
An RStudio project was created for the MYBPC3 analysis. The project contains the folders `data`, `results`, and `figures`.

### 2. NCBI Gene retrieval
The human MYBPC3 gene was searched and verified using NCBI. The gene summary was saved as `results/gene_summary.txt`.

### 3. Reference protein retrieval
The assigned transcript, NM_000256.3, was used to identify the matching reference protein. The corresponding protein accession was **NP_000247.2**. The protein FASTA sequence was saved in `data/MYBPC3_protein.fasta`.

### 4. Protein sequence analysis
The FASTA sequence was read in R and the frequencies of all 20 standard amino acids were calculated. The results were saved as `results/amino_acid_composition.csv`.

The MYBPC3 reference protein contains **1,274 amino acids**. The three most frequent amino acids were:
- Valine (V): 112
- Alanine (A): 102
- Glutamic acid (E): 100

A bar plot was created as `figures/01_amino_acid_composition.png`.

### 5. ClinVar variant analysis
ClinVar Variation ID 42540 was examined for the assigned variant. The variant is:

**NM_000256.3:c.1504C>T (p.Arg502Trp)**

The nucleotide change replaces C with T and changes arginine (R) at amino-acid position 502 to tryptophan (W). The variant is classified as a missense variant.

### 6. Variant position
A linear protein map was created to show the position of amino acid 502 within the 1,274-amino-acid MYBPC3 protein. The figure was saved as `figures/02_variant_position.png`.

### 7. Key R commands

The analysis was performed in R using the `rentrez` and `ggplot2` packages.

```r
library(rentrez)
library(ggplot2)

gene_symbol <- "MYBPC3"
transcript_accession <- "NM_000256.3"
clinvar_id <- "42540"
variant_protein <- "p.Arg502Trp"
```

The human MYBPC3 gene was searched in NCBI:

```r
gene_hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)
```

The assigned transcript was retrieved from NCBI Nucleotide and the matching protein accession was extracted from the GenBank record:

```r
gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

m <- regexpr('/protein_id="NP_[0-9]+[.][0-9]+"', gb)
protein_accession <- sub(
  '.*"([^"]+)".*',
  '\\1',
  regmatches(gb, m)
)
```

The reference protein FASTA sequence was then retrieved:

```r
protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)
```

The amino-acid sequence was read and the 20 standard amino acids were counted:

```r
fasta_lines <- readLines(fasta_file)
sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]
standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))

results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)
```

The amino-acid composition plot was generated with `ggplot2`:

```r
fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = paste(gene_symbol, "amino-acid composition"),
    x = "Amino acid",
    y = "Number of residues"
  )

ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)
```

The position of the R502W variant was plotted on the 1,274-amino-acid reference protein:

```r
position <- 502

fig2 <- ggplot() +
  geom_segment(
    aes(x = 1, xend = protein_length, y = 0, yend = 0),
    linewidth = 2
  ) +
  geom_point(
    aes(x = position, y = 0),
    size = 4
  ) +
  annotate(
    "text",
    x = position,
    y = 0.25,
    label = variant_protein,
    size = 4
  ) +
  xlim(0, protein_length) +
  ylim(-0.6, 0.6) +
  theme_void() +
  labs(
    title = paste(gene_symbol, "reference protein"),
    x = "Amino-acid position"
  )

ggsave(
  "figures/02_variant_position.png",
  fig2,
  width = 8,
  height = 2.5,
  dpi = 300
)
```

## Results

MYBPC3 encodes cardiac myosin-binding protein C (cMyBP-C), a protein associated with the thick filament of striated muscle. It interacts with components of the contractile machinery and contributes to regulation of cardiac muscle contraction.

The p.Arg502Trp variant changes a conserved arginine residue to the larger amino acid tryptophan. Position 502 is located within the Ig-like C2-type 3 domain of the protein.

ClinVar currently reports an aggregate classification of **Pathogenic/Likely pathogenic**, with multiple submitters contributing evidence and no conflicts in the aggregate classification.


## Interpretation and Discussion

MYBPC3 encodes cardiac myosin-binding protein C (cMyBP-C), a protein associated with the thick filament of striated muscle. It interacts with components of the contractile machinery and contributes to the regulation of cardiac muscle contraction. The reference MYBPC3 protein analyzed in this project is 1,274 amino acids long, and the assigned variant occurs at amino-acid position 502.

The assigned variant, NM_000256.3:c.1504C>T (p.Arg502Trp), is a missense variant in which arginine (R) is replaced by tryptophan (W) at position 502. Position 502 is located within the Ig-like C2-type 3 domain of MYBPC3. The change therefore affects the amino-acid sequence of an existing protein rather than producing a premature stop codon. The amino-acid composition analysis describes the overall composition of the reference protein, but it cannot by itself determine whether Arg502 is functionally important or whether the R502W variant is pathogenic.

Clinical evidence supports an association between R502W and hypertrophic cardiomyopathy. ClinVar reports an aggregate classification of Pathogenic/Likely pathogenic, with multiple submissions contributing evidence and no conflicts in the aggregate classification. The variant has been reported in affected individuals and families, including evidence of disease segregation. This clinical evidence supports an association between the variant and cardiomyopathy, although clinical association alone does not describe the exact molecular mechanism.

Functional evidence provides additional support. ClinVar submissions report studies indicating that the variant alters MYBPC3 protein function and may affect processes relevant to cardiac muscle contraction and relaxation. These findings provide experimental support for a functional effect of the variant. However, computational predictions should be considered separately from experimental observations. Computational tools can support a hypothesis about the effect of an amino-acid substitution, but they cannot independently prove that a variant causes disease.

Overall, the evidence supports a pathway in which the MYBPC3 R502W substitution alters cMyBP-C function, which may affect regulation of cardiac sarcomere contraction and contribute to hypertrophic cardiomyopathy. However, the exact molecular mechanism should be interpreted cautiously because different types of evidence address different biological questions. A useful validation experiment would be to compare engineered cardiac muscle cells carrying R502W with normal control cells and measure contractile force and calcium sensitivity. This would directly test whether the variant produces measurable changes in cardiac muscle contractility.

## Proposed Validation Experiment

A useful experiment would be to compare engineered cardiac muscle cells carrying the R502W variant with normal control cells. Contractile force and calcium sensitivity could be measured to determine whether the variant produces abnormal cardiac muscle contractility.

## Limitations

This analysis is based on reference sequence information, database records, clinical reports, and available functional evidence. Protein sequence composition alone cannot establish pathogenicity. Computational predictions also cannot replace experimental validation.

## Files

- `analysis.R` — R analysis script
- `data/MYBPC3_protein.fasta` — reference protein sequence
- `results/gene_summary.txt` — NCBI gene information
- `results/protein_accession.txt` — reference protein accession
- `results/amino_acid_composition.csv` — amino-acid frequencies
- `figures/01_amino_acid_composition.png` — amino-acid composition plot
- `figures/02_variant_position.png` — position of p.Arg502Trp

## Database Sources

- NCBI Gene
- NCBI Nucleotide
- NCBI Protein
- NCBI ClinVar
- UniProt

**Access date:** October 27, 2026


## References

1. NCBI ClinVar. **MYBPC3 NM_000256.3:c.1504C>T (p.Arg502Trp), Variation ID 42540.** Current aggregate classification: Pathogenic/Likely pathogenic.

2. Carrier, L. et al. (2010). **Short communication: The cardiac myosin binding protein C Arg502Trp mutation: a common cause of hypertrophic cardiomyopathy.** Journal of Molecular and Cellular Cardiology.

3. Zhang, X. et al. (2014). **NMR structure of the C3 domain of human cardiac myosin binding protein-C with the hypertrophic cardiomyopathy-related mutation R502W.** Biochemistry.

4. Mavacamten shows broad benefit in human and mouse models of MYBPC3-related hypertrophic cardiomyopathy.

## Database Links

- NCBI Gene — MYBPC3, Gene ID 4607
- NCBI Nucleotide — NM_000256.3
- NCBI Protein — NP_000247.2
- NCBI ClinVar — Variation ID 42540
- UniProt — Q14896, cardiac myosin-binding protein C

**Database access date:** October 27, 2026
