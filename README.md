# Host Gene Classifier

`Classification_genes.py` classifies host genes from a CSV file using:

- MyGene.info annotations
- Gene Ontology (GO) terms
- deterministic fallback rules based on gene symbols and biotypes

The annotation file remains an annotation-only cache. By default, the classified data is written back into the input CSV with added host-gene columns.

## What it does

Given an input CSV and the column containing host gene Ensembl IDs, the script:

1. Reads the input table.
2. Normalizes host gene Ensembl IDs.
3. Uses an annotation file next to the script by default:
   `<input_stem>.mygene_annotations.csv`
4. Reuses that annotation file when it exists and downloads only missing genes from MyGene.info.
5. Downloads `go-basic.obo` if needed.
6. Classifies each unique host gene into a broad functional category.
7. Appends classification columns back onto the input CSV.
8. Overwrites the input CSV by default, or writes to a separate file when `--output` is provided.

## Classification categories

- `Ribosomal protein`
- `Ribosome biogenesis`
- `Translation`
- `RNA processing`
- `Transcription`
- `Genome maintenance`
- `Transport`
- `Development`
- `Signaling`
- `Other coding`
- `lncRNA`
- `Intergenic`

## Requirements

- Python 3
- `pandas`
- `goatools`
- internet access when downloading:
  - MyGene.info annotations
  - Gene Ontology `go-basic.obo`

Install dependencies with:

```bash
pip install pandas goatools
```

## Input

Minimum required input:

- a CSV file
- a column containing host gene Ensembl IDs

Example input:

```csv
Host_Gene
ENSG00000142599
ENSG00000048707
ENSG00000116138
```

## Annotation file

By default the script uses this annotation file path:

```text
<script_dir>/<input_stem>.mygene_annotations.csv
```

That file is reused across runs. If it already contains some genes, the script downloads only the missing ones and saves the updated annotation file back to the same path.

If you provide `--annotation-file`, it should follow the same MyGene-style structure used by the script output cache, with a recognizable gene ID column such as `gene_id`, `query`, or `_id`.

Use `--force-download` to ignore the existing annotation file contents and rebuild annotations for the genes present in the input CSV.

## Usage

Basic usage:

```bash
python3 Classification_genes.py input.csv --host-gene-column Host_Gene
```

Use a specific annotation file:

```bash
python3 Classification_genes.py input.csv \
  --host-gene-column Host_Gene \
  --annotation-file custom_annotations.csv
```

Force a fresh annotation download:

```bash
python3 Classification_genes.py input.csv \
  --host-gene-column Host_Gene \
  --force-download
```

Write to a specific output file:

```bash
python3 Classification_genes.py input.csv \
  --host-gene-column Host_Gene \
  --output results.csv
```

## Command-line options

```text
usage: Classification_genes.py [-h] --host-gene-column HOST_GENE_COLUMN
                               [--output OUTPUT]
                               [--annotation-file ANNOTATION_FILE]
                               [--force-download]
                               [--refresh-go]
                               [--email EMAIL]
                               input_csv
```

Main options:

- `--host-gene-column`: column in the input CSV containing host gene Ensembl IDs
- `--output`: output CSV path; when omitted, the input CSV is overwritten in place
- `--annotation-file`: annotation CSV to reuse and update
- `--force-download`: ignore existing annotations in the annotation file and download them again
- `--refresh-go`: re-download `go-basic.obo`
- `--email`: optional email to pass to MyGene.info

## Output

By default, the script overwrites the input CSV in place.

If `--output` is provided, the script writes exactly to that path.

Example:

```text
python3 Classification_genes.py input.csv --host-gene-column Host_Gene --output results.csv
```

The output keeps the original input rows and appends host-gene columns such as:

- `host_gene_id_resolved`
- `host_gene_symbol`
- `host_gene_description`
- `host_gene_biotype`
- `host_gene_taxid`
- `host_gene_annotation_notfound`
- `host_gene_type`
- `host_gene_type_rule`
- `host_gene_go_bp_ids`
- `host_gene_go_bp_terms`
- `host_gene_go_mf_ids`
- `host_gene_go_mf_terms`
- `host_gene_go_cc_ids`
- `host_gene_go_cc_terms`

The annotation cache file remains separate and contains annotation-only fields.

## Classification logic summary

The classifier applies these rules in broad order:

1. `Intergenic` for empty or explicit intergenic entries
2. `lncRNA` when the gene biotype matches long non-coding RNA biotypes
3. `Ribosomal protein` using strict ribosomal GO or naming rules
4. GO-based assignment using Biological Process and Molecular Function anchor terms
5. Symbol-based fallback rules when GO evidence is missing or uninformative
6. `Other coding` if no stronger signal is found

`host_gene_type_rule` records the rule that produced the final class.

## Notes

- The script expects Ensembl-style host gene IDs.
- Version suffixes such as `.1`, `.2`, and similar are stripped during normalization.
- Duplicate host genes are classified once and then merged back into the input table.
- If MyGene.info cannot resolve a gene, the output keeps the row and marks `host_gene_annotation_notfound = True`.
