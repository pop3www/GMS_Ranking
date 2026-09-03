# Expression input validation

Overall status: **PASS**

Sample sheet: `config/sample_sheet.csv`
Expected samples: 32
Normalize simple repeated non-ENSG IDs: True

## Counts
Path: `data/raw/RawCountFile_rsemgenes.txt`
Status: **PASS**
Genes: 36604
Samples: 32
Sample set matches sample sheet: True
Duplicated genes: 0
Negative values: 0
Fractional values: 39434

## TPM
Path: `data/raw/TPMCountFile_rsemgenes.txt`
Status: **PASS**
Genes: 36604
Samples: 32
Sample set matches sample sheet: True
Duplicated genes: 0
Negative values: 0
Counts/TPM same gene order: True
Counts/TPM same gene set: True
Counts/TPM common genes: 36604

## Notes
- Count matrix contains fractional values. This is compatible with RSEM expected counts; for DESeq2, handle integer conversion explicitly in the DE script and document it.
- Normalized TPM gene IDs: {'LYT2_LYT2': 'LYT2', 'Mcherry_Mcherry': 'Mcherry', 'MYCER-ER_MYCER-ER': 'MYCER-ER'}

