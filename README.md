# Mapping orthoDB to RefSeq RNA

## **Species**

| Name   | dog  | elephant | horse | human | monkey | platypus | rabbit | mouse | cow  |
|--------|------|----------|-------|-------|--------|----------|--------|-------|------|
| TAX_ID | 9615 | 9785     | 9796  | 9606  | 9544   | 9258     | 9986   | 10090 | 9913 |

### Download the required data

```bash
bash prepare_species.sh
```

### Map any species

```bash
python analyse.py <species_name> <species_tax_id>
```

Run flags:

- `--detailed-report` generate detailed report.
- `--deep-search` activate deep search option. (Explained in [Methods](./methods.md))

### Check mapping quality

```bash
python check_mapping_quality.py <species_name> <species_tax_id>
```

> You can add species as much as you can.
>  
> Mapping is not guaranteed to be **PERFECT**.
>  
> Not all transcripts are mapped.
>
> Running with [Pypy](https://pypy.org/) Will speed up the process significantly.
