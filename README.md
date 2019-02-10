# refseq_orthodb

Mapping RefSeq RNA to OrthoDB

 
# Run example:

**Species**

| Name   | dog  | elephant | horse | human | monkey | platypus | rabbit | mouse | cow  |
|--------|------|----------|-------|-------|--------|----------|--------|-------|------|
| TAX_ID | 9615 | 9785     | 9796  | 9606  | 9544   | 9258     | 9986   | 10090 | 9913 |

### Download the required data
```bash
bash prepare_species.sh
```
### Process any species
```bash
python analyse.py <species_name> <species_tax_id>
```

> You can add species as much as you can.
 
> Mapping is not guaranteed to be **PERFECT**.

> Not all reads are mapped.
> 
