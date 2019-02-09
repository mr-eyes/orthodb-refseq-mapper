# refseq_orthodb
Mapping RefSeq RNA to OrthoDB

# Run example (Horse refseq):
```bash
python analyse.py ref_seq_data/horse_headers.txt orthodb_data/horse_odb10v0_genes.tab
python check_mapping_quality.py mapped.json orthodb_data/horse_odb10v0_genes.tab
```
