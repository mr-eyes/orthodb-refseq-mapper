# **Methods**

## Essential dictionaries construction

- Parse `gene2refseq.txt` and construct dictionary `{geneSymbol:NCBI_ID}`

- Parse the `odb*_genes.tab` file to construct two dicts:
  - *odb_ncbi_ogid* : Mapping NCBI ID to Orthologus Gene Unique ID
  - *odb_ogid_desc* : Mapping Orthologus Gene Unique ID to Gene Description
  
**Parsing is performed as following:**
> As stated in  [OrthoDB README](https://v100.orthodb.org/download/README.txt) the file  `odb10v0_genes.tab` has column names `[uq_ortho_gene_id, tax_id, prot_seq_id, synonyms(), description]`
So, If the key and value ∈ (columns names list) will just directly parse it from the `odb10v0_genes.tab` file, If ∉ to the columns list, will try to parse it from another file like in the following example.

**dict:** *odb_ncbi_ogid* `NCBI ID` as a key, and to get the NCBI from the `Gene Symbol`:

- Check the corresponding `NCBI ID` to `Gene Symbol`  from the `gene2refseq.txt` file.
    a. If found, just add it to the dictionary and continue.
    b. If the `Gene Symbol` not found in the `gene2refseq.txt` file deep search option could be selected (Slows down the process as it make online query)

-Parsing the `gene2refseq.txt` file to construct a **dict:** *ncbi_to_refseq* that map `NCBI IDs` to `RefSeq` `transcript ID(s)` , it's a one:many relationship.
> Now we have all `NCBI IDs` we extracted from the `odb10v0_genes.tab` file with all it's corresponding `Transcripts IDs` that are found in the refseq `FASTA` file.

## Mapping Heuristics

 Let's short-naming the files, "`odb10v0_genes.tab`" > odb_genes & `{species}_rna.fna.gz` > `refseq_fa`

 First, We check for the relation between the `NCBI IDs` in the `odb_genes` and its corresponding `Transcript IDs` in the `refseq_fa` file.

 The relation could be categorized into one of 3 types, for each `NCBI ID` **:**

**Case 1. [1:1] Relationship**   If the `NCBI ID` occurred in a **single** `odb_genes` record, and it has just **one** corresponding `Transcript ID` in the `refseq_fa`

**Case 2. [1:M] Relationship** If the `NCBI ID` occurred in **single** `odb_genes` record, and it has **multiple** corresponding `Transcript IDs` in the **refseq_fa**

**Case 3. [M:M] Relationship** If the `NCBI ID` occurred in **multiple** `odb_genes` records,   and it has **multiple** corresponding `Transcript IDs` in the **refseq_fa**

### Case 1.
Just map the `Transcript ID` to the `odb gene ID`

### Case 2.

 - Check if there is an "isoform" in any of the `odb_genes` descriptions
	 - If **YES**
		 - Map only the transcripts that are having corressponding `transcript variants`
	 - If **NO**
		 - Map all the transcripts to the same  `odb gene ID`.

### Case 3. (Matching by gene description)
- Iterating over each corressponding transcript header:
	 -  Check if there is an "isoform" in any of the `odb_genes` descriptions
		 - If **YES**
			 - Replace all "transcript variants" with "isoforms" in transcript headers 
		- Iterating over every gene description in all corresponding `og genes IDs`
			- Check if the refseq transcript's gene description exist in the og gene description.
				- If **YES** perform mapping.