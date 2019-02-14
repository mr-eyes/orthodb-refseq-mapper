# **Methods**

## Samples of mentioned files


## Essential dictionaries construction

- refseq dictionaries: Parse the `gene2refseq.txt` to construct two dicts:
  - *gene_name_to_ncbigid* : `{geneSymbol:NCBI_ID}`
  - *ncbi_to_refseq* : map `NCBI IDs` to `RefSeq` `transcript ID(s)` , it's a one:many relationship.

- ODB dictionaries: 
  - *odb_ogid_desc* : Mapping `Orthologus Gene Unique ID` to `Gene Description` by parsing the `odb10v0_genes.tab` file
  - *odb_ncbi_ogid* : Mapping `NCBI ID` to `Orthologus Gene Unique ID` as following:
      1. parsing the `odb10v0_genes.tab` file to get `Gene Symbol` and `Orthologus Gene Unique ID`
      2. Quary the *gene_name_to_ncbigid* to get the `NCBI ID`, if not found
      3. deep search option could be selected (Slows down the process as it make online query)
  

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
