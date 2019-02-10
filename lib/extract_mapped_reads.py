import json
from Bio import SeqIO


map_json_file = "mapped.json"
fasta_file = "ref_seq_data/GCF_002863925.1_EquCab3.0_rna.fna"
output_file = "mapped_reads.fa"
names_file = "mapped_reads.fa.names"

# Reading resulted mapping
with open(map_json_file) as f:
    mapped = dict(json.load(f))

mapped_ids = map(str, list(mapped.keys()))

all_sequences = SeqIO.parse(open(fasta_file), 'fasta')
names = open(names_file, "w")
with open(output_file, "w") as handle:
    for fasta in all_sequences:
        tr_id = fasta.description.split()[0]
        if tr_id in mapped_ids:
            description = tr_id
            names.write(description + "\t" + description + "\n")
            fasta.description = fasta.id = fasta.name = description
            SeqIO.write(fasta, handle, "fasta")

        
names.close()
