import json
from lib import orthodb_class as ortho
import sys
from glob import glob

if len(sys.argv) < 2:   
    exit("python check_mapping_quality.py <species_dir> <taxonomy_id>")

species = sys.argv[1]
taxonomy_id = int(sys.argv[2])

orthodb_file = glob("species/" + species + "/orthodb/" + species + "*_genes.tab")[0]
gene2refseq_file = glob("species/" + species + "/*gene2refseq.txt")[0]
map_json_file = glob("species/" + species + "/map/*json")[0]



# mapped json file, orthodb source, ref mapping file

# Reading reference mapping

ncbi2refseq = {}
with open(gene2refseq_file, 'r') as gen2refseq:
    for line in gen2refseq:
        line = line.strip().split(" ")
        ncbi = line[0]
        refseq = line[1]

        if ncbi not in ncbi2refseq:
            ncbi2refseq[ncbi] = [refseq]
        else:
            ncbi2refseq[ncbi].append(refseq)


# Reading resulted mapping
with open(map_json_file) as f:
    mapped = dict(json.load(f))

# Reading OrthoDB og_id to NCBI
ODB = ortho.OrthoDB(gene2refseq_file)
odb_ogid_ncbi = ODB.odb_genes_info(path=orthodb_file, tax_id=9796, key="og_id", value="ncbi_gid")

correct = 0
wrong = 0

for tr_id, og_id in mapped.items():
    tr_id = str(tr_id)
    og_id = str(og_id)

    ncbi_id = odb_ogid_ncbi[og_id][0]
    refseq_ids = ncbi2refseq[ncbi_id]

    if tr_id in refseq_ids:
        correct += 1
    else:
        wrong += 1


print ("Correct: %d\nWrong: %d" % (correct, wrong))
