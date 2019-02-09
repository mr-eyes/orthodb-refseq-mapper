import json
import orthodb_class as ortho
import sys

# mapped json file, orthodb source, ref mapping file

# Reading reference mapping

if len(sys.argv) < 3:
    exit("python check_mapping_quality.py <map*json> <orthodb_genes.tab>")


map_json_file = sys.argv[1]
orthodb_file = sys.argv[2]



ncbi2refseq = {}
with open("horse_gen2refseq.txt", 'r') as gen2refseq:
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
ODB = ortho.OrthoDB()
odb_ogid_ncbi = ODB.odb_genes_info(
    path=orthodb_file, tax_id=9796, key="og_id", value="ncbi_gid")

correct = 0
wrong = 0

for tr_id, og_id in mapped.items():
    ncbi_id = odb_ogid_ncbi[og_id][0]
    refseq_ids = ncbi2refseq[ncbi_id]

    if tr_id in refseq_ids:
        correct += 1
    else:
        wrong += 1


print ("Correct: %d\nWrong: %d" % (correct, wrong))
