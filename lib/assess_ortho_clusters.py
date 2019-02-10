import json

map_json_file = "mapped.json"
og2genes_file = "orthodb_data/horse_odb10v0_OG2genes.tab"

# Reading resulted mapping
with open(map_json_file) as f:
    mapped = dict(json.load(f))

# Read OG2Genes File
og2genes = {}
genes2og = {}
gg = []
with open(og2genes_file, 'r') as og:
    for line in og:
        line = line.split()
        gene = line[0]
        og = line[1]

        gg.append(og)

print(len(gg))
print(len(set(gg)))

from collections import Counter

x = dict(Counter(gg))

for k,v in x.items():
    if v == 2:
        print(k)
        break

