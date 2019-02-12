from lib import orthodb_class as ortho
import re
import json
import sys
from Bio import SeqIO
from glob import glob

if len(sys.argv) < 2:   
    exit("python analyse.py <species_dir> <taxonomy_id>")

species = sys.argv[1]
taxonomy_id = int(sys.argv[2])
deep_search = 0

if len(sys.argv) == 4:
    if sys.argv == "-d":
        deep_search = True

refseq_file = glob("species/" + species + "/refseq/*gz")[0]
orthodb_file = glob("species/" + species + "/orthodb/" + species + "*_genes.tab")[0]
gene2refseq_file = glob("species/" + species + "/*gene2refseq.txt")[0]
og2genes_file = glob("species/" + species + "/orthodb/" + species + "*_OG2genes.tab")[0]


class Report():
    orthodb_genes_parsing = {}
    def __init__(self):
        pass

REPORT= Report()


def header_tr_line(path):
    import gzip
    res = {}

    if ".gz" in path:
        with gzip.open(path, 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                tr_id = record.id
                line = record.description
                res[tr_id] = line
    else:
        with open(path, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                tr_id = record.id
                line = record.description
                res[tr_id] = line

    return res


def write_mapped_reads(final_map, species, path):
    global og2genes
    mapped_ids = final_map.keys()
    names_file = "species/%s/map/%s_mapped_reads.fa.names" % (species, species)
    mapped_reads_file = "species/%s/map/%s_mapped_reads.fa" % (species, species)
    print("Parsing refseq fasta file...")
    unfound_in_og2genese = 0
    names = open(names_file, "w")

    if ".gz" in path:
        import gzip
        all_seqs = SeqIO.parse(gzip.open(path, 'rt'), "fasta")
    else:
        all_seqs = SeqIO.parse(open(path, 'rt'), "fasta")

    with open(mapped_reads_file, "w") as handle:
        for fasta in all_seqs:
            tr_id = fasta.description.split()[0]
            if tr_id in mapped_ids:
                odb_gene_id = final_map[tr_id]
                if odb_gene_id not in og2genes:
                    unfound_in_og2genese += 1
                    odb_group_ids = "NULL"

                else:
                    odb_group_ids = ";".join(og2genes[odb_gene_id])

                description = tr_id + "|" + odb_gene_id + "|" + odb_group_ids + "|" + species
                names.write(description + "\t" + description + "\n")
                fasta.description = fasta.id = fasta.name = description
                SeqIO.write(fasta, handle, "fasta")

    print("[LOG] %d gene_ids not found in OG2Genes" % (unfound_in_og2genese))


## Parsing OG2Genes
og2genes = {}
with open(og2genes_file, 'r') as f:
    for line in f:
        line = line.strip().split()
        group_id = line[0]
        gene_id = line[1]
        if gene_id in og2genes:
            og2genes[gene_id].append(group_id)
        else:
            og2genes[gene_id] = [group_id]



ODB = ortho.OrthoDB(gene2refseq_file)
ODB.debug = True

if deep_search:
    ODB.activate_deep_search()

odb_ncbi_ogid = ODB.odb_genes_info(path=orthodb_file, tax_id=taxonomy_id, key="ncbi_id", value="odb_gene_id")
odb_ogid_desc = ODB.odb_genes_info(path=orthodb_file, tax_id=taxonomy_id, key="odb_gene_id", value="description")

rfsq_tr_header = header_tr_line(refseq_file)
matched_ncbi_ids = set()


ncbi2refseq = {}
with open(gene2refseq_file, 'r') as gen2refseq:
    ncbi_ids = odb_ncbi_ogid.keys()
    refseq_ids = rfsq_tr_header.keys()

    for line in gen2refseq:
        line = line.strip().split(" ")
        ncbi = line[0]
        refseq = line[1]
        
        if ncbi not in ncbi_ids:
            continue

        if refseq not in refseq_ids:
            continue

        if ncbi not in ncbi2refseq:
            ncbi2refseq[ncbi] = [refseq]
        else:
            ncbi2refseq[ncbi].append(refseq)

final_map = {}

_cases={"1:1":0,"1:m":0,"m:m":0}

#check if the ID exist one time or more
for ncbi_id, refseq_ids in ncbi2refseq.items():
    at_least_one_isoform = False
    max_no_of_separators = 0
    _ncbi_freq = len(odb_ncbi_ogid[ncbi_id])  # How many this NCBI_ID occurred?
    # How many refseq IDs mapped from this NCBI ID?
    _refseqs_freq = len(ncbi2refseq[ncbi_id])

    # Check if the NCBI ID occurred just once.
    if _ncbi_freq == 1:
        if _refseqs_freq == 1:
            # Good case, 1:1
            _cases["1:1"] += 1
            og_id = odb_ncbi_ogid[ncbi_id][0]
            tr_id = refseq_ids[0]

            if tr_id not in final_map:
                final_map[tr_id] = og_id
                matched_ncbi_ids.add(ncbi_id)

        else:
            # Oh! , It has occurred many times in the OrthoDB
            # This is Many to 1 Relationshit
            _cases["1:m"] += 1

            isoform = False
            og_id = odb_ncbi_ogid[ncbi_id][0]
            __desc = odb_ogid_desc[og_id][0].lower()

            # hmmm, OrthoDB may got specific about single isoform, let's check that.
            if "isoform" in __desc:
                isoform = True
                # Last comma to prevent mismatching X1 with X10 ;)
                __desc = "isoform " + __desc.split("isoform")[-1].strip() + ","

            if isoform:
                for _id in refseq_ids:
                    line = rfsq_tr_header[_id]
                    line = line.replace("transcript variant", "isoform")
                    if __desc in line.lower():
                        tr_id = line.split()[0].strip()
                        final_map[tr_id] = og_id
                        matched_ncbi_ids.add(ncbi_id)

            # Not specific? Ok, no problems.
            else:
                for _refseq_id in refseq_ids:
                    final_map[_refseq_id] = og_id
                    matched_ncbi_ids.add(ncbi_id)

    # O_O Many to many relationshittt, let's map with only the gene description.
    else:
        # Many to many
        _cases["m:m"] += 1

        # get all og_ids for the ncbi_id
        og_ids = odb_ncbi_ogid[ncbi_id]
        desc_to_ogid = {}

        for _id in og_ids:
            __desc = odb_ogid_desc[_id][0].lower()
            if "isoform" in __desc:
                # Last comma to prevent mismatching X1 with X10 ;)
                __desc = "isoform " + __desc.split("isoform")[-1].strip() + ","

            desc_to_ogid[__desc] = _id

        refseq_desc = {}

        __isoform = False
        if "isoform" in " ".join(desc_to_ogid.keys()):
            __isoform = True
            at_least_one_isoform = True

        for _id in refseq_ids:
            # Will need more branching
            line = rfsq_tr_header[_id]
            _no_of_separators = len(line.split(","))

            if _no_of_separators > max_no_of_separators:
                max_no_of_separators = _no_of_separators

            if __isoform:
                line = line.replace("transcript variant", "isoform")

            tr_id = line.split()[0].strip()
            
            for __desc, og_id in desc_to_ogid.items():

                if "isoform" in __desc:
                    __desc = "isoform " + __desc.split("isoform")[-1].strip()
                    if __desc in line.lower():
                        final_map[tr_id] = og_id
                        matched_ncbi_ids.add(ncbi_id)
                        

                else:
                    if _no_of_separators < max_no_of_separators:
                        if __desc.lower() in line.lower():
                            final_map[tr_id] = og_id
                            matched_ncbi_ids.add(ncbi_id)
            
                    

## Check how many reads mapped
if len(final_map) != len(rfsq_tr_header):
    print("%d Transcripts matched, %d missing from total %d" %
            (len(final_map), len(rfsq_tr_header) - len(final_map), len(rfsq_tr_header)))
    
    print("%d NCBI IDs mapped from total of %d IDs" % (len(matched_ncbi_ids), len(odb_ncbi_ogid.keys())))
else:
    print("All records matched for the NCBI ID: %s" % (ncbi_id))

output_json = "species/%s/map/%s_" % (species, species)
f = open(output_json + "mapped.json", "w")
f.write(json.dumps(final_map, indent=4, sort_keys=True))
f.close()


## Extract matched reads
print("Writing mapped reads fasta and names files..")
write_mapped_reads(final_map, species, refseq_file)




## Print the unmatched to a separate file: unmatched_transcripts.txt
# unmatched = set(rfsq_tr_header.keys()) - set(final_map.keys())

# with open("unmatched_transcripts.txt", 'w') as unmatched_file:
#     for _tr_id in unmatched:
#         unmatched_file.write(rfsq_tr_header[_tr_id] + "\n")


print ("-----------")
print (_cases)
