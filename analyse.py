import orthodb_class as ortho
import refseq
from tqdm import tqdm
import re
import json
import sys
from tqdm import tqdm

refseq_file = "test_cases/refseq.txt"
orthodb_file = "test_cases/odb.tsv"

# refseq_file = "ref_seq_data/horse_headers.txt"
# orthodb_file = "orthodb_data/horse_odb10v0_genes.tab"


def header_tr_geneName(path):
    res = {}
    with open(path, 'r') as refseq:
        for line in refseq:
            tr_id = line.split(" ")[0].replace(">", "")
            gene_name = re.sub('^.*\((.*?)\)[^\(]*$', '\g<1>', line)
            res[tr_id] = gene_name

    return res


def header_tr_line(path):
    res = {}
    with open(path, 'r') as refseq:
        for line in refseq:
            tr_id = line.split(" ")[0].replace(">", "")
            res[tr_id] = line.strip()

    return res


ODB = ortho.OrthoDB()
odb_ncbi_ogid = ODB.odb_genes_info(path=orthodb_file, tax_id=9796, key="ncbi_gid", value="og_id")
odb_ogid_desc = ODB.odb_genes_info(path=orthodb_file, tax_id=9796, key="og_id", value="description")

rfsq_tr_geneName = header_tr_geneName(refseq_file)
rfsq_tr_header = header_tr_line(refseq_file)

ncbi2refseq = {}
with open("horse_gen2refseq.txt", 'r') as gen2refseq:
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

        else:
            # Oh! , It has occurred many times in the OrthoDB
            # This is Many to 1 Relationshit
            _cases["1:m"] += 1

            isoform = False
            og_id = odb_ncbi_ogid[ncbi_id][0]
            __desc = odb_ogid_desc[og_id][0]

            # hmmm, OrthoDB may got specific about single isoform, let's check that.
            if "isoform" in __desc:
                isoform = True
                # Last comma to prevent mismatching X1 with X10 ;)
                __desc = "isoform " + __desc.split("isoform")[-1].strip() + ","

            if isoform:
                for _id in refseq_ids:
                    line = rfsq_tr_header[_id]
                    line = line.replace("transcript variant", "isoform")
                    if __desc in line:
                        tr_id = line.split()[0].replace(">", "").strip()
                        final_map[tr_id] = og_id

            # Not specific? Ok, no problems.
            else:
                for _refseq_id in refseq_ids:
                    final_map[_refseq_id] = og_id

    # O_O Many to many relationshittt, let's map with only the gene description.
    else:
        # Many to many
        _cases["m:m"] += 1

        # get all og_ids for the ncbi_id
        og_ids = odb_ncbi_ogid[ncbi_id]
        desc_to_ogid = {}

        for _id in og_ids:
            __desc = odb_ogid_desc[_id][0]
            if "isoform" in __desc:
                # Last comma to prevent mismatching X1 with X10 ;)
                __desc = "isoform " + __desc.split("isoform")[-1].strip() + ","

            desc_to_ogid[__desc] = _id

        refseq_desc = {}

        __isoform = False
        if "isoform" in " ".join(desc_to_ogid.keys()):
            __isoform = True

        for _id in refseq_ids:
            # Will need more branching
            line = rfsq_tr_header[_id]
            if __isoform:
                line = line.replace("transcript variant", "isoform")

            tr_id = line.split()[0].replace(">", "").strip()
            #tr_desc = line.split(",")[-2].strip()

            for og_desc, og_id in desc_to_ogid.items():
                if og_desc in line:
                    final_map[tr_id] = og_id
                    break

## Check how many reads mapped
if len(final_map) != len(rfsq_tr_header):
    print("%d Transcripts matched, %d missing from total %d" %
            (len(final_map), len(rfsq_tr_header) - len(final_map), len(rfsq_tr_header)))
else:
    print("All records matched for the NCBI ID: %s" % (ncbi_id))

f = open("mapped.json", "w")
f.write(json.dumps(final_map, indent=4, sort_keys=True))
f.close()

## Print the unmatched to a separate file: unmatched_transcripts.txt

unmatched = set(rfsq_tr_header.keys()) - set(final_map.keys())

with open("unmatched_transcripts.txt", 'w') as unmatched_file:
    for _tr_id in unmatched:
        unmatched_file.write(rfsq_tr_header[_tr_id] + "\n")


print ("-----------")
print (_cases)
