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
deep_search = False
detailed_report = False

if len(sys.argv) > 3:
    if "--deep-search" in sys.argv:
        deep_search = True
    
    if "--detailed-report" in sys.argv:
        detailed_report = True

refseq_file = glob("species/" + species + "/refseq/*gz")[0]
orthodb_file = glob("species/" + species + "/orthodb/" + species + "*_genes.tab")[0]
gene2refseq_file = glob("species/" + species + "/*gene2refseq.txt")[0]
og2genes_file = glob("species/" + species + "/orthodb/" + species + "*_OG2genes.tab")[0]

print("Processing %s" % (species))

class Report():
    orthodb_genes_parsing = []
    analysis_summary = []
    missing_ogids = set()
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
    global refseq2ncbi
    
    mapped_ids = final_map.keys()
    names_file = "species/%s/map/%s_mapped_reads.fa.names" % (species, species)
    mapped_reads_file = "species/%s/map/%s_mapped_reads.fa" % (species, species)
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
                    REPORT.missing_ogids.add(odb_gene_id)
                    odb_group_ids = "NULL"
                    continue # Don't write it to disk.

                else:
                    odb_group_ids = ";".join(og2genes[odb_gene_id])

                description = tr_id + "|" + odb_gene_id + "|" + odb_group_ids + "|" + species
                names.write(description + "\t" + odb_gene_id + "|" + refseq2ncbi[tr_id][0] + "\n")
                fasta.description = fasta.id = fasta.name = description
                SeqIO.write(fasta, handle, "fasta")

    # print("[LOG] %d gene_ids not found in OG2Genes" % (unfound_in_og2genese))


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
    print("Deep Search is activated..")
    ODB.activate_deep_search()

print("Parsing files...")

odb_ncbi_ogid = ODB.odb_genes_info(path=orthodb_file, tax_id=taxonomy_id, key="ncbi_id", value="odb_gene_id")
odb_ogid_desc = ODB.odb_genes_info(path=orthodb_file, tax_id=taxonomy_id, key="odb_gene_id", value="description")

# REPORTING
REPORT.orthodb_genes_parsing.append(ODB._OrthoDB__report)


rfsq_tr_header = header_tr_line(refseq_file)
matched_ncbi_ids = set()


ncbi2refseq = {}
refseq2ncbi = {}

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

        #NCBI 2 REFSEQ
        if ncbi not in ncbi2refseq:
            ncbi2refseq[ncbi] = [refseq]
        else:
            ncbi2refseq[ncbi].append(refseq)
        
        #REFSEQ 2 NCBI
        if ncbi not in refseq2ncbi:
            refseq2ncbi[refseq] = [ncbi]
        else:
            refseq2ncbi[refseq].append(ncbi)

final_map = {}

_cases={"1:1":0,"1:m":0,"m:m":0}

print("Mapping...")
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
            # This is 1 to Many Relationshit
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
_total_no_matched_og_ids = len(set(final_map.values()))
_total_no_matched_transcrtipts = len(final_map)
_total_no_og_ids = len(odb_ogid_desc)
_total_no_transcripts = len(rfsq_tr_header)
_unmatched_transcripts = _total_no_transcripts - _total_no_og_ids
_matched_genes = len(matched_ncbi_ids)
_total_genes = len(odb_ncbi_ogid.keys())
_unmatched_genes = _total_genes - _matched_genes

REPORT.analysis_summary.append("%d Mapped Ortho genes from total of %d Ortho Genes" % (_total_no_matched_og_ids, _total_no_og_ids))
REPORT.analysis_summary.append("%d Mapped Refseq Transcripts from total of %d Transcripts" % (_total_no_matched_transcrtipts, _total_no_transcripts))
REPORT.analysis_summary.append("%d matched Genes from total of %d Genes" % (_matched_genes, _total_genes))



output_json = "species/%s/map/%s_" % (species, species)
f = open(output_json + "mapped.json", "w")
f.write(json.dumps(final_map, indent=4, sort_keys=True))
f.close()


## Extract matched reads
print("Writing mapped reads fasta and names files..")
write_mapped_reads(final_map, species, refseq_file)

output_report = "species/%s/map/report/%s_report_" % (species, species)

# Report Summary

with open(output_report + "summary.log", 'w') as summary_report:
    summary_report.write("~~ %s Mapping Summary ~~\n" % (species))
    summary_report.write("OrthoDB Files parsing summary\n")
    for record in REPORT.orthodb_genes_parsing[0]["summary"]:
        summary_report.write(record + "\n")
    
    summary_report.write("------------------\n")
    summary_report.write("Mapping Summary \n")

    for record in REPORT.analysis_summary:
        summary_report.write(record + "\n")


if detailed_report:
    print("Generating detailed report")
    
    # Print the unmatched to a separate file: unmatched_transcripts.txt
    unmatched = set(rfsq_tr_header.keys()) - set(final_map.keys())

    with open(output_report + "unmatched_transcripts.log", 'w') as unmatched_file:
        for _tr_id in unmatched:
            unmatched_file.write(rfsq_tr_header[_tr_id] + "\n")

    with open(output_report + "unfound_genes.log", 'w') as unfound:
        unfound.write("Failed to retreieve the NCBI IDs of these gene symbols \n\n")
        
        for _gene_symbol in REPORT.orthodb_genes_parsing[0]["unfound_geneSymbols"]:
            unfound.write(_gene_symbol + "\n")
    

    with open(output_report + "missing_ogIDs.log", 'w') as missings:
        missings.write("List of Ortho Genes not found in OG2Genes.tab file \n")
        for og in REPORT.missing_ogids:
            missings.write(og + "\n")