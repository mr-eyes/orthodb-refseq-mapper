class OrthoDB:
    def __init__(self):
        pass

    def odb_genes_info(self, path, tax_id, key, value):
        """
        File: odb*genes.tab
        info = [og_id, tax_id, prot_seq_id, uniprot_id, ensemble_gene_name, ncbi_gid, description]
        return : dict {uq_ortho_gene_id : info}
        """

        info_location = {
            "og_id": 0,
            "tax_id": 1,
            "prot_seq_id": 2,
            "uniprot_id": 3,
            "ensemble_gene_name": 4,
            "ncbi_gid": 5,
            "description": 6
        }

        if key not in info_location.keys() or value not in info_location.keys():
            print (
                "Please select from [uq_ortho_gene_id, tax_id, prot_seq_id, uniprot_id, ensemble_gene_name, ncbi_gid, description]")
            return 0

        result = {}
        key_idx = info_location[key]
        val_idx = info_location[value]
        total = 0
        found = 0

        with open(path, "r") as f:
            next(f)  # skip header
            for line in f:
                line = line.replace("\n", "").split("\t")
                if str(tax_id) not in line[1]:
                    continue

                feature_key = line[key_idx]
                feature_value = line[val_idx]
                total += 1

                if len(feature_value) > 2:
                    _value_ = feature_value
                    found += 1
                else:
                    _value_ = "-"

                if feature_key not in result:
                    result[feature_key] = [_value_]
                else:
                    result[feature_key].append(_value_)

        all_values = [x for v in result.values() for x in v]
        print ("(Tax_ID %d) %d %s found (%d uniq) in %d record" %
               (tax_id, found, value, len(set(all_values)), total))
        return result
