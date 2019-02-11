import mygene

class OrthoDB:
    gene_name_to_ncbigid = {}
    _op_thorough = 0
    debug = 0
    
    def __init__(self, gene2refseq_file):
        self.mg = mygene.MyGeneInfo()
        with open(gene2refseq_file, 'r') as f:
            for line in f:
                line = line.strip().split(" ")
                ncbi = line[0]

                gene_name = line[2]
                if gene_name not in self.gene_name_to_ncbigid:
                    self.gene_name_to_ncbigid[gene_name] = [ncbi]
                else:
                    self.gene_name_to_ncbigid[gene_name].append(ncbi)

    def activate_deep_search(self):
        self._op_thorough = True

    def get_ncbi_id(self, gene_symbol, tax_id):
        result = self.mg.query(gene_symbol, species=tax_id, fields="_id")["hits"]
        if len(result):
            result = result[0]["_id"]
            if self.debug:
                print("GENE: %s NCBI_ID: %s" % (gene_symbol, result[0]["_id"]))

            return result
        return False
    
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
            "description": -1
        }

        if key not in info_location.keys() or value not in info_location.keys():
            print ("Please select from [uq_ortho_gene_id, tax_id, prot_seq_id, uniprot_id, ensemble_gene_name, ncbi_gid, description]")
            return 0

        result = {}
        key_idx = info_location[key]
        val_idx = info_location[value]
        total = 0
        found = 0
        __report = {"unfound_geneSymbols":set()}

        with open(path, "r") as f:
            for line in f:
                total += 1

                if len(line) < 5:
                    continue

                line = line.replace("\n", "").split("\t")

                if str(tax_id) not in line[1]:
                    continue

                # "Check if the key is NCBI ID, and because it's not available, attempt to get ncbi_id from {gene_name_to_ncbigid}"
                if key_idx == 5:
                    _gene_name = line[3]
                    
                    if _gene_name not in self.gene_name_to_ncbigid:
                        # if Option: Deep search is activated try to query it online
                        if self._op_thorough:
                            _received_ncbi_id = self.get_ncbi_id(_gene_name, tax_id)

                            if _received_ncbi_id: # Check it's not empty
                                self.gene_name_to_ncbigid[_gene_name] = _received_ncbi_id
                            else:
                                # Still can't find it? reporting the unfound geneSymbols
                                __report["unfound_geneSymbols"].add(_gene_name)
                                continue

                        else:
                            # Option: Deep search not activated
                            __report["unfound_geneSymbols"].add(_gene_name) # reporting the unfound geneSymbol
                            continue 
                       
                    # Gene Symbol found
                    feature_key = self.gene_name_to_ncbigid[_gene_name][0]
                
                else:
                    feature_key = line[key_idx]

            # Check if the value is NCBI ID, and because it's not available, attempt to get ncbi_id from {gene_name_to_ncbigid}
                if val_idx == 5:
                    _gene_name = line[3]

                    if _gene_name not in self.gene_name_to_ncbigid:
                        # if Option:Thorough is activated try to query it online
                        if self._op_thorough:
                            _received_ncbi_id = self.get_ncbi_id(_gene_name, tax_id)

                            if _received_ncbi_id: # Check it's not empty
                                self.gene_name_to_ncbigid[_gene_name] = _received_ncbi_id
                            else:
                                # Still can't find it? reporting the unfound geneSymbol
                                __report["unfound_geneSymbols"].add(_gene_name)
                                continue

                        else:
                            # Option:Thorough not activated
                            __report["unfound_geneSymbols"].add(_gene_name) # reporting the unfound geneSymbol
                            continue 
                       
                    # Gene Symbol found
                    feature_value = self.gene_name_to_ncbigid[_gene_name][0]

                else:
                    feature_value = line[val_idx]
                

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
        if self.debug:
            print ("(Tax_ID %d) %d %s found (%d uniq) in %d record" % (tax_id, found, value, len(set(all_values)), total))
        
        return result
