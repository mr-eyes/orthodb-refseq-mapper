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
                print("GENE: %s NCBI_ID: %s" % (gene_symbol, result))

            return result
        return False

    def deep_search(self, gene_name):
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
            "description": -1
        }

        if key not in info_location.keys() or value not in info_location.keys():
            print ("Please select from [uq_ortho_gene_id, tax_id, prot_seq_id, uniprot_id, ensemble_gene_name, ncbi_gid, description]")
            return 0

        result = {}
        key_idx = info_location[key]
        val_idx = info_location[value]
        total = 0
        values_found = 0
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
                    found = False

                    # It may have ";" (multiple synonymus) check and iterate then try to find
                    if ";" in _gene_name:
                            multiple_synonymus = 1
                            _gene_names_list = set(_gene_name.split(";"))
                            _intersection = list(_gene_names_list.intersection(set(self.gene_name_to_ncbigid.keys())))
                            
                            if len(_intersection) == 1:
                                _gene_name = _intersection[0] # The first and only gene and it's found
                                found = True
                                
                    else:
                        multiple_synonymus = 0
                    
                    
                    # If the gene name is not found in the gene2refseq file.
                    if _gene_name not in self.gene_name_to_ncbigid:
                        if not multiple_synonymus:  # Single gene_name
                                _gene_names_list = [_gene_name]
                        
                        
                        # Deep search is activated try to query it online
                        if self._op_thorough and not found:
                            # if multiple (contains ;), _gene_names_list already set before.

                            for _gene in _gene_names_list:
                                    _received_ncbi_id = self.get_ncbi_id(_gene, tax_id)
                                    if _received_ncbi_id:  # Check it's not empty
                                        self.gene_name_to_ncbigid[_gene] = _received_ncbi_id # Add to the dict for future use
                                        found = True
                                        break # Just break, we've found it's NCBI ID

                            # Still not found? report to the missings. Skip the current iteration
                            if not found:                            
                                __report["unfound_geneSymbols"].add(_gene_name)
                                continue


                        # Option: Deep search not activated
                        # report to the missings. Skip the current iteration
                        else:
                            __report["unfound_geneSymbols"].add(_gene_name) # reporting the unfound geneSymbol
                            continue
                       
                    # Iteration reached here? then the Gene Symbol has been found
                    feature_key = self.gene_name_to_ncbigid[_gene_name][0]
                
                # The requested feature is not NCBI ID.
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
                            # Option:Thorough not activated, so just report the missing.
                            __report["unfound_geneSymbols"].add(_gene_name) # reporting the unfound geneSymbol
                            continue 
                       
                    # Gene Symbol found
                    feature_value = self.gene_name_to_ncbigid[_gene_name][0]

                else:
                    feature_value = line[val_idx]
                

                if len(feature_value) > 2:
                    _value_ = feature_value
                    values_found += 1
                    
                else:
                    # No real value
                    continue

                if feature_key not in result:
                    result[feature_key] = [_value_]
                else:
                    result[feature_key].append(_value_)

        all_values = [x for v in result.values() for x in v]
        if self.debug:
            #_unit_report =  "Creation of dict(%s:%s) from a total of %d records\n" % (key, value, total)
            _unit_report = "%d %ss mapped to %d %ss, distinct values count: %d\n" % (len(result), key.upper() , len(all_values), value.upper(), len(set(all_values)))
            print(_unit_report)
            
        return result
