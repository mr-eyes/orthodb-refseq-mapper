import gzip
from Bio import SeqIO
from tqdm import tqdm

class RefSeq:

    def __init__(self, file_name):
        self.file_name = file_name

    def get_transcript_ids(self):
        tr_ids = []
        with open(self.file_name,'r') as f:
            for line in f:
                if line[0] != ">":
                    continue

                tr_ids.append(line.split()[0].replace(">", "").strip())
        
        return tr_ids
                
                

