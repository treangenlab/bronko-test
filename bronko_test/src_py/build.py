from collections import defaultdict
import pyfastx
import pickle

from src_py.utils import nucleotide_to_canonical_binary, nucleotide_to_binary
from src_py.lcb import assign_buckets

class Bucket:
    def __init__(self, location, idx, ref_base, canonical, seq_name):
        self.location = location
        self.idx = idx
        self.ref_base = ref_base
        self.canonical = canonical
        self.seq_name = seq_name

def build(fasta, k, out='out'):
    """Builds and index from a fasta file

    Args:
        fasta (_type_): _description_
        k (_type_): _description_
    """
    
    ### General Statistics
    num_seqs = 0
    num_bases = 0
    
    ### Sequence Information
    seq_info = {}
        
    kmer_dict = defaultdict(list)

    ### Read in the fasta file using fastx
    for seq in pyfastx.Fasta(fasta):
        sequence = seq.seq
        seq_info[seq.name] = len(sequence)
        num_seqs += 1
        num_bases += len(sequence)
        
        # sequence = sequence[:20]
        for i in range(0, len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmer_bin, canon = nucleotide_to_canonical_binary(kmer)
            buckets = assign_buckets(kmer_bin, k)
            
            for j, bucket in enumerate(buckets):
                kmer_dict[bucket].append(Bucket(i, j, kmer[j], canon, seq.name))
                
    for bucket in list(kmer_dict.keys()):
        if kmer_dict[bucket][0].location == 620:
            print(bucket, kmer_dict[bucket][0].location, kmer_dict[bucket][0].idx)
        
    ## General Dictionary With Information on the Names and Lengths of the Sequences (seq_info) and the buckets (kmer_dict)                                       
    output = {
        "kmers":kmer_dict,
        "sequences":seq_info
    }
    
    with open(f'{out}/iv_kmer.pkl', 'wb') as file:
        pickle.dump(output, file)
        
    # print(num_seqs, num_bases)
    
    