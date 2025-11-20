from collections import defaultdict
import pyfastx
import pickle
import numpy as np
import subprocess
import time
import os

from src_py.utils import nucleotide_to_canonical_binary, binary_to_nucleotide
from src_py.lcb import assign_buckets

def screen(fastq:str, 
           index:str, 
           k:int, 
           out:str,
           threads:int=3,
           min_kmer_count:int=10, 
           strand_filter=True,
           n_per_strand=2,
           kmer_edges=True,
           edge_filter_len=2,
           debug=False):
    """ Screens a fastq file against an index by counting the kmers and applying a (1-2) bucketing function to rapidly map
    kmers and their counts to local regions

    Args:
        fastq (str): input path of fastq file
        index (str): input path of the bronko index to be used
        k (int): kmer size
        out (str): output location
        min_kmer_count (int, optional): Minimum number of kmers to set for kmer counting. Defaults to 10.
        strand_filter (bool, optional): Filter minor variation by strand frequencies. Defaults to True.
        n_per_strand (int, optional): Min number of kmers that much support each strand at a particular location to be called. Defaults to 2.
        debug (bool, optional): Print extra debugging information. Defaults to False.
    """

    start_time = time.time()
    
    if edge_filter_len > 3: 
        raise ValueError("edge filter too large")
    
    with open(index, 'rb') as file:
        input = pickle.load(file)
        kmer_index = input["kmers"]
        seq_info = input["sequences"]
            
    ### General Statistics
    num_reads = 0
    num_bases = 0
    total_kmers = 0
    
    ## kmer counts
    kmers = defaultdict(int) ## TODO: don't save info, just read directly and compute immediately
           
    # if os.path.exists(f'{out}/mer_counts.fa'): 
    #     subprocess.run([
    #         'rm', f'{out}/mer_counts*'
    #     ], check=True)

    
    ### Read in the fasta file using jellyfish
    # subprocess.run([
    #     'jellyfish', 'count',
    #     '-m', str(k),
    #     '-s', '100M',
    #     '-t', str(10),
    #     '-L', str(min_kmer_count),
    #     fastq], check=True)
    
    if not os.path.exists("/home/Users/rdd4/tmp"):
        subprocess.run([
            "mkdir", "/home/Users/rdd4/tmp"], check=True)
    
    subprocess.run([
        'kmc',
        f'-k{k}'
        '-m2'
        f'-t{threads}',
        '-b',
        f'-ci{min_kmer_count}',
        '-cs100000',
        fastq,
        f"{out}/mer_counts.res",
        f"/home/Users/rdd4/tmp"], check=True)
    
    print(f"Kmers Counted in {time.time()-start_time}s")
    
    # subprocess.run([
    #     'jellyfish', 'dump',
    #     'mer_counts.jf',
    #     '-o', f'{out}/mer_counts.fa'], check=True)
    
    subprocess.run([
        'kmc_tools', 'transform',
        f"{out}/mer_counts.res",
        "dump", f"{out}/mer_counts.txt"
        ], check=True)
    
    print(f"Kmers Counted + Dumped in {time.time()-start_time}s")
    
    with open(f"{out}/mer_counts.txt", 'r') as f:
        for line in f:
            data = line.strip().split()
            kmer, count = data[0], int(data[1]) 
            kmer_bin, canon = nucleotide_to_canonical_binary(kmer)
            kmers[kmer] = {"canonical_kmer":kmer_bin, "n": count, "rc":canon}

    ## reading from jellyfish
    # for seq in pyfastx.Fasta(f"{out}/mer_counts.fa"):

    #     kmer = seq.seq.upper() # TODO: make formal fix with uncertain nucleotides (N), etc
    #     n = int(seq.name)
        
    #     kmer_bin, canon = nucleotide_to_canonical_binary(kmer)
    #     kmers[kmer] = {"canonical_kmer":kmer_bin, "n": int(seq.name), "rc":canon}
    
    print(f"Kmers Counted + Dumped + Read in {time.time()-start_time}s")
    print(f'Unique kmers above count {min_kmer_count}: {len(kmers)}')
    
    # output = [[0]*4 for _ in range(30000)]
    output = {seq: {"counts":np.zeros((seq_length, 4), dtype=int), "ref":np.full((seq_length), '', dtype='str')} for seq, seq_length in seq_info.items()}
    output_rev = {seq: {"counts":np.zeros((seq_length, 4), dtype=int), "ref":np.full((seq_length), '', dtype='str')} for seq, seq_length in seq_info.items()} ##TODO: revisit whether there is a better way to store the reverse information
    output_counts = {seq: {"counts":np.zeros((seq_length, 4), dtype=int), "ref":np.full((seq_length), '', dtype='str')} for seq, seq_length in seq_info.items()}
    output_rev_counts = {seq: {"counts":np.zeros((seq_length, 4), dtype=int), "ref":np.full((seq_length), '', dtype='str')} for seq, seq_length in seq_info.items()} ##TODO: revisit whether there is a better way to store the reverse information

    unmapped_kmers = 0
    uniquely_mapped_kmers = 0
    doubly_mapped_kmers = 0
    perfectly_mapped_kmers = 0
    rc_kmers = 0
        
        
    seq_to_check = 'PV249929.1'
    pos_to_check = 1394
        
    for kmer, kmer_info in list(kmers.items()):
        try:
            n = kmer_info["n"]
        except:
            print(kmer, kmer_info)
        rc = kmer_info["rc"]
        kmer_bin = kmer_info["canonical_kmer"]
        
        
        buckets = assign_buckets(kmer_bin, k)
        buckets_matched = 0   
        
        if kmer_edges:
            buckets = buckets[edge_filter_len:-edge_filter_len-1]
        
        for bucket in buckets:
            if len(kmer_index[bucket])==1: ##only unique hits (not any buckets that are mapped to multiple locations), should play with this though
                bucket_info = kmer_index[bucket][0] ##get first bucket 
                genome_pos = bucket_info.location
                sequence = bucket_info.seq_name
                nuc_x = bucket_info.idx
                ref = bucket_info.ref_base
                ref_canon = bucket_info.canonical
                
                
                if ref_canon:
                    pos = k - nuc_x - 1
                    extracted_bits = (kmer_bin >> (2 * (k - pos - 1))) & 0b11
                    extracted_bits ^= 0b11
                    if rc:
                        output_counts[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] += 1
                        if output[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] < n: ## take the maximum of any hits
                            output[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] = n
                            output[sequence]["ref"][genome_pos + nuc_x] = ref
                    else:
                        output_rev_counts[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] += 1
                        if output_rev[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] < n: ## take the maximum of any hits
                            output_rev[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] = n
                            output_rev[sequence]["ref"][genome_pos + nuc_x] = ref
                    
                else:
                    pos = nuc_x
                    extracted_bits = (kmer_bin >> (2 * (k - pos - 1))) & 0b11
                    if rc:
                        output_rev_counts[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] += 1
                        if output_rev[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] < n: ## take the maximum of any hits
                            output_rev[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] = n
                            output_rev[sequence]["ref"][genome_pos + nuc_x] = ref
                    else:
                        output_counts[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] += 1
                        if output[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] < n: ## take the maximum of any hit
                            output[sequence]["counts"][genome_pos + nuc_x, int(extracted_bits)] = n
                            output[sequence]["ref"][genome_pos + nuc_x] = ref

                if debug:     
                    if kmer == "GCTGCATGTTTGTATAACG":
                        print("forward",binary_to_nucleotide(kmer_bin, k), n, genome_pos, nuc_x, bin(extracted_bits), rc, bucket)
                    if kmer == "CGTTATACAACCATGCAGC":
                        print("reverse", binary_to_nucleotide(kmer_bin, k), n, genome_pos, nuc_x, bin(extracted_bits), rc, bucket)
                    
                    if sequence == seq_to_check and genome_pos + nuc_x + 1 == pos_to_check:
                        print(kmer, binary_to_nucleotide(kmer_bin, k), n, genome_pos, nuc_x, pos, binary_to_nucleotide(extracted_bits, 1), rc, ref_canon, bucket)

                buckets_matched += 1

            elif len(kmer_index[bucket]) > 1 and debug:
                print(f"Found multiple buckets for bucket: {bucket}")
                for bucket_info in kmer_index[bucket]:
                    print(bucket_info.location, bucket_info.idx)

            
        if buckets_matched == 0:
            unmapped_kmers += 1
        elif buckets_matched == 1:
            uniquely_mapped_kmers += 1
        elif buckets_matched < k:
            doubly_mapped_kmers += 1
        elif buckets_matched == k:
            perfectly_mapped_kmers += 1
            
    # if sorted(output.keys()) != sorted(output_rev.keys()):
    #     raise ValueError("Reverse strands do not cover same segments as output")
                
    print(f'\nUnmapped Kmers: {unmapped_kmers}/{len(kmers)}={unmapped_kmers/len(kmers)*100}%')     
    print(f'Uniquely Mapped Kmers: {uniquely_mapped_kmers}/{len(kmers)}={uniquely_mapped_kmers/len(kmers)*100}%')     
    print(f'Doubly Mapped Kmers: {doubly_mapped_kmers}/{len(kmers)}={doubly_mapped_kmers/len(kmers)*100}%')     
    print(f'Perfect Mapped Kmers: {perfectly_mapped_kmers}/{len(kmers)}={perfectly_mapped_kmers/len(kmers)*100}%')     
    print(f'RC Kmers: {rc_kmers}/{len(kmers)}={rc_kmers/len(kmers)*100}%')     
                    
    print("\nWriting Output") 
    ## bronko_variants.tsv --> tsv file with counts at variant positions
    ## bronko.vcf --> vcf format of the same variants
    with open(f'{out}/bronko.tsv', "w") as f, open(f'{out}/bronko.vcf', "w") as f_vcf:
        
        # write tsv header
        f.write("reference\tindex\tref\tA\tC\tG\tT\ta\tc\tg\tt\n")
        
        # Write VCF header
        f_vcf.write("##fileformat=VCFv4.2\n")
        for seq_name, length in seq_info.items():
            f_vcf.write(f"##contig=<ID={seq_name},length={length}>\n")
        f_vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f_vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        f_vcf.write("##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Strand Bias P-value\">\n")
        f_vcf.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"FwdRef,RevRef,FwdAlt,RevAlt\">\n")
        f_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        bases = ['A', 'C', 'G', 'T']
        base_to_index = {b: i for i, b in enumerate(bases)}
        
        for (seq, seq_data), (_, seq_data_rev), (_, seq_counts), (_, seq_counts_rev) in zip(output.items(), output_rev.items(), output_counts.items(), output_rev_counts.items()): ##loop over every sequence
            
            map_info = seq_data["counts"]
            map_info_rev = seq_data_rev["counts"]
            ref_forward = seq_data["ref"]
            ref_reverse = seq_data_rev["ref"]
            count_info = seq_counts["counts"]
            count_info_rev = seq_counts_rev["counts"]

            for i, (row, row_rev, ref, ref_rev, count, count_rev) in enumerate(zip(map_info, map_info_rev, ref_forward, ref_reverse, count_info, count_info_rev)): ##loop through each index and add that line
                ref = ref if ref != '' else ref_rev
                pos = i + 1
                
                ## Generate count table
                # row_total = [(row[base] if count[base] > 2 else 0) + (row_rev[base] if count_rev[base] > 2 else 0) for base in range(len(row))]
                if strand_filter:
                    row_total = [(row[base]+row_rev[base] if (count[base] >= n_per_strand and count_rev[base] >= n_per_strand) else 0) for base in range(len(row))]
                else:
                    row_total = [(row[base]+row_rev[base]) for base in range(len(row))]
                    
                f.write(f"{seq}\t{pos}\t{ref}\t" + "\t".join(map(str, row)) + "\t" + "\t".join(map(str, row_rev)) + "\n")

                
                base_counts = dict(zip(bases, row_total))
                                
                if debug and seq == seq_to_check and i + 1 == pos_to_check:
                    print(row, row_rev)     
                    
                if pos <= k or pos >= len(map_info) - k + 1: ##filter out ends of sequences
                    continue
                
                total_depth = sum(base_counts.values())
                if total_depth == 0:
                    continue
                
                ref_count = base_counts[ref]
                alt_candidates = [b for b in bases if b != ref and base_counts.get(b, 0) > 0]

                for alt in alt_candidates:
                    alt_count = base_counts[alt]
                    af = alt_count / total_depth
                    if af < 0.01:
                        continue

                    # DP4: fwd_ref, rev_ref, fwd_alt, rev_alt
                    ref_idx = base_to_index[ref]
                    alt_idx = base_to_index[alt]
                    dp4 = [
                        row[ref_idx],
                        row_rev[ref_idx],
                        row[alt_idx],
                        row_rev[alt_idx]
                    ]
                    info = f"DP={total_depth};AF={af:.5f};SB=0;DP4={','.join(map(str, dp4))}"
                    f_vcf.write(f"{seq}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
                    

            
            
    ## Print out general counts table  
    ## bronko.tsv --> all kmers   
    with open(f'{out}/bronko_counts.tsv', "w") as f:
        f.write("reference\tindex\tref\tA\tC\tG\tT\ta\tc\tg\tt\n")
        for (seq, seq_data), (_, seq_data_rev) in zip(output_counts.items(), output_rev_counts.items()): ##loop over every sequence
            map_info = seq_data["counts"]
            map_info_rev = seq_data_rev["counts"]
            ref_forward = seq_data["ref"]
            ref_reverse = seq_data_rev["ref"]
            for i, (row, row_rev, ref, ref_rev) in enumerate(zip(map_info, map_info_rev, ref_forward, ref_reverse)): ##loop through each index and add that line
                ref = ref if ref != '' else ref_rev
                f.write(f"{seq}\t{i+1}\t{ref}\t" + "\t".join(map(str, row)) + "\t" + "\t".join(map(str, row_rev)) + "\n")
                if debug and seq == seq_to_check and i + 1 == pos_to_check:
                    print(row, row_rev)

    