import numpy as np
import random
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio import SeqIO
import seaborn as sns

NUCLEOTIDES  = {'A', 'C', 'G', 'T'}

def write_random_genome(n:int, nucleotides=NUCLEOTIDES):
    """If want to simply write a completely random genome, instead of inputting a fasta file
    
    Writes a genome of length n (segmented genomes unavailable for now)

    Args:
        n (int): length of genome
        nucleotides (set, optional): set of characters to pull from for sequence. Defaults to NUCLEOTIDES.
        
    Returns:
        string: sequence of length n containing nucleotides
    """
    
    return ''.join(random.choices(list(nucleotides), k=n))

def input_genome(fasta:str):
    """Reads in a fasta file and returns the genome sequence (if segmented, will all be appended together)

    Args:
        fasta (str): location of fasta file
        
    Returns:
        string: sequence of genome within fasta file
    """
    seq = ""
    for record in SeqIO.parse(fasta, "fasta"):
        seq = seq + record.seq
    
    return str(seq)

def mutate_genome(genome:str, mutation_rate:float=0.01, generations:int=1, rates:list=[0.9, 0.05, 0.05]):
    """Returns evolved version of genome at given mutation rate over generations
    
    Args:
        genome(str): base genome
        mutation_rate(float): what percentage of genome will change after 1 generation (float between 0.0-1.0)
        generations (int): number of generations
        rates (list[tuple, tuple, tuple]): percentage of mutations that are mutations, insertions, and deletions
    """
    num_nt = len(genome)
    
    for _ in range(generations):
        for _ in range(int(mutation_rate*num_nt)):
            
            loc = np.random.randint(0, num_nt)
            mutation_type = np.random.choice(['MUTATION', 'INSERTION', 'DELETION'], p=rates)
            if mutation_type == 'MUTATION':
                genome = genome[:loc] + random.choice(list(NUCLEOTIDES)) + genome[loc+1:]
            elif mutation_type == 'INSERTION':
                insertion_length = int(1+np.random.exponential(0.5))
                genome = genome[:loc] + ''.join(random.choices(list(NUCLEOTIDES), k=insertion_length)) + genome[loc+1:]
            elif mutation_type == 'DELETION':
                deletion_length = int(1+np.random.exponential(0.5))
                genome = genome[:loc-int(deletion_length/2)] + genome[loc+1+int(deletion_length/2):]
    return genome

def get_strains(genome:str, n_strains:int=1, min_seq_identity:float=0.99, rates:list=[1.0, 0.0, 0.0]):
    """Will get multiple strains (of varying nt identity) to the seed genome (only SNP variation)
    
    genome: seed genome
    n_strains: number of strains that you want to generate
    min_seq_identity: how far the genomes can be from the seed (0.8 == 80% nt same) 
    
    """
    mutation_rates = np.random.rand(n_strains) * (1-min_seq_identity)
    genomes = [mutate_genome(genome, mutation_rate=mutation_rate, rates=rates) for mutation_rate in mutation_rates]
    return genomes


def sequence_genome(genome:str, read_length:int=150, depth:int=30, sequencing_error_rate=0.01, plot=False):
    """Returns read set after 'sequencing' genome, to certain average depth with sequencing error rate"""
    genome_length = len(genome)
    num_reads = int(genome_length / read_length * depth)
    seq_locs = np.random.randint(0, genome_length-read_length, size=num_reads)
    
    reads = []
    for i in seq_locs:
        read = genome[i:i+read_length]
        read_w_error = mutate_genome(read, sequencing_error_rate, rates=[1.0, 0.0, 0.0])
        reads.append(read_w_error)
        
    if plot:
        depths = [0]*genome_length
        for seq_loc in seq_locs:
            depths[seq_loc: seq_loc+read_length] = [x+1 for x in depths[seq_loc: seq_loc+read_length]]
        
        plt.bar(range(genome_length), depths)
        
    return reads


def metagenomic_sequencing(genome:str, n_strains:int=10, seq_similarity:float=0.99, read_length:int=150, min_depth:float=1.0, max_depth:float=100.0, sequencing_error_rate:int=0.01, strain_depths:list=[], output=False):
    """Performs 'metagenomic' sequencing on set of n_strains
    
    Args:
        genome (str): base genome to mutate into strains and sequence
        n_strains (int): the number of distinct lineages (for intrahost variation this is probably low, with high seq similarity)
        seq_similarity (float): how similar the strains will be to the base genome
        read_length (int): length of all reads (will be updated to have some variation)
        min_depth (float): minimum sequencing depth for any strain (>0)
        max_depth (float): depth of the original sequence (no other strains will be more abundant)
        sequencing_error_rate (flaot): error rate induced into reads
    """
    strains = get_strains(genome, n_strains, seq_similarity) + [genome]    
    
    if len(strain_depths) == n_strains:
        depths = strain_depths
    else:
        depths = np.random.lognormal(1, 1, size=n_strains)
        depths = np.append([max(int(min_depth), int(depth)/100 * max_depth) for depth in depths], max_depth)
    print(depths)
    
    reads = []
    for i, (strain, depth) in enumerate(zip(strains, depths)):
        reads += sequence_genome(strain, read_length=read_length, depth=depth, sequencing_error_rate=sequencing_error_rate)
        if output:
            print(f'Sequenced strain {i} at depth {depth}')
    return reads, strains, depths


def write_fasta(sequences, filename, ids=[]):
    """Writes a list of DNA sequences to a FASTA file.

    Args:
        sequences: A list of DNA sequence strings.
        ids: A list of sequence identifiers (one for each sequence).
        filename: The name of the output FASTA file.
    """
    if len(ids)==0:
        ids = range(len(sequences))
        
    with open(filename, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(">" + str(ids[i]) + "\n")
            f.write(seq + "\n")
            
def write_fastq(sequences, filename, ids=[]):
    """Writes a list of DNA sequences to a FASTQ file with best quality scores.

    Args:
        sequences: A list of DNA sequence strings.
        ids: A list of sequence identifiers (one for each sequence).
        filename: The name of the output FASTQ file.
    """
    if len(ids) == 0:
        ids = range(len(sequences))
    
    with open(filename, 'w') as f:
        for i, seq in enumerate(sequences):
            # Write the FASTQ entry
            f.write("@" + str(ids[i]) + "\n")  # Sequence identifier
            f.write(seq + "\n")                # Sequence
            f.write("+\n")                     # Separator
            f.write('I' * len(seq) + "\n")     # Quality score (all 'I' for best quality)
