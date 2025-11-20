def nucleotide_to_binary(kmer: str) -> int:
    """Convert a DNA k-mer string to a binary representation using 2 bits per base."""
    base_to_bits = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
    
    binary = 0
    for base in kmer:
        binary = (binary << 2) | base_to_bits[base]
    
    return binary

def binary_to_nucleotide(binary: int, k: int) -> str:
    """Convert a binary representation to a DNA k-mer string."""
    bits_to_base = ['A', 'C', 'G', 'T']
    
    kmer = []
    for _ in range(k):
        # Extract the last 2 bits
        base = binary & 0b11  
        kmer.append(bits_to_base[base])
        binary >>= 2

    # Reverse since the bases are extracted in reverse order
    return ''.join(reversed(kmer))

def reverse_complement_64(kmer_bin: int, kmer_size: int) -> int:
    """Compute the reverse complement of a 2-bit encoded k-mer packed in an int.

    Assumes kmer_bin uses 2*kmer_size bits. Mimics the bitwise operations in the C++ ReverseComp64 function.
    """
    # Step 1: Complement the k-mer (note: only flip the lower 2*k bits)
    mask = (1 << (2 * kmer_size)) - 1
    res = ~kmer_bin & mask

    # Step 2: Reverse the bits in 2-bit chunks (bit twiddling from C++)
    res = ((res >> 2) & 0x3333333333333333) | ((res & 0x3333333333333333) << 2)
    res = ((res >> 4) & 0x0F0F0F0F0F0F0F0F) | ((res & 0x0F0F0F0F0F0F0F0F) << 4)
    res = ((res >> 8) & 0x00FF00FF00FF00FF) | ((res & 0x00FF00FF00FF00FF) << 8)
    res = ((res >> 16) & 0x0000FFFF0000FFFF) | ((res & 0x0000FFFF0000FFFF) << 16)
    res = ((res >> 32) & 0x00000000FFFFFFFF) | ((res & 0x00000000FFFFFFFF) << 32)

    # Step 3: Shift to align the result to the rightmost bits
    res >>= 2 * (32 - kmer_size)

    return res


def nucleotide_to_binary(kmer: str) -> int:
    base_to_bits = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11, 'N':0b00} ## TODO: fix Ns
    binary = 0
    for base in kmer:
        binary = (binary << 2) | base_to_bits[base]
    return binary

def nucleotide_to_canonical_binary(kmer: str) -> tuple[int, bool]:
    kmer_bin = nucleotide_to_binary(kmer)
    revcomp_bin = reverse_complement_64(kmer_bin, len(kmer))
    return (min(kmer_bin, revcomp_bin), revcomp_bin < kmer_bin)

def binary_rc(kmer: int, k:int):
    """Gives reverse complement of kmer"""
    bitmask = (1 << k*2) - 1
    return ~kmer & bitmask

# def nucleotide_to_canonical_binary(seq: str) -> bytes:
#     mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
#     binary = ''.join(mapping[nt] for nt in seq.upper())
#     if binary[0] == '1':
#         binary = ''.join('1' if b == '0' else '0' for b in binary)
#     return int(binary, 2)