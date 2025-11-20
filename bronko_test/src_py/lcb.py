def assign_buckets(kmer: int, k: int):
    
    buckets = [0]* k
    num_A = [0] * k
    val = [0] * k
    mu = [0] * k

    # Initialize
    mask = 3 << ((k - 1) * 2)  # 3lu << ((n-1)<<1)
    p = 1 << ((k - 1) * 2)     # 4^(n-1)
    cur = kmer & mask

    val[0] = kmer - cur
    mu[0] = p + ((cur >> 2) * (k - 1)) if cur else val[0]
    sum_mu = mu[0]
    # print(f"num_A: {num_A}, mask {bin(mask)}:, cur: {bin(cur)}, p:{bin(p)}, val: {val}")
    
    # Iterate through the k-mer
    for i in range(1, k):
        num_A[i] = num_A[i - 1] + (0 if cur else 1) # add to the number of As if the current base is A, else don't

        mask >>= 2               # Move mask to the right by 2 bits
        cur = kmer & mask           # Extract the current base
        p >>= 2                  # Reduce power by 4
        val[i] = val[i - 1] - cur
        mu[i] = p + ((cur >> 2) * (k - i - 1)) if cur else val[i]
        sum_mu += mu[i]
        
        # print(f"num_A: {num_A}, mask {bin(mask)}:, cur: {bin(cur)}, p:{bin(p)}, val[i]: {val}")

    # Initialize positions for storing in the buckets
    mask = 3 << ((k - 1) * 2)
    j = 0
    
    
    # Store the bucket values for that kmer
    for i in range(k):
        cur = kmer & mask
        mask >>= 2

        p = sum_mu - mu[i] + val[i] - num_A[i] * cur + 1 + num_A[i]

        buckets[j] = p
        j += 1

    return buckets

