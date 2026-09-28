import argparse
import os
import time

from src_py import io
from src_py.build import build
from src_py.screen import screen
from src_py.bench import bench

import logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s', 
    datefmt='%Y-%m-%d %H:%M:%S'
)

DEFAULT_KMER_SIZE = 15

def main():
    parser = argparse.ArgumentParser(description="Rapid Intrahost Variation Calling")
    subparsers = parser.add_subparsers(dest='command')
    
    build_command = subparsers.add_parser('build', help='Build an index from a reference')
    build_command.add_argument('-i', type=str, help='Path to fasta file')
    build_command.add_argument('-o', type=str, default='out', help='Path to output folder')
    build_command.add_argument('-k', default=DEFAULT_KMER_SIZE, type=int, help=f'K-mer size (default={DEFAULT_KMER_SIZE})')
    
    screen_command = subparsers.add_parser('screen', help='Run sequencing data against an index')
    screen_command.add_argument('-i', type=str, help='Path to fastq file')
    screen_command.add_argument('-f', type=str, help='Path to index')
    screen_command.add_argument('-o', type=str, default='out', help='Path to output')
    screen_command.add_argument('-k', default=DEFAULT_KMER_SIZE, type=int, help=f'K-mer size (default={DEFAULT_KMER_SIZE})')
    
    bench_command = subparsers.add_parser('bench', help='Run build+screen, and compare to bowtie2')
    bench_command.add_argument('-fa', type=str, help='Fasta reference')
    bench_command.add_argument('-fq', nargs="*", help='Single end fastq data')
    bench_command.add_argument('--r1', nargs="*", help='Paired end R1 fastq data')
    bench_command.add_argument('--r2', nargs="*", help='Paired end R2 fastq data')
    bench_command.add_argument('-o', type=str, default='out', help='Path to output')
    bench_command.add_argument('--min-af', type=float, default=0.03, help='Min allele frequency to bench')
    bench_command.add_argument('-k', default=DEFAULT_KMER_SIZE, type=int, help=f'K-mer size (default={DEFAULT_KMER_SIZE})')
    bench_command.add_argument('--no-rerun', default=True, action='store_false', help=f'Only re-run bronko (do not re-run ivar+lofreq). Warning will delete timing info')
    bench_command.add_argument('--ivar-all', dest='ivar_require_pass', default=True, action='store_false', help='Compare against every iVar row, including ones iVar flagged PASS=FALSE (default: PASS=TRUE only)')
    bench_command.add_argument('-t', '--threads', default=10, help=f'Number of threads')

    args = parser.parse_args()
    k = args.k
    out = args.o
    
    
    start = time.time()
    
    if args.command == 'build':
        seq_data = args.i

        print(k, seq_data)
        build(seq_data, k, out)
        
        print(f'Build done in {time.time()-start}')

    if args.command == 'screen':
        seq_data = args.i

        index = args.f
        screen(seq_data, index, k, out)
        print(f'Screen done in {time.time()-start}s')
        
    if args.command == 'bench':
        if not args.fq or not (args.r1 and args.r2):
            if args.fq or (args.r1 and args.r2 and len(args.r1) == len(args.r2)):
                bench(args.fq, args.r1, args.r2, args.fa, out, k=k, min_af=args.min_af, re_run=args.no_rerun, threads=args.threads, ivar_require_pass=args.ivar_require_pass)
            else:
                print(args.r1, args.r2)
                print("Number of r1 and r2 do not match")
        else:
            print("No reads provided")
        


if __name__=="__main__":
    main()