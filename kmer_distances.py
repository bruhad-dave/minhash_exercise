## importing
from importlib.resources import path
from itertools import combinations
import os
import argparse
import mmh3 ## python wrapper of MurmurHash3 (https://pypi.org/project/mmh3/); tried this one as it was mentinoed in the instructions + stuck with this one since it's much faster.
import pandas as pd
import Bio.SeqIO as sio
#from simhash import Simhash ## simhash had the intuitive feature of generating similar hashes for similar inputs, but proved quie slow hashing large lists.

## argparse
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minhash", action="store_true", help="First of two fasta files to process.")
parser.add_argument("-d", "--fasta_dir", type=str, help="A directory containing fasta files to process.")
parser.add_argument("-f1", "--fasta_1", type=str, help="First of two fasta files to process. Incompatible with -d/--fasta_dir.")
parser.add_argument("-f2", "--fasta_2", type=str, help="Second of two fasta files to process. Incompatible with -d/--fasta_dir.")
parser.add_argument("-k", "--kmer_size", type=int, help="Length of kmer to generate. Requires -m/--minhash.")
parser.add_argument("-s", "--sketch_size", type=int, nargs='?', const=1, default=1000, help="Number of hashes to include in sketch. Defaults to 1000. Requires -m/--minhash.")
parser.add_argument("-c", "--cyclic", type=bool, nargs='?', const=1, default=False, help="Is the input cyclic? Defaults to False. If file contains multiple fasta records, this will be set to False for that file.")
args = parser.parse_args()

## get args
minhash = args.minhash
fasta_dir = args.fasta_dir
fasta_1 = args.fasta_1
fasta_2 = args.fasta_2
kmer_size = args.kmer_size
sketch_size = args.sketch_size
cyclic = args.cyclic

## raise exception if fasta_dir provided with fasta_1 or fasta_2
if fasta_dir and (fasta_1 or fasta_2):
    raise Exception("\n\n>>> Please provide EITHER a directory containing fasta files OR exactly two fasta files!\n\n")
## raise exception if only one fasta file provided
if (fasta_1 and not fasta_2) or (fasta_2 and not fasta_1):
    raise Exception("\n\n>>> Please provide exactly two fasta files!\n\n")

## checkpoint -- checking valid arguments
print(f"Do hashing? {minhash}, \nFasta Dir: {fasta_dir}, \nFile 1: {fasta_1}, \nFile 2: {fasta_2}, \nK-mer length: {kmer_size}, \nSketch size: {sketch_size}, \nCyclic? {cyclic}")

## generate kmers from sequence
def get_kmers(seq, k, cyclic:False):
    """Generate kmers of length k from a given string.
    Function adapted from a de Bruijn graph tutorial via the Eaton Lab (https://eaton-lab.org/slides/genomics/answers/nb-10.2-de-Bruijn.html)

    Args:
        seq (string): The string from which to generate kmers
        k (int): length of kmer
        cyclic (Bool): Is the given string cyclic/circular?

    Returns:
        Dict: Returns a dictionary containing unique kmers as keys and corresponding counts as values.
    """
    kmers = {}
    for i in range (0, len(seq)):
        kmer = seq[i:i+k] ## set up kmer substring
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += seq[:k-length] ## start including nucleotides from the "start" of the seq if cyclic
        else:
            if len(kmer) != k:
                continue ## don't include substrings near the end of the seq if not cyclic
        if kmer in kmers: ## update kmers and counts in kmers dict
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers

## calculate jaccard index of a pair of lists
def get_jaccard_index(first_list, second_list):
    """Calculate Jaccard distance between two sequences, using a list of kmers from each sequence as input.

    Args:
        first_list (list): The first of two lists to get Jaccard distance for.
        second_list (list): The second of two lists to get Jaccard distance for.

    Returns:
        Float: Jaccard index for given strings.
    """
    union = len(set(first_list) | set(second_list))
    intersection = len(set(first_list) & set(second_list))
    jaccard_index = intersection/union

    return jaccard_index

## read in fasta files; check how many records in each file -- complete genome vs contigs;
## generate list of kmers for each file
def read_fasta_and_generate_kmers(file):
    """Read in a fasta file, check if it contains one record (then cyclic = args.cyclic) or multiple records (then cyclic = False),
    then generate kmers from sequence(s) in the given file.

    Args:
        file (filepath): The fasta file to process.

    Returns:
        list: list of kmers generated from the sequence(s) in the given file.
    """
    id_list = []
    kmer_list = []
    for record in sio.parse(file, "fasta"):
        id_list.append(record.id)

    if len(id_list) > 1:
        for record in sio.parse(file, "fasta"):
            print("Working with sequence of length ", len(record.seq), "in file ", file, "with cyclic = ", False) ## checkpoint -- should give multiple lengths per file
            kmer_dict = get_kmers(str(record.seq), kmer_size, False)
            kmer_list_multiple = [str(kmer) for kmer in list(kmer_dict.keys())] ## temporary list, gets overwritten for each record
            kmer_list += kmer_list_multiple ## master list, updates with items from each of the lists above
    else:
        for record in sio.parse(file, "fasta"):
            print("Working with sequence of length ", len(record.seq), "in file ", file, "with cyclic = ", cyclic) ## checkpoint -- should give single length value per file
            kmer_dict = get_kmers(str(record.seq), kmer_size, cyclic)
            kmer_list = [str(kmer) for kmer in list(kmer_dict.keys())]
    print(f"#kmers in {file} = {len(kmer_list)}") ## checkpoint -- should be a number around 2,000,000 for these samples
    return kmer_list

def read_folder_and_generate_kmers(some_dir, minhash):
    file_dict = {}
    for fasta in os.listdir(some_dir):
        full_path = os.path.join(os.path.abspath(some_dir), fasta)
        ## make a dict containing file basename as the key and a list of kmers as the first element in an otherwise empty list of values.
        file_dict[fasta] = [read_fasta_and_generate_kmers(os.path.abspath(full_path))]
        ## if minhash is set to true, also create hash list and sketch for each set of kmers.
        if minhash:
            for key, values in file_dict.items():
                file_dict[key].append(sorted([mmh3.hash(kmer) for kmer in values[0]]))
            for key, values in file_dict.items():
                file_dict[key].append(sorted(values[1])[0:sketch_size])
    return file_dict


if not fasta_dir: ## if two fasta files given as input
    kmers_file1 = read_fasta_and_generate_kmers(fasta_1)
    kmers_file2 = read_fasta_and_generate_kmers(fasta_2)

    jacc_index_nohash = get_jaccard_index(kmers_file1, kmers_file2)
    print("Jaccard index for sequences in the given files is: ", jacc_index_nohash)
    print("Jaccard distance for sequences in the given files is: ", (1-jacc_index_nohash))
    if minhash:
        hashes_file1 = sorted([mmh3.hash(kmer) for kmer in kmers_file1])
        hashes_file2 = sorted([mmh3.hash(kmer) for kmer in kmers_file2])
        sketch_file1 = hashes_file1[0:sketch_size]
        sketch_file2 = hashes_file2[0:sketch_size]
        jacc_index_withhash = get_jaccard_index(sketch_file1, sketch_file2)
        print("Jaccard index after hashing+sketching for sequences in the given files is: ", jacc_index_withhash)
        print("Jaccard distance after hashing+sketching for sequences in the given files is: ", (1-jacc_index_withhash))
else: ## if user chooses to provide a directory containing two or more fasta files
    print("Processing a directory...")
    file_kmers_dict = read_folder_and_generate_kmers(fasta_dir, minhash)

    sample_list = list(file_kmers_dict.keys())
    pairs = list(combinations([sample for sample in sample_list], 2))

    for pair in pairs:
        print(f"Jaccard index for {pair[0]} & {pair[1]} is: ", get_jaccard_index(file_kmers_dict[pair[0]][0], file_kmers_dict[pair[1]][0]))
        if minhash:
            print(f"Jaccard index for {pair[0]} & {pair[1]} after minhash  is: ", get_jaccard_index(file_kmers_dict[pair[0]][2], file_kmers_dict[pair[1]][2]))

