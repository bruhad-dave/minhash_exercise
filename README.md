# minhash_exercise
Code for a simple minhash implementation.  

## Usage
### Compare two fasta files:
```
## without hashing  
python3 ./kmer_distances.py -f1 <path/to/fasta1.fa> -f2 <path/to/fasta2.fa> -k kmer-size  

## with minhash  
python3 ./kmer_distances.py -f1 <path/to/fasta1.fa> -f2 <path/to/fasta2.fa> -k kmer-size -m  

## see how minhash distances compare to full distance at different sketch sizes  
python3 ./kmer_distances.py -f1 <path/to/fasta1.fa> -f2 <path/to/fasta2.fa> -k kmer-size -m -v True  
```

### Compare fasta files in a directory (ideally three or more):
```
## without hashing  
python3 ./kmer_distances.py -d <path/to/fasta/dir/> -k kmer-size  

## with minhash  
python3 ./kmer_distances.py -d <path/to/fasta/dir/> -k kmer-size -m  

## make a neighbour-joining tree with full distances and one with minhash distances using the given sketch size  
python3 ./kmer_distances.py -d <path/to/fasta/dir/> -k kmer-size -m -t True  
## if path/to/fasta/dir/ contains only two fasta files, the script will do everything except make the tree
```
### All options:
```
-m/--minhash : Use minhash? Flips to True if option is used, False if nothing specified. Takes no arguments.  
-d/--fasta_dir : Path to directory containing fasta files (which can be .fa/.fasta/.faa/.fna). Not compatible with -f1 and -f2.  
-f1/fasta_1 : Path to the first of two fasta files to compare. Incompatible with -d.  
-f2/fasta_2 : Path to the second of two fasta files to compare. Incompatible with -d.  
-k/--kmer_size : Size of kmer to use (int). Required, no default.  
-s/--sketch_size : Sketch size to use for minhash (int). Defaults to 1000.  
-c/--cyclic : Consider sequenes in fasta inputs as cyclic (bool)? Defaults to False. Always set to False if a fasta file has >1 record even if user sets to True.  
-t/--make_tree : Make a neighbour-joining tree (bool)? Defaults to False. Only works when -d used unless specified directory contains only 2 files.  
-v/--var_sketchsize : Plot Jaccard distance at varying sketch sizes and compare to full distance (bool). Defaults to False. Only works when -f1 & -f2 are used.  
```

### Imports:
- combinations (from itertools), os, argparse  
- numpy, pandas  
- scikit.bio (to make the NJ tree)  
- matplotlib (to make plot for minhash distances with varying sketch sizes)  
- mmh3 (hash function)  
- Bio.SeqIO (to work with fasta files)  

- [x] Steps 1 through 5  
- [x] Output sketches to file  
- [x] A plot to represent that minhash distances with progressively larger sketch sizes converge with full distances  
- [x] Neighbour-joining tree  
