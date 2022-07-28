# CONVERGRAPH
Highly experimental tool to create a mutation co-occurrence graph for viewing in tools like GEPHI. Ulimate goal is to find convergently evolved shared mutations.

## Build the code

- [Install Rust](https://www.rust-lang.org/tools/install)
- From the project directory: `cargo build --release`
- Optionally, you can try to speed it up by specifying the architecture before compiling. Examples: 
  - `export RUSTFLAGS="-C target-cpu=haswell"`
  - `export RUSTFLAGS="-C target-cpu=native"`
  - `export RUSTFLAGS="-C target-cpu=x86-64-v3"`

## Run the code
Current way to run:

```{bash}
# Mac OS
zcat < data/input-data.1.txt.gz  |time -l target/release/convergraph  > data/graph.dot

# Linux
zcat data/input-data.1.txt.gz  |time -v target/release/convergraph  > data/graph.dot
```

## Description of the Method

Algorithm Input: 
- A set of aligned amino acid sequences spanning at least two clades or lineages (MSA)
- An outgroup sequence in the same aligned coordinate space that is ancestral to the data in the MSA
- Parameters
  - *conservation threshold* - the lower bound frequency required to exempt sites from the calculation
  - *minimum co-occurence support* - the count at this edges are removed the graph
  - *minimum co-occurrence frequency* - the frequency at which edges are removed the graph

Creating the initial graph:
1. Calculate frequencies for all amino acids at each site in the MSA. 
2. Remove sites that are *conserved*, ie, sites where the most frequent residue exceeds the conservation threshold (default 97%)
3. For each query sequence, calculate all amino acid substitutions for sites found in (2) versus the supplied outgroup sequence.
4. Add substitutions as nodes in an undirected graph.
5. Add edges where substitutions co-occur in the same query sequence.
6. Weight edges by the co-occurrence counts.

Pruning the graph:
1. Remove edges that fail to meet the minimum co-occurrence threshold (default: 4)
2. Remove edges that fail to meet the minimum co-occurrence frequency (default: 10%)
3. Remove orphaned nodes (no edges)
4. Output pruned graph
