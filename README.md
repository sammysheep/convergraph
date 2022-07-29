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
USAGE:
    convergraph [OPTIONS] --reference-file <Reference File> [QUERY_FILE (or use stdin)]

ARGS:
    <QUERY_FILE (or use stdin)>    

OPTIONS:
    -c, --conservation-threshold <conservation threshold>
            [default: 0.97]

    -f, --minimum-cooccrrence-frequency <minimum co-occrrence frequency>
            [default: 0.1]

    -q, --query-has-header
            

    -r, --reference-file <Reference File>
            

    -s, --minimum-coocurrence-support <minimum co-ocurrence support>
            [default: 4]

    -V, --version
            Print version information
```

## Description of the Method

Algorithm Input: 
- A set of aligned amino acid sequences spanning at least two clades or lineages (eg, an MSA or PSAs with shared reference). Usage refers to as `query file`. May have a header or not specified by `-q`. Records are unix line-delimited.
- An outgroup sequence in the same aligned coordinate space that is ancestral to the data in the *query* file. This is the *reference* file. Specified using: `-r file`.
- Parameters
  - *conservation threshold* - the lower bound frequency required to exempt sites from the calculation
  - *minimum co-occurence support* - the count at this edges are removed the graph
  - *minimum co-occurrence frequency* - the frequency at which edges are removed the graph

Creating the initial graph:
1. Calculate frequencies for all amino acids at each site in the *query file*. 
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
