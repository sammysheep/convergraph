# CONVERGRAPH

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

## Algorithm

TO-DO