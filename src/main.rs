// Filename:         convergraph
// Description:      Highly experimental tool to create a mutation co-occurrence
//                   graph for viewing in tools like GEPHI. Ulimate goal is to
//                   find convergently evolved shared mutations.
//
// Date dedicated:   2022-07-27
// Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
//
// Citation:         Unpublished
//
// =============================================================================
//
//                            PUBLIC DOMAIN NOTICE
//
//  This source code file or script constitutes a work of the United States
//  Government and is not subject to domestic copyright protection under 17 USC §
//  105. This file is in the public domain within the United States, and
//  copyright and related rights in the work worldwide are waived through the CC0
//  1.0 Universal public domain dedication:
//  https://creativecommons.org/publicdomain/zero/1.0/
//
//  The material embodied in this software is provided to you "as-is" and without
//  warranty of any kind, express, implied or otherwise, including without
//  limitation, any warranty of fitness for a particular purpose. In no event
//  shall the Centers for Disease Control and Prevention (CDC) or the United
//  States (U.S.) government be liable to you or anyone else for any direct,
//  special, incidental, indirect or consequential damages of any kind, or any
//  damages whatsoever, including without limitation, loss of profit, loss of
//  use, savings or revenue, or the claims of third parties, whether or not CDC
//  or the U.S. government has been advised of the possibility of such loss,
//  however caused and on any theory of liability, arising out of or in
//  connection with the possession, use or performance of this software.
//
//  Please provide appropriate attribution in any work or product based on this
//  material.

use std::io;

use clap::Parser;
use csv::ReaderBuilder;
use petgraph::data::Build;
use serde::Deserialize;
use std::path::PathBuf;

const ALPHA_LENGTH: usize = (b'Z' - b'A') as usize + 1 + 3;
const AA_DELETE: usize = ALPHA_LENGTH - 3;
const AA_STOP: usize = ALPHA_LENGTH - 2;
const AA_ELSE: usize = ALPHA_LENGTH - 1;

fn read_records(amino_acid_sequence: &mut Vec<Vec<u8>>, has_header: bool) {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(io::stdin());

    if has_header {
        let _ = rdr.headers().expect("Input file is missing headers!\n");
    }

    for result in rdr.deserialize() {
        let record: Result<Record, _> = result;
        match record {
            Ok(r) => amino_acid_sequence.push(r.aa_aln.as_bytes().to_vec()),
            Err(e) => panic!("{}", e),
        }
    }
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct Record {
    cds_id: String,
    accession: String,
    date_first_seen: String,
    strain_count: String,
    country_first_seen: String,
    aa_aln: String,
    cds_aln: String,
}

#[derive(Parser)]
#[clap(author = "Samuel S. Shepard, Centers for Disease Control and Prevention (vfn4@cdc.gov)")]
#[clap(version = "0.1")]
#[clap(
    about = "Highly experimental tool to create a mutation co-occurrence graph for viewing in tools like GEPHI. Ulimate goal is to find convergently evolved shared mutations."
)]
#[clap(long_about = None)]
struct Cli {
    #[clap(long, short, value_parser, value_name = "Reference File")]
    reference_file: PathBuf,
    #[clap(
        default_value_t = 4,
        long,
        short = 's',
        value_parser,
        value_name = "minimum co-ocurrence support"
    )]
    minimum_coocurrence_support: u32,
    #[clap(
        default_value_t = 0.10,
        long,
        short = 'f',
        value_parser,
        value_name = "minimum co-occrrence frequency"
    )]
    minimum_cooccrrence_frequency: f32,
    #[clap(
        default_value_t = 0.97,
        long,
        short,
        value_parser,
        value_name = "conservation threshold"
    )]
    conservation_threshold: f32,
    #[clap(short, long, short, action)]
    query_has_header: bool,
}

fn main() {
    let args = Cli::parse();

    // MT019531 / WUHAN
    let ref_sequence = std::fs::read_to_string(args.reference_file).expect("Bad reference file");
    let ref_sequence = ref_sequence.as_bytes();
    let minimum_coocurrence_support = args.minimum_coocurrence_support;
    let minimum_cooccrrence_frequency: f32 = args.minimum_cooccrrence_frequency;
    let conservation_threshold: f32 = args.conservation_threshold;

    let mut sequences: Vec<Vec<u8>> = vec![];
    read_records(&mut sequences, args.query_has_header);

    // Maximum length of sequences, should be uniform but this is required to avoid issues
    let seq_len = sequences
        .iter()
        .map(|s| s.len())
        .max()
        .unwrap_or(0 as usize);

    let number_sequences = sequences.len();
    eprintln!("Data are {number_sequences} x {seq_len}");

    // Amino acid count table
    let mut counts: Vec<[u32; ALPHA_LENGTH]> = vec![[0u32; ALPHA_LENGTH]; seq_len];

    for s in sequences.iter() {
        for (i, b) in s.iter().enumerate() {
            match b {
                b'A'..=b'Z' => counts[i][(b - b'A') as usize] += 1,
                b'a'..=b'z' => counts[i][(b - b'a') as usize] += 1,
                b'-' => counts[i][AA_DELETE] += 1,
                b'*' => counts[i][AA_STOP] += 1,
                _ => counts[i][AA_ELSE] += 1,
            }
        }
    }

    // Calculate conservation and filter
    let mut valid_positions: Vec<usize> = vec![];
    for (i, d) in counts.iter().enumerate() {
        let sum: u32 = d.iter().sum();

        let mut max_value = 0;
        let mut max_index = None;

        for index in 0..d.len() {
            if counts[i][index] > max_value {
                max_value = counts[i][index];
                max_index = Some(index);
            }
        }

        if let Some(index) = max_index {
            let aa = index_to_aa(index);

            let freq: f32 = counts[i][index] as f32 / sum as f32;
            if freq < conservation_threshold {
                valid_positions.push(i);
                let p = i + 1;
                eprintln!("{p:0>4} / {aa}: {freq:.4} ({sum})");
            }
        }
    }

    use itertools::Itertools;
    use petgraph::graphmap::GraphMap;

    // Add nodes and edges to graph
    let mut aa_mut_net: GraphMap<NodeSub, u32, petgraph::Undirected> = GraphMap::new();
    for s in sequences {
        let mut nodes: Vec<NodeSub> = Vec::new();
        for ptr in valid_positions.iter() {
            let i = *ptr;
            if i >= s.len() || i >= ref_sequence.len() {
                continue;
            }

            if ref_sequence[i] != s[i] {
                let aa_sub = NodeSub::new(i, ref_sequence[i], s[i]);
                aa_mut_net.add_node(aa_sub);
                nodes.push(aa_sub);
            }
        }

        nodes
            .iter()
            .tuple_combinations::<(_, _)>()
            .for_each(|(a, b)| {
                if aa_mut_net.contains_edge(*a, *b) {
                    let e = *aa_mut_net.edge_weight(*a, *b).unwrap();
                    aa_mut_net.update_edge(*a, *b, e + 1);
                } else {
                    aa_mut_net.add_edge(*a, *b, 1_u32);
                }
            });
    }

    // remove nodes with lack of support or frequency
    let nodes: Vec<NodeSub> = aa_mut_net.nodes().collect();
    for node in nodes.iter() {
        let edges: Vec<(NodeSub, NodeSub, u32)> = aa_mut_net
            .edges(*node)
            .map(|(a, b, w)| (a, b, *w))
            .collect();

        for (a, b, w) in edges {
            let f = w as f32 / number_sequences as f32;
            if w < minimum_coocurrence_support || f < minimum_cooccrrence_frequency {
                aa_mut_net.remove_edge(a, b);
            }
        }
    }

    // remove unconnected nodes
    for node in nodes.iter() {
        let neighbor_count = aa_mut_net.neighbors(*node).count();
        if neighbor_count == 0 {
            aa_mut_net.remove_node(*node);
        }
    }

    // print results
    use petgraph::dot::Dot;
    use regex::Regex;

    // hack output: TO-DO, improve
    let re = Regex::new(r#"--\s*\d+\s*\[\s*label\s*=\s*"(\d+)""#).unwrap();
    let dot = format!("{:?}", Dot::new(&aa_mut_net));
    let dot = re.replace_all(&dot, "$0, weight=$1");
    println!("{dot}");

    /*
        let dot = re.replace_all(&dot, |caps: &Captures| {
        format!(
            "{}, weight={}",
            &caps[0],
            caps[1].parse::<f32>().unwrap_or_default() / number_sequences as f32
        )
    }); */
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash)]
struct NodeSub {
    index: u32,
    ancestral: u8,
    derived: u8,
}

impl std::fmt::Debug for NodeSub {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}",
            self.ancestral as char,
            (self.index + 1).to_string(),
            self.derived as char
        )
    }
}

impl std::fmt::Display for NodeSub {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}",
            self.ancestral as char,
            (self.index + 1).to_string(),
            self.derived as char
        )
    }
}

impl NodeSub {
    fn new(index: usize, reference: u8, query: u8) -> Self {
        NodeSub {
            index: index as u32,
            ancestral: reference,
            derived: query,
        }
    }
}

fn index_to_aa(index: usize) -> char {
    match index {
        AA_DELETE => '-',
        AA_STOP => '*',
        AA_ELSE => '?',
        _ => (index as u8 + b'A') as char,
    }
}
