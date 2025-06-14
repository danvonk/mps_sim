mod circuit;
mod config;
mod fingerprint;
mod gate;
mod gate_scheduler;
mod options;
mod parser;
mod simulator;

use circuit::Circuit;
use options::Options;
use structopt::StructOpt;
use config::{Complex, constants};
use crate::simulator::Compactifiable;
use crate::fingerprint::Fingerprint;

pub fn print_complex(c: &Complex) -> String {
    if c.im > -constants::ZERO_THRESHOLD {
        format!("{:.8}+{:.8}i", c.re, c.im.abs(),)
    } else {
        format!("{:.8}-{:.8}i", c.re, c.im.abs(),)
    }
}

fn main() {
    let options = Options::from_args();
    let source = std::fs::read_to_string(&options.input).unwrap();
    let program = match parser::parse_program(&source) {
        Ok(program) => program,
        Err(_) => {
            panic!("Failed to parse program.");
        }
    };

    let circuit = match Circuit::new(program) {
        Ok(circuit) => circuit,
        Err(err) => {
            panic!("Failed to build circuit.");
        }
    }
    .decompose();

    let result = simulator::mps_simulator::run(circuit).compactify();

    let mut fingerprint = Fingerprint::new(10);

    for (bidx, weight) in result {
        fingerprint.insert(bidx.clone(), weight);
    }

    println!("computed fingerprint:");
    fingerprint
        .iter()
        .enumerate()
        .for_each(|(idx, (bidx, weight))| {
            println!(
                "fp{idx} {:0width$} {}",
                bidx,
                print_complex(&weight),
                width = 64,
            );
        });
}
