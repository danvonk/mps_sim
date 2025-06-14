use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "mps_sim", about = "MPS Simulator")]
pub struct Options {
    #[structopt(
        parse(from_os_str),
        name = "input",
        short = "i",
        long = "input",
        help = "path to the input qasm file"
    )]
    pub input: PathBuf,

    #[structopt(
        parse(from_os_str),
        name = "output",
        short = "o",
        long = "output",
        help = "output file path in which to dump densities. if not specified, densities are not printed"
    )]
    pub output: Option<PathBuf>,
}