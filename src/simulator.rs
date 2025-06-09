pub mod mps_simulator;

use crate::config::Complex;


pub struct Basis {
    bits: u64
}

pub trait Compactifiable {
    fn compactify(self) -> Box<dyn Iterator<Item = (Basis, Complex)>>;
}