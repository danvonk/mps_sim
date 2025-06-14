pub mod mps_simulator;
pub mod state;
pub mod state_expander;

use std::fmt::{self, Display, Formatter};

use crate::config::Complex;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub struct BasisIdx {
    bits: u64
}

pub trait Compactifiable {
    fn compactify(self) -> Box<dyn Iterator<Item = (BasisIdx, Complex)>>;
}

impl BasisIdx {
    pub fn new(bits: &str) -> Self {
        Self {
            bits: u64::from_str_radix(bits, 2).unwrap(),
        }
    }

    fn get(&self, qi: usize) -> bool {
        self.bits & (1 << qi) != 0
    }

    fn flip(&self, qi: usize) -> Self {
        Self {
            bits: self.bits ^ (1 << qi),
        }
    }

    fn zeros() -> Self {
        Self { bits: 0 }
    }

    fn set(&self, qi: usize) -> Self {
        Self {
            bits: self.bits | (1 << qi),
        }
    }

    fn unset(&self, qi: usize) -> Self {
        Self {
            bits: self.bits & !(1 << qi),
        }
    }

    fn swap(&self, qi1: usize, qi2: usize) -> Self {
        let tmp = ((self.bits >> qi1) ^ (self.bits >> qi2)) & 1;

        Self {
            bits: self.bits ^ (tmp << qi1) ^ (tmp << qi2),
        }
    }

    pub fn from_idx(idx: usize) -> Self {
        Self { bits: idx as u64 }
    }

    pub fn as_idx(&self) -> usize {
        self.bits as usize
    }

    pub fn empty_key(_num_qubits: usize) -> Self {
        Self { bits: (1 << 63) }
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        self.bits.to_be_bytes().to_vec()
    }
}

impl Display for BasisIdx {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match f.width() {
            Some(width) => format!("{:0width$b}", self.bits, width = width).fmt(f),
            None => format!("{:b}", self.bits).fmt(f),
        }
    }
}

 

