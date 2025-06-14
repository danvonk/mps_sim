use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::config::{Complex};
use crate::simulator::BasisIdx;

#[derive(PartialEq)]
struct Entry(BasisIdx, Complex);

impl Eq for Entry {}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // FIXME: bucket magnitude?
        if self.1.norm() < other.1.norm() {
            Some(Ordering::Less)
        } else if self.1.norm() > other.1.norm() {
            Some(Ordering::Greater)
        } else {
            self.0.as_bytes().partial_cmp(&other.0.as_bytes())
        }
    }
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.1.norm() < other.1.norm() {
            Ordering::Less
        } else if self.1.norm() > other.1.norm() {
            Ordering::Greater
        } else {
            other.0.as_bytes().cmp(&self.0.as_bytes())
        }
    }
}

pub struct Fingerprint {
    num_entries: usize,
    entries: BTreeSet<Entry>,
}

impl Fingerprint {
    pub fn new(num_entries: usize) -> Self {
        Self {
            num_entries,
            entries: BTreeSet::new(),
        }
    }

    pub fn insert(&mut self, bidx: BasisIdx, weight: Complex) {
        self.entries.insert(Entry(bidx, weight));
        if self.entries.len() > self.num_entries {
            self.entries.pop_first();
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (BasisIdx, Complex)> + '_ {
        // reverse the comparison here so that the output fingerprint is sorted
        // like normal binary numbers, for readability
        self.entries
            .iter()
            .rev()
            .map(|Entry(bidx, weight)| (bidx.clone(), *weight))
    }
}
