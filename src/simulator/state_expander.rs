use crate::{gate::Gate, simulator::state::MPS};

pub struct ExpandResult<'a> {
    pub state: &'a mut MPS,
    pub num_nonzeros: usize,
    pub num_gate_apps: usize,
}

pub fn expand<'a>(
    gates: Vec<&'a Gate>,
    num_qubits: usize,
    prev_num_nonzeros: usize,
    mps: &'a mut MPS,
) -> ExpandResult<'a> {
    let mut gate_apps = 0;

    for g in gates {
        mps.apply_gate(g);
        gate_apps += 1;
    }

    ExpandResult {
        state: mps,
        num_gate_apps: gate_apps,
        num_nonzeros: 0,
    }
}