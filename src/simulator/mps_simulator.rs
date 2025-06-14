use crate::config::{Complex, Real};
use nalgebra::*;
use crate::profile;

use crate::circuit::Circuit;
use crate::gate_scheduler::{GateScheduler, create_gate_scheduler};
use crate::simulator::{state::MPS, state_expander::ExpandResult, state_expander::expand};

#[macro_export]
macro_rules! profile {
    ($($exp:expr)+) => {
        {
            let _instant = std::time::Instant::now();
            let _result = {
                $($exp)+
            };
            let _duration = _instant.elapsed();

            (_duration, _result)
        }
    }
}

pub fn run(circuit: Circuit) -> MPS {
    let num_gates = circuit.num_gates();
    let num_qubits = circuit.num_qubits;
    let mut num_nonzeros = 1;
    let mut prev_num_nonzeros = 1;

    let mut num_gates_visited = 0;
    let mut state = MPS::singleton(num_qubits);
    let mut num_gate_apps = 0;

    let mut gate_scheduler = create_gate_scheduler(&circuit);

    let (duration, _) = profile!(loop {
        let these_gates = gate_scheduler
            .pick_next_gates()
            .into_iter()
            .map(|idx| &circuit.gates[idx])
            .collect::<Vec<_>>();

        log::debug!("applying gates: {:?}", these_gates);

        if these_gates.is_empty() {
            break;
        }

        let num_gates_visited_here = these_gates.len();

        let (
            duration,
            ExpandResult {
                num_nonzeros: new_num_nonzeros,
                num_gate_apps: num_gate_apps_here,
                ..
            },
        ) = profile!(expand(these_gates, num_qubits, prev_num_nonzeros, &mut state));

        let density = {
            let max_num_states: u64 = 1 << num_qubits;
            num_nonzeros as Real / max_num_states as Real
        };

        let throughput = (num_gate_apps_here as Real / 1e6) / duration.as_secs_f32();

        println!(
            "gate: {:<3} density: {:.8} nonzero: {:>10} hop: {:<2} time: {:.4}s throughput: {:.2}M gates/s",
            num_gates_visited,
            density,
            num_nonzeros,
            num_gates_visited_here,
            duration.as_secs_f32(),
            throughput
        );

        num_gates_visited += num_gates_visited_here;
        num_gate_apps += num_gate_apps_here;
        prev_num_nonzeros = num_nonzeros;
        num_nonzeros = new_num_nonzeros;
    });

    let final_density = {
        let max_num_states: u64 = 1 << num_qubits;
        num_nonzeros as f64 / max_num_states as f64
    };

    println!(
        "gate: {:<2} density: {:.8} nonzero: {:>10}\ngate app count: {}, time: {}s",
        num_gates_visited,
        final_density,
        num_nonzeros,
        num_gate_apps,
        duration.as_secs_f32()
    );

    assert!(num_gates_visited >= num_gates);
    state
}
