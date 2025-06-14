use crate::{circuit::Circuit, config::{GateIdx, QubitIdx}, gate::Gate};

pub struct GateScheduler<'a> {
    frontier: Vec<GateIdx>,
    num_gates: usize,
    num_qubits: usize,
    gate_touches: Vec<&'a [QubitIdx]>,
    gate_is_branching: Vec<bool>,
    max_branching_stride: usize,
    disable_gate_fusion: bool,
}

pub fn okay_to_visit(
    num_gates: usize,
    gate_touches: &[&[QubitIdx]],
    frontier: &[GateIdx],
    gi: GateIdx,
) -> bool {
    gi < num_gates && gate_touches[gi].iter().all(|qi| frontier[*qi] == gi)
}

pub fn mark_as_visit(
    num_gates: usize,
    gate_touches: &[&[QubitIdx]],
    frontier: &mut [GateIdx],
    gi: GateIdx,
) {
    log::debug!("visiting gate: {}", gi);
    assert!(okay_to_visit(num_gates, gate_touches, frontier, gi));
    for qi in gate_touches[gi] {
        let next = next_touch(num_gates, gate_touches, *qi, gi + 1);

        frontier[*qi] = next;
        log::debug!("updated frontier[{}] to {}", qi, frontier[*qi]);
    }
}

pub fn next_touch(
    num_gates: usize,
    gate_touches: &[&[QubitIdx]],
    qi: QubitIdx,
    gi: GateIdx,
) -> GateIdx {
    if gi >= num_gates {
        num_gates
    } else if gate_touches[gi].contains(&qi) {
        gi
    } else {
        next_touch(num_gates, gate_touches, qi, gi + 1)
    }
}

impl<'a> GateScheduler<'a> {
    pub fn pick_next_gates(&mut self) -> Vec<GateIdx> {
        if self.disable_gate_fusion {
            match self.visit_nonbranching() {
                Some(gi) => vec![gi],
                None => match self.visit_branching() {
                    Some(gi) => vec![gi],
                    None => vec![],
                },
            }
        } else {
            let mut num_branching_so_far = 0;
            let mut next_gates = Vec::<GateIdx>::new();

            while num_branching_so_far < self.max_branching_stride {
                next_gates.append(&mut self.visit_maximal_nonbranching_run());

                if let Some(next_gate) = self.visit_branching() {
                    num_branching_so_far += 1;
                    next_gates.push(next_gate);
                } else {
                    break;
                }
            }

            // each gate in next_gates should be marked as already visited
            assert!(next_gates.iter().all(|gi| !okay_to_visit(
                self.num_gates,
                &self.gate_touches,
                &self.frontier,
                *gi
            )));

            log::debug!("next gates: {:?}", next_gates);

            next_gates
        }
    }

    pub fn new(
        num_gates: usize,
        num_qubits: usize,
        gate_touches: Vec<&'a [QubitIdx]>,
        gate_is_branching: Vec<bool>,
        disable_gate_fusion: bool,
    ) -> Self {
        log::debug!(
            "initializing greedy nonbranching gate scheduler with {} gates and {} qubits",
            num_gates,
            num_qubits
        );
        let scheduler = Self {
            frontier: (0..num_qubits)
                .map(|qi| next_touch(num_gates, &gate_touches, qi, 0))
                .collect(),
            num_gates,
            num_qubits,
            gate_touches,
            gate_is_branching,
            max_branching_stride: 2,
            disable_gate_fusion,
        };

        assert_eq!(scheduler.frontier.len(), num_qubits);
        assert_eq!(scheduler.gate_touches.len(), num_gates);
        assert_eq!(scheduler.gate_is_branching.len(), num_gates);

        log::debug!("initial frontier: {:?}", scheduler.frontier);

        scheduler
    }

    fn visit_nonbranching(&mut self) -> Option<GateIdx> {
        let candidate = self
            .frontier
            .iter()
            .filter(|gi| {
                *gi < &self.num_gates
                    && !self.gate_is_branching[**gi]
                    && okay_to_visit(self.num_gates, &self.gate_touches, &self.frontier, **gi)
            })
            .nth(0)
            .cloned();

        if let Some(gi) = candidate {
            mark_as_visit(self.num_gates, &self.gate_touches, &mut self.frontier, gi);
        }
        candidate
    }

    fn visit_maximal_nonbranching_run(&mut self) -> Vec<GateIdx> {
        let mut non_branching_gates = Vec::new();

        loop {
            let mut selection = Vec::<GateIdx>::new();

            for qi in 0..self.num_qubits {
                loop {
                    let next_gi = self.frontier[qi];
                    if next_gi >= self.num_gates
                        || self.gate_is_branching[next_gi]
                        || !okay_to_visit(
                            self.num_gates,
                            &self.gate_touches,
                            &self.frontier,
                            next_gi,
                        )
                    {
                        break;
                    } else {
                        assert!(okay_to_visit(
                            self.num_gates,
                            &self.gate_touches,
                            &self.frontier,
                            next_gi
                        ));
                        mark_as_visit(
                            self.num_gates,
                            &self.gate_touches,
                            &mut self.frontier,
                            next_gi,
                        );
                        selection.push(next_gi);
                    }
                }
            }

            if selection.is_empty() {
                break;
            } else {
                non_branching_gates.append(&mut selection);
            }
        }
        non_branching_gates
    }

    fn visit_branching(&mut self) -> Option<GateIdx> {
        let result = self
            .frontier
            .iter()
            .filter(|gi| {
                *gi < &self.num_gates
                    && self.gate_is_branching[**gi]
                    && okay_to_visit(self.num_gates, &self.gate_touches, &self.frontier, **gi)
            })
            .nth(0)
            .copied();

        if let Some(gi) = result {
            mark_as_visit(self.num_gates, &self.gate_touches, &mut self.frontier, gi);
        }

        result
    }
}


pub fn create_gate_scheduler(
    circuit: &Circuit,
) -> Box<GateScheduler> {
    let num_gates = circuit.num_gates();
    let num_qubits = circuit.num_qubits;
    let gate_touches = circuit
        .gates
        .iter()
        .map(|gate| gate.touches.as_slice())
        .collect();
    let gate_is_branching = circuit
        .gates
        .iter()
        .map(|gate| gate.is_branching())
        .collect();

    Box::new(GateScheduler::new(
        num_gates,
        num_qubits,
        gate_touches,
        gate_is_branching,
        false
    ))
}
