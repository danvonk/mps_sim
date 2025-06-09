
use crate::config::{QubitIdx, Real};

#[derive(Debug, Clone)]
pub enum GateDefn {
    CCX {
        control1: QubitIdx,
        control2: QubitIdx,
        target: QubitIdx,
    },
    CPhase {
        control: QubitIdx,
        target: QubitIdx,
        rot: Real,
    },
    CSwap {
        control: QubitIdx,
        target1: QubitIdx,
        target2: QubitIdx,
    },
    CX {
        control: QubitIdx,
        target: QubitIdx,
    },
    CZ {
        control: QubitIdx,
        target: QubitIdx,
    },
    FSim {
        left: QubitIdx,
        right: QubitIdx,
        theta: Real,
        phi: Real,
    },
    Hadamard(QubitIdx),
    PauliY(QubitIdx),
    PauliZ(QubitIdx),
    Phase {
        rot: Real,
        target: QubitIdx,
    },
    RX {
        rot: Real,
        target: QubitIdx,
    },
    RY {
        rot: Real,
        target: QubitIdx,
    },
    RZ {
        rot: Real,
        target: QubitIdx,
    },
    S(QubitIdx),
    Sdg(QubitIdx),
    SqrtX(QubitIdx),
    SqrtXdg(QubitIdx),
    Swap {
        target1: QubitIdx,
        target2: QubitIdx,
    },
    T(QubitIdx),
    Tdg(QubitIdx),
    U {
        target: QubitIdx,
        theta: Real,
        phi: Real,
        lambda: Real,
    },
    X(QubitIdx),
    Other {
        name: String,
        params: Vec<Real>,
        args: Vec<QubitIdx>,
    },
}


impl GateDefn {
    fn decompose_ccx(defn: &GateDefn) -> Vec<GateDefn> {
        match defn {
            GateDefn::CCX {
                control1,
                control2,
                target,
            } => vec![
                GateDefn::Hadamard(*target),
                // CNOT(control2 -> target)
                GateDefn::CX {
                    control: *control2,
                    target: *target,
                },
                GateDefn::Tdg(*target),
                // CNOT(control1 -> target)
                GateDefn::CX {
                    control: *control1,
                    target: *target,
                },
                GateDefn::T(*target),
                // CNOT(control2 -> target)
                GateDefn::CX {
                    control: *control2,
                    target: *target,
                },
                GateDefn::Tdg(*target),
                // CNOT(control1 -> target)
                GateDefn::CX {
                    control: *control1,
                    target: *target,
                },
                GateDefn::T(*control2),
                GateDefn::T(*target),
                GateDefn::Hadamard(*target),
            ],
            _ => vec![],
        }
    }

    pub fn decompose_cswap(gate: &GateDefn) -> Vec<GateDefn> {
        match *gate {
            GateDefn::CSwap {
                control,
                target1,
                target2,
            } => {
                let mut decomp = vec![GateDefn::CX {
                    control: target1,
                    target: target2,
                }];
                decomp.append(
                    &mut GateDefn::CCX {
                        control1: control,
                        control2: target2,
                        target: target1,
                    }
                    .decompose_gate(),
                );
                decomp.push(GateDefn::CX {
                    control: target1,
                    target: target2,
                });

                decomp
            }
            _ => vec![],
        }
    }

    pub fn decompose_gate(&self) -> Vec<GateDefn> {
        match self {
            GateDefn::CCX { .. } => GateDefn::decompose_ccx(self),
            GateDefn::CSwap { .. } => GateDefn::decompose_cswap(self),
            _ => vec![self.clone()],
        }
    }
}
