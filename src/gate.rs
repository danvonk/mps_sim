mod unitary;
mod defn;

use crate::config::{QubitIdx};
use defn::GateDefn;


#[derive(Debug)]
pub struct Gate {
    pub defn: GateDefn,
    pub touches: Vec<QubitIdx>,
}

fn create_touches(defn: &GateDefn) -> Vec<QubitIdx> {
    match *defn {
        GateDefn::Hadamard(qi)
        | GateDefn::PauliY(qi)
        | GateDefn::PauliZ(qi)
        | GateDefn::Phase { target: qi, .. }
        | GateDefn::S(qi)
        | GateDefn::Sdg(qi)
        | GateDefn::SqrtX(qi)
        | GateDefn::SqrtXdg(qi)
        | GateDefn::T(qi)
        | GateDefn::Tdg(qi)
        | GateDefn::X(qi) => vec![qi],
        GateDefn::CPhase {
            control, target, ..
        }
        | GateDefn::CZ { control, target }
        | GateDefn::CX { control, target } => vec![control, target],
        GateDefn::CCX {
            control1,
            control2,
            target,
        } => vec![control1, control2, target],
        GateDefn::FSim { left, right, .. } => vec![left, right],
        GateDefn::RX { target, .. } | GateDefn::RY { target, .. } | GateDefn::RZ { target, .. } => {
            vec![target]
        }
        GateDefn::CSwap {
            control,
            target1,
            target2,
        } => vec![control, target1, target2],
        GateDefn::Swap { target1, target2 } => vec![target1, target2],
        GateDefn::U { target, .. } => vec![target],
        GateDefn::Other { .. } => vec![],
    }
}



impl Gate {
    pub fn new(defn: GateDefn) -> Self {
        let touches = create_touches(&defn);
        Self { defn, touches }
    }
}