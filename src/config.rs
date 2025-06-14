use num_complex::Complex32;

pub type QubitIdx = usize;
pub type GateIdx = usize;
pub type Real = f32;
pub type Complex = Complex32;

pub mod constants {
    pub const RECP_SQRT_2: super::Real = std::f32::consts::FRAC_1_SQRT_2;
    pub const ZERO_THRESHOLD: super::Real = 0.00000001;

    pub const BOND_DIM_THRESHOLD: usize = 100;
}
