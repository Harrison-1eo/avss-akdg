mod fp64;
// mod bigint;

pub trait Field: Sized
    + std::ops::Neg
    + std::ops::Add
    + std::ops::AddAssign
    + std::ops::Sub
    + std::ops::SubAssign
    + std::ops::Mul
    + std::ops::MulAssign
    + std::cmp::PartialEq
    + std::fmt::Display {
        fn random_element() -> Self;
        fn inverse(&self) -> Self;
        fn pow(&self, n: u64) -> Self;
        fn is_zero(&self) -> bool;
        fn get_unity_root() -> Self;
}
