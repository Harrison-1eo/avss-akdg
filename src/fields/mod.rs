pub mod fp64;
// mod bigint;

pub trait Field: Sized + Clone + Copy
    + std::ops::Neg<Output = Self>
    + std::ops::Add<Output = Self>
    + std::ops::AddAssign
    + std::ops::Sub<Output = Self>
    + std::ops::SubAssign
    + std::ops::Mul<Output = Self>
    + std::ops::MulAssign
    + std::cmp::PartialEq
    + std::fmt::Display {
        fn from_int(x: u64) -> Self;
        fn random_element() -> Self;
        fn inverse(&self) -> Self;
        fn pow(&self, n: u64) -> Self;
        fn is_zero(&self) -> bool;
        fn get_generator(order: u64) -> Self;
}
