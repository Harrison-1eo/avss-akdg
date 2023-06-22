pub mod fp64;
pub mod mersenne61_ext;

pub trait Field:
    Sized
    + Clone
    + Copy
    + std::ops::Neg<Output = Self>
    + std::ops::Add<Output = Self>
    + std::ops::AddAssign
    + std::ops::Sub<Output = Self>
    + std::ops::SubAssign
    + std::ops::Mul<Output = Self>
    + std::ops::MulAssign
    + std::cmp::PartialEq
    + std::fmt::Display
    + std::fmt::Debug
{
    const LOG_ORDER: u64;
    const ROOT_OF_UNITY: Self;

    fn from_int(x: u64) -> Self;
    fn random_element() -> Self;
    fn inverse(&self) -> Self;
    fn is_zero(&self) -> bool;

    fn get_generator(order: usize) -> Self {
        if (order & (order - 1)) != 0 || order > (1 << Self::LOG_ORDER) {
            panic!("invalid order");
        }
        let mut res = Self::ROOT_OF_UNITY;
        let mut i = 1u64 << Self::LOG_ORDER;
        while i > order as u64 {
            res *= res;
            i >>= 1;
        }
        res
    }

    fn pow(&self, mut n: u64) -> Self {
        let mut ret = Self::from_int(1);
        let mut base = self.clone();
        while n != 0 {
            if n % 2 == 1 {
                ret *= base;
            }
            base *= base;
            n >>= 1;
        }
        ret
    }
}

mod field_tests {
    use super::*;

    pub fn add_and_sub<T: Field>() {
        for _i in 0..100 {
            let a = T::random_element();
            let b = T::random_element();
            let c = a + b - a;
            assert!(b == c)
        }
    }

    pub fn mult_and_inverse<T: Field>() {
        for _i in 0..100 {
            let a = T::random_element();
            let b = a.inverse();
            assert_eq!(a * b, T::from_int(1));
            assert_eq!(b * a, T::from_int(1));
        }
    }

    pub fn assigns<T: Field>() {
        for _i in 0..10 {
            let mut a = T::random_element();
            let aa = a;
            let b = T::random_element();
            a += b;
            assert_eq!(a, aa + b);
            a -= b;
            assert_eq!(a, aa);
            a *= b;
            assert_eq!(a, aa * b);
            a *= b.inverse();
            assert_eq!(a, aa);
            assert!((-a + a).is_zero());
        }
    }

    pub fn pow_and_generator<T: Field>() {
        assert_eq!(T::get_generator(1), T::from_int(1));
        let x = T::get_generator(1 << 32);
        assert_eq!(x.pow(1 << 32), T::from_int(1));
        assert_ne!(x.pow(1 << 31), T::from_int(1));
    }
}
