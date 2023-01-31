#[derive(Debug, Clone, Copy)]
pub struct Bigint<const N: usize> {
    data: [u64; N]
}

impl<const N: usize> Bigint<N> {
    pub fn from_int(x: u64) -> Bigint<N> {
        let mut res = Bigint {
            data: [0; N]
        };
        res.data[0] = x;
        res
    }

    fn from_str(s: &str) -> Bigint<N> {
        unimplemented!()
    }

    fn one() -> Bigint<N> {
        Bigint::from_int(1)
    }

    fn is_zero(&self) -> bool {
        for i in self.data {
            if i != 0 {
                return false;
            }
        }
        true
    }
}

impl<const N: usize> std::fmt::Display for Bigint<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.data[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_zero() {
        let a: Bigint<1> = Bigint::from_int(1);
        assert!(!a.is_zero());
        let a: Bigint<1> = Bigint::from_int(0);
        assert!(a.is_zero());
    }
}
