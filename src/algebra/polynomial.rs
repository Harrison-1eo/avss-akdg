use super::coset::Coset;
use super::field::Field;

pub struct Polynomial<T: Field> {
    coefficients: Vec<T>,
}

impl<T: Field> Polynomial<T> {
    pub fn new(mut coefficients: Vec<T>) -> Polynomial<T> {
        let zero = T::from_int(0);
        while *coefficients.last().unwrap() == zero {
            coefficients.pop();
        }
        Polynomial { coefficients }
    }

    pub fn coefficients(&self) -> &Vec<T> {
        &self.coefficients
    }

    pub fn random_polynomial(degree: usize) -> Polynomial<T> {
        Polynomial {
            coefficients: (0..degree).map(|_| Field::random_element()).collect(),
        }
    }

    pub fn degree(&self) -> usize {
        let n = self.coefficients.len();
        if n == 0 {
            0
        } else {
            n - 1
        }
    }

    pub fn evaluation_at(&self, x: T) -> T {
        let mut res = Field::from_int(0);
        for i in self.coefficients.iter().rev() {
            res *= x;
            res += *i;
        }
        res
    }

    pub fn evaluation_over_coset(&self, coset: &Coset<T>) -> Vec<T> {
        coset.fft(&self.coefficients)
    }
}

struct VanishingPolynomial<T: Field> {
    degree: usize,
    shift: T,
}

impl<T: Field> VanishingPolynomial<T> {
    fn new(coset: &Coset<T>) -> VanishingPolynomial<T> {
        let degree = coset.num_elements();
        VanishingPolynomial {
            degree,
            shift: coset.shift().pow(degree as u64),
        }
    }

    // The n roots of the equation x^n - a^n = 0 are a*w_n^0, ..., a*w_n*{n-1}
    // Thus, f(x) = (x - a*w_n^0)...(x - a*w_n^{n-1}) = x^n - a^n
    fn evaluation_at(&self, x: T) -> T {
        x.pow(self.degree as u64) - self.shift
    }
}

#[cfg(test)]
mod test {
    use super::super::field::fp64::Fp64;
    use super::*;

    #[test]
    fn evaluation() {
        let coset = Coset::new(32, Fp64::random_element());
        let all_elements = coset.all_elements();
        let poly = Polynomial::random_polynomial(32);
        let eval = poly.evaluation_over_coset(&coset);
        for i in 0..coset.num_elements() {
            assert_eq!(eval[i], poly.evaluation_at(all_elements[i]));
        }
        let poly = VanishingPolynomial::new(&coset);
        for i in coset.all_elements() {
            assert_eq!(Fp64::from_int(0), poly.evaluation_at(i));
        }
    }
}
