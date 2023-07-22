pub mod prover;
pub mod verifier;

use util::algebra::field::Field;

#[derive(Debug, Clone)]
pub struct Tuple<T: Field> {
    pub a: T,
    pub b: T,
    pub c: T,
}

impl<T: Field> Tuple<T> {
    pub fn verify(&self, beta: T, folding_param: T) -> bool {
        let v = beta * (self.a + self.b) + folding_param * (self.a - self.b);
        v * T::INVERSE_2 == self.c * beta
    }
}