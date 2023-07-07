use std::{cell::RefCell, rc::Rc};

use crate::rolling_fri::QueryResult;
use crate::{
    algebra::{coset::Coset, field::Field},
    one2many::{verifier::One2ManyVerifier, RandomOracle},
};

#[derive(Debug, Clone)]
pub struct Tuple<T: Field> {
    pub a: T,
    pub b: T,
    pub c: T,
}

impl<T: Field> Tuple<T> {
    pub fn verify(&self, beta: T, folding_param: T) -> bool {
        let v = self.a + self.b + folding_param * (self.a - self.b) * beta.inverse();
        v * T::from_int(2).inverse() == self.c
    }
}

pub struct AvssParty<T: Field, const N: usize> {
    pub verifier: Rc<RefCell<One2ManyVerifier<T, N>>>,
    open_point: Vec<T>,
    oracle: Rc<RefCell<RandomOracle<T>>>,
    tuples: Vec<Tuple<T>>,
}

impl<T: Field + 'static, const N: usize> AvssParty<T, N> {
    pub fn new(
        interpolate_coset: &Coset<T>,
        open_point: Vec<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
    ) -> AvssParty<T, N> {
        AvssParty {
            verifier: Rc::new(RefCell::new(One2ManyVerifier::new(
                interpolate_coset,
                oracle,
            ))),
            open_point,
            tuples: vec![],
            oracle: oracle.clone(),
        }
    }

    pub fn add_tuple(&mut self, tuple: &Tuple<T>) {
        self.tuples.push(tuple.clone());
    }

    pub fn verify_tuple(&mut self) -> bool {
        assert_eq!(self.tuples.len(), self.open_point.len());
        let beta = self.oracle.borrow().beta();
        for (tuple, folding_param) in self.tuples.iter().zip(self.open_point.iter()) {
            if !tuple.verify(beta, folding_param.clone()) {
                return false;
            }
        }
        let len = self.tuples.len();
        for i in 0..len {
            let a = self.tuples[i].a;
            let b = self.tuples[i].b;
            if i == 0 {
                self.verifier
                    .borrow_mut()
                    .set_map(Box::new(move |v: T, x: T, _challenge: T| {
                        (v - b) * (x + beta).inverse() + v + (v - a) * (x - beta).inverse()
                    }) as Box<dyn Fn(T, T, T) -> T>);
            } else {
                let c = self.tuples[i - 1].c;
                self.verifier
                    .borrow_mut()
                    .set_map(Box::new(move |v: T, x: T, _challenge: T| {
                        (v - c) * (x - beta * beta).inverse()
                            + (v - b) * (x + beta).inverse()
                            + v
                            + (v - a) * (x - beta).inverse()
                    }) as Box<dyn Fn(T, T, T) -> T>);
            }
        }
        true
    }

    pub fn verify(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
    ) -> bool {
        self.verifier
            .borrow()
            .verify(folding_proofs, function_proofs)
    }
}
