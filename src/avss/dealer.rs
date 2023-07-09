use std::cell::RefCell;
use std::rc::Rc;

use super::party::AvssParty;
use crate::algebra::coset::Coset;
use crate::random_oracle::RandomOracle;
use crate::util::QueryResult;
use crate::{
    algebra::{field::Field, polynomial::MultilinearPolynomial},
    one2many::prover::One2ManyProver,
};

pub struct Dealer<T: Field> {
    prover: One2ManyProver<T>,
    evaluations: Vec<T>,
}

impl<T: Field + 'static> Dealer<T> {
    fn fold(values: &Vec<T>, parameter: T, mut shift_inv: T, generator_inv: T) -> Vec<T> {
        let mut res = vec![];
        let len = values.len() / 2;
        for i in 0..len {
            let x = values[i];
            let nx = values[i + len];
            let new_v = (x + nx) + parameter * (x - nx) * shift_inv;
            res.push(new_v * T::from_int(2).inverse());
            shift_inv *= generator_inv;
        }
        res
    }

    fn batch_folding(
        polynomial: &MultilinearPolynomial<T>,
        folding_parameter: &Vec<Vec<T>>,
        coset: &Coset<T>,
    ) -> (Vec<Vec<(Vec<T>, Box<dyn Fn(T, T, T) -> T>)>>, Vec<T>) {
        let mut shift_inv = coset.shift().inverse();
        let mut generator_inv = coset.generator().inverse();
        let mut res = vec![vec![(
            coset.fft(polynomial.coefficients()),
            Box::new(move |v: T, _x: T, _c: T| v) as Box<dyn Fn(T, T, T) -> T>,
        )]];
        let mut evaluations = vec![];
        let total_round = folding_parameter.len();
        for (round, i) in folding_parameter.iter().enumerate() {
            if round < total_round - 1 {
                let mut polies = vec![];
                for (index, j) in i.iter().enumerate() {
                    polies.push((
                        Self::fold(
                            &res[round][index % res[round].len()].0,
                            *j,
                            shift_inv,
                            generator_inv,
                        ),
                        Box::new(move |v: T, _x: T, _c: T| v) as Box<dyn Fn(T, T, T) -> T>,
                    ));
                }
                res.push(polies);
                shift_inv *= shift_inv;
                generator_inv *= generator_inv;
            } else {
                for (index, j) in i.iter().enumerate() {
                    let values = &res[round][index % res[round].len()].0;
                    let x = values[0];
                    let nx = values[values.len() / 2];
                    let new_v = (x + nx) + (*j) * (x - nx) * shift_inv;
                    evaluations.push(new_v * T::INVERSE_2);
                }
            }
        }
        (res, evaluations)
    }

    pub fn new(
        polynomial: &MultilinearPolynomial<T>,
        interpolate_coset: &Coset<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
        folding_parameter: &Vec<Vec<T>>,
    ) -> Self {
        let total_round = polynomial.variable_num();
        let (functions, evaluations) =
            Self::batch_folding(polynomial, folding_parameter, interpolate_coset);
        Dealer {
            evaluations,
            prover: One2ManyProver::new(total_round, interpolate_coset, functions, oracle),
        }
    }

    pub fn commit_functions(&self, avss_party: &Vec<AvssParty<T>>) {
        let verifiers = avss_party.iter().map(|x| x.verifier.clone()).collect();
        self.prover.commit_functions(&verifiers);
    }

    pub fn commit_foldings(&self, avss_party: &Vec<AvssParty<T>>) {
        let verifiers = avss_party.iter().map(|x| x.verifier.clone()).collect();
        self.prover.commit_foldings(&verifiers);
    }

    pub fn send_evaluations(&self, avss_party: &mut Vec<AvssParty<T>>) {
        for i in 0..avss_party.len() {
            avss_party[i].set_share(self.evaluations[i]);
        }
    }

    pub fn prove(&mut self) {
        self.prover.prove();
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        self.prover.query()
    }
}
