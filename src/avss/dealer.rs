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
    evaluations: Vec<MultilinearPolynomial<T>>,
}

impl<T: Field + 'static> Dealer<T> {
    fn fold(values: &Vec<T>, parameter: T, coset: &Coset<T>) -> Vec<T> {
        let len = values.len() / 2;
        let res = (0..len)
            .into_iter()
            .map(|i| {
                let x = values[i];
                let nx = values[i + len];
                let new_v = (x + nx) + parameter * (x - nx) * coset.element_inv_at(i);
                new_v * T::INVERSE_2
            })
            .collect();
        res
    }

    fn batch_folding(
        total_round: usize,
        polynomial: &MultilinearPolynomial<T>,
        folding_parameter: &Vec<Vec<T>>,
        coset: &Vec<Coset<T>>,
    ) -> (Vec<Vec<Vec<T>>>, Vec<MultilinearPolynomial<T>>) {
        let mut res = vec![vec![(coset[0].fft(polynomial.coefficients()))]];
        let variable_num = polynomial.variable_num();
        let mut evaluations = vec![];
        for round in 0..total_round {
            let len = res[round].len();
            if round < total_round - 1 {
                let mut evaluations = vec![];
                for (index, j) in folding_parameter[round].iter().enumerate() {
                    let next_evaluation =
                        Self::fold(&res[round][index & (len - 1)], *j, &coset[round]);
                    evaluations.push(next_evaluation);
                }
                res.push(evaluations);
            } else {
                for (index, j) in folding_parameter[round].iter().enumerate() {
                    let next_evaluation =
                        Self::fold(&res[round][index & (len - 1)], *j, &coset[round]);
                    let mut coefficients = coset[round + 1].ifft(&next_evaluation);
                    coefficients.truncate(1 << (variable_num - total_round));
                    evaluations.push(MultilinearPolynomial::new(coefficients));
                }
            }
        }
        (res, evaluations)
    }

    pub fn new(
        total_round: usize,
        polynomial: &MultilinearPolynomial<T>,
        interpolate_coset: &Vec<Coset<T>>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
        folding_parameter: &Vec<Vec<T>>,
    ) -> Self {
        let (functions, evaluations) = Self::batch_folding(
            total_round,
            polynomial,
            folding_parameter,
            interpolate_coset,
        );
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
            avss_party[i].set_share(&self.evaluations[i % self.evaluations.len()]);
        }
    }

    pub fn prove(&mut self) {
        self.prover.prove();
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        self.prover.query()
    }
}
