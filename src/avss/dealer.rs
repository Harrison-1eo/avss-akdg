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

pub struct Dealer<T: Field, const N: usize> {
    prover: One2ManyProver<T, N>,
    evaluations: Vec<T>,
}

impl<T: Field + 'static, const N: usize> Dealer<T, N> {
    fn batch_folding(
        polynomial: &MultilinearPolynomial<T>,
        folding_parameter: &Vec<Vec<T>>,
    ) -> Vec<Vec<MultilinearPolynomial<T>>> {
        let mut res = vec![vec![polynomial.clone()]];
        for (round, i) in folding_parameter.iter().enumerate() {
            let mut polies = vec![];
            for (index, j) in i.iter().enumerate() {
                polies.push(res[round][index % res[round].len()].folding(*j));
            }
            res.push(polies);
        }
        res
    }

    pub fn calculate_functions(
        functions: &Vec<Vec<MultilinearPolynomial<T>>>,
        coset: &Coset<T>,
    ) -> Vec<Vec<(Vec<T>, Box<dyn Fn(T, T, T) -> T>)>> {
        let mut coset = coset.clone();
        let mut res = vec![];
        for i in 0..(functions.len() - 1) {
            let mut function_round = vec![];
            for j in 0..functions[i].len() {
                let coeff = functions[i][j].coefficients();
                if i == 0 {
                    function_round.push((
                        coset.fft(coeff),
                        Box::new(move |v: T, _x: T, _c: T| v) as Box<dyn Fn(T, T, T) -> T>,
                    ));
                } else {
                    function_round.push((
                        coset.fft(coeff),
                        Box::new(move |v: T, _x: T, _c: T| v) as Box<dyn Fn(T, T, T) -> T>,
                    ));
                }
            }
            coset = coset.pow(2);
            res.push(function_round);
        }
        res
    }

    pub fn new(
        polynomial: &MultilinearPolynomial<T>,
        interpolate_coset: &Coset<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
        folding_parameter: &Vec<Vec<T>>,
    ) -> Self {
        assert_eq!(polynomial.variable_num(), N);
        let functions = Self::batch_folding(&polynomial, folding_parameter);
        assert_eq!(functions.len(), polynomial.variable_num() + 1);
        // let beta = oracle.borrow_mut().generate_beta();

        // let tuples = Self::calculate_tuples(&functions, beta);
        // assert_eq!(tuples.len(), polynomial.variable_num());
        // let _last = tuples.last().unwrap();
        let evaluations = functions
            .last()
            .unwrap()
            .iter()
            .map(|x| x.coefficients()[0])
            .collect();
        let functions = Self::calculate_functions(&functions, &interpolate_coset);
        Dealer {
            // tuples,
            evaluations,
            prover: One2ManyProver::new(interpolate_coset, functions, oracle),
        }
    }

    pub fn commit_functions(&self, avss_party: &Vec<AvssParty<T, N>>) {
        let verifiers = avss_party.iter().map(|x| x.verifier.clone()).collect();
        self.prover.commit_functions(&verifiers);
    }

    pub fn commit_foldings(&self, avss_party: &Vec<AvssParty<T, N>>) {
        let verifiers = avss_party.iter().map(|x| x.verifier.clone()).collect();
        self.prover.commit_foldings(&verifiers);
    }

    pub fn send_evaluations(&self, avss_party: &mut Vec<AvssParty<T, N>>) {
        for i in 0..avss_party.len() {
            avss_party[i].set_share(self.evaluations[i]);
        }
    }
    // pub fn send_tuples(&self, avss_party: &mut Vec<AvssParty<T, N>>) {
    //     for i in &self.tuples {
    //         for (number, party) in avss_party.iter_mut().enumerate() {
    //             party.add_tuple(&i[number % i.len()]);
    //         }
    //     }
    // }

    pub fn prove(&mut self) {
        self.prover.prove();
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        self.prover.query()
    }
}
