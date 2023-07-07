use std::cell::RefCell;
use std::rc::Rc;

use super::party::{AvssParty, Tuple};
use crate::algebra::coset::Coset;
use crate::one2many::RandomOracle;
use crate::rolling_fri::QueryResult;
use crate::{
    algebra::{field::Field, polynomial::MultilinearPolynomial},
    one2many::prover::One2ManyProver,
};

pub struct Dealer<T: Field, const N: usize> {
    prover: One2ManyProver<T, N>,
    tuples: Vec<Vec<Tuple<T>>>,
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

    fn calculate_tuples(
        functions: &Vec<Vec<MultilinearPolynomial<T>>>,
        beta: T,
    ) -> Vec<Vec<Tuple<T>>> {
        let mut tuples = vec![];
        for i in 0..(functions.len() - 1) {
            let mut round_tuple = vec![];
            for j in 0..functions[i + 1].len() {
                let k = j % functions[i].len();
                round_tuple.push(Tuple {
                    a: functions[i][k].evaluate_as_polynomial(beta),
                    b: functions[i][k].evaluate_as_polynomial(-beta),
                    c: functions[i + 1][j].evaluate_as_polynomial(beta * beta),
                });
            }
            tuples.push(round_tuple);
        }
        tuples
    }

    pub fn calculate_functions(
        functions: &Vec<Vec<MultilinearPolynomial<T>>>,
        coset: &Coset<T>,
        tuples: &Vec<Vec<Tuple<T>>>,
        beta: T,
    ) -> Vec<Vec<(Vec<T>, Box<dyn Fn(T, T, T) -> T>)>> {
        let mut coset = coset.clone();
        let mut res = vec![];
        for i in 0..(functions.len() - 1) {
            let mut function_round = vec![];
            for j in 0..functions[i].len() {
                let coeff = functions[i][j].coefficients();
                let a = tuples[i][j].a;
                let b = tuples[i][j].b;
                if i == 0 {
                    function_round.push((
                        coset.fft(coeff),
                        Box::new(move |v: T, x: T, _c: T| {
                            (v - b) * (x + beta).inverse() + v + (v - a) * (x - beta).inverse()
                        }) as Box<dyn Fn(T, T, T) -> T>,
                    ));
                } else {
                    let c = tuples[i - 1][j].c;
                    function_round.push((
                        coset.fft(coeff),
                        Box::new(move |v: T, x: T, _challenge: T| {
                            (v - c) * (x - beta * beta).inverse()
                                + (v - b) * (x + beta).inverse()
                                + v
                                + (v - a) * (x - beta).inverse()
                        }) as Box<dyn Fn(T, T, T) -> T>,
                    ));
                }
            }
            coset = coset.pow(2);
            res.push(function_round);
        }
        res
    }

    pub fn new(
        polynomial: MultilinearPolynomial<T>,
        interpolate_coset: &Coset<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
        folding_parameter: &Vec<Vec<T>>,
    ) -> Self {
        assert_eq!(polynomial.variable_num(), N);
        let functions = Self::batch_folding(&polynomial, folding_parameter);
        assert_eq!(functions.len(), polynomial.variable_num() + 1);
        let beta = oracle.borrow_mut().generate_beta();

        let tuples = Self::calculate_tuples(&functions, beta);
        assert_eq!(tuples.len(), polynomial.variable_num());
        let _last = tuples.last().unwrap();

        let functions = Self::calculate_functions(&functions, &interpolate_coset, &tuples, beta);
        Dealer {
            tuples,
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

    pub fn send_tuples(&self, avss_party: &mut Vec<AvssParty<T, N>>) {
        for i in &self.tuples {
            for (number, party) in avss_party.iter_mut().enumerate() {
                party.add_tuple(&i[number % i.len()]);
            }
        }
    }

    pub fn prove(&mut self) {
        self.prover.prove();
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        self.prover.query()
    }
}
