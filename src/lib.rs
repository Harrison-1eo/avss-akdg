pub mod algebra {
    pub mod coset;
    pub mod field;
    pub mod polynomial;
}
pub mod merkle_tree;
pub mod one2many;
pub mod util;
pub mod avss {
    pub mod dealer;
    pub mod party;
}
pub mod random_oracle;

const TERMINATE_ROUND: usize = 5;

use std::{cell::RefCell, rc::Rc};

use crate::algebra::field::Field;
use crate::algebra::polynomial::MultilinearPolynomial;
use crate::avss::dealer::Dealer;
use crate::avss::party::AvssParty;
use crate::random_oracle::RandomOracle;

use crate::algebra::coset::Coset;
use crate::algebra::field::mersenne61_ext::Mersenne61Ext;
use crate::util::split_n;

pub fn deal(log_n: usize, code_rate: usize) {
    let log_t = log_n - 2;
    let oracle = Rc::new(RefCell::new(RandomOracle::new()));
    let mut interpolate_cosets = vec![Coset::new(
        1 << (log_t + code_rate),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..log_t {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let polynomial = MultilinearPolynomial::random_polynomial(log_t);

    let x_shift = Mersenne61Ext::random_element();
    let coset_x = Coset::new(1 << log_n, x_shift);
    let mut folding_parameter = vec![];
    let v = split_n((1 << log_t) - 1);
    for i in &v {
        folding_parameter.push(coset_x.pow(*i).all_elements());
    }
    let mut parties = vec![];
    for i in 0..(1 << log_n) {
        let mut open_point = vec![];
        for j in 0..log_t {
            open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
        }
        parties.push(AvssParty::new(
            log_t - TERMINATE_ROUND,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_t - TERMINATE_ROUND,
        &polynomial,
        &interpolate_cosets,
        &oracle,
        &folding_parameter,
    );
    dealer.send_evaluations(&mut parties);
    assert_eq!(
        parties[715].share(),
        polynomial.evaluate(parties[715].open_point())
    );
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    oracle.borrow_mut().generate_queries(20);
    let (folding, function) = dealer.query();
    let mut folding715 = vec![];
    let mut function715 = vec![];
    for i in 0..(log_t - TERMINATE_ROUND) {
        if i < log_t - TERMINATE_ROUND - 1 {
            folding715.push(folding[i][715 % folding[i].len()].clone());
        }
        function715.push(function[i][715 % function[i].len()].clone());
    }
    assert!(parties[715].verify(folding715, function715));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_avss() {
        deal(12, 3);
    }
}
