#[macro_use]
extern crate criterion;

use criterion::*;

use std::{cell::RefCell, rc::Rc};

use frolling::algebra::field::Field;
use frolling::algebra::polynomial::MultilinearPolynomial;
use frolling::avss::dealer::Dealer;
use frolling::avss::party::AvssParty;
use frolling::random_oracle::RandomOracle;

use frolling::algebra::coset::Coset;
use frolling::algebra::field::mersenne61_ext::Mersenne61Ext;
use frolling::util::split_n;

fn deal(log_t: usize, code_rate: usize) {
    let oracle = Rc::new(RefCell::new(RandomOracle::new()));
    let interpolate_coset = Coset::new(1 << (log_t + code_rate), Mersenne61Ext::random_element());
    let polynomial = MultilinearPolynomial::random_polynomial(log_t);

    let x_shift = Mersenne61Ext::random_element();
    let coset_x = Coset::new(1 << (log_t + 2), x_shift);
    let mut folding_parameter = vec![];
    let v = split_n((1 << log_t) - 1);
    assert_eq!(v.len(), log_t);
    for i in &v {
        folding_parameter.push(coset_x.pow(*i).all_elements());
    }
    let mut parties = vec![];
    for i in 0..(1 << (log_t + 2)) {
        let mut open_point = vec![];
        for j in 0..log_t {
            open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
        }
        parties.push(AvssParty::new(&interpolate_coset, open_point, &oracle));
    }
    let mut dealer = Dealer::new(&polynomial, &interpolate_coset, &oracle, &folding_parameter);
    dealer.send_evaluations(&mut parties);
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    oracle.borrow_mut().generate_queries(20);
    let (folding, function) = dealer.query();
    let mut folding715 = vec![];
    let mut function715 = vec![];
    for i in 0..log_t {
        if i < log_t - 1 {
            folding715.push(folding[i][715 % folding[i].len()].clone());
        }
        function715.push(function[i][715 % function[i].len()].clone());
    }
    assert!(parties[715].verify(folding715, function715));
    assert_eq!(
        parties[715].share(),
        polynomial.evaluate(parties[715].open_point())
    );
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("create proof");
    group.sample_size(10);

    group.bench_function("with_setup", move |b| {
        b.iter(|| {
            deal(12, 3);
        })
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
