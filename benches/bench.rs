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

fn vss_deal(log_n: usize, code_rate: usize, terminate_round: usize) {
    let log_t = log_n - 1;
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
            log_t - terminate_round,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_t - terminate_round,
        &polynomial,
        &interpolate_cosets,
        &oracle,
        &folding_parameter,
    );
    dealer.send_evaluations(&mut parties);
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    oracle.borrow_mut().generate_queries(33);
    dealer.query();
}

fn vss_verify(c: &mut Criterion, log_n: usize, code_rate: usize, terminate_round: usize) {
    let log_t = log_n - 1;
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
            log_t - terminate_round,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_t - terminate_round,
        &polynomial,
        &interpolate_cosets,
        &oracle,
        &folding_parameter,
    );
    dealer.send_evaluations(&mut parties);
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    oracle.borrow_mut().generate_queries(33);
    let (folding, function) = dealer.query();
    let mut folding0 = vec![];
    let mut function0 = vec![];
    for i in 0..(log_t - terminate_round) {
        if i < log_t - terminate_round - 1 {
            folding0.push(folding[i][0].clone());
        }
        function0.push(function[i][0].clone());
    }
    let mut group = c.benchmark_group("verify proof");
    group.sample_size(10);
    group.bench_function(format!("vss verify {}", log_n), move |b| {
        b.iter(|| {
            parties[0].verify(folding0.clone(), function0.clone());
        })
    });
}

fn bench_vss_deal(c: &mut Criterion) {
    let mut group = c.benchmark_group("create proof");
    group.sample_size(10);
    for i in 5..21 {
        let terminate_round = if i < 9 { i - 4 } else { 5 };
        group.bench_function(format!("vss prove {}", i), move |b| {
            b.iter(|| {
                vss_deal(i, 3, terminate_round);
            })
        });
    }
}

fn bench_vss_verify(c: &mut Criterion) {
    for i in 5..21 {
        let terminate_round = if i < 9 { i - 4 } else { 5 };
        vss_verify(c, i, 3, terminate_round);
    }
}

fn avss_deal(log_n: usize, code_rate: usize, terminate_round: usize) {
    let log_t = log_n - 2;
    let log_d = log_t * 2;
    let oracle = Rc::new(RefCell::new(RandomOracle::new()));
    let mut interpolate_cosets = vec![Coset::new(
        1 << (log_t * 2 + code_rate),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..log_d {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let polynomial = MultilinearPolynomial::random_polynomial(log_d);

    let x_shift = Mersenne61Ext::random_element();
    let coset_x = Coset::new(1 << log_n, x_shift);
    let mut folding_parameter = vec![];
    let v = split_n((1 << log_t) - 1);
    for i in &v {
        folding_parameter.push(coset_x.pow(*i).all_elements());
    }
    let y_shift = Mersenne61Ext::random_element();
    let coset_y = Coset::new(1 << log_n, y_shift);
    let last_len = folding_parameter.last().unwrap().len();
    for i in &v {
        folding_parameter.push(
            coset_y
                .pow(*i)
                .all_elements()
                .iter()
                .map(|x| (0..last_len).into_iter().map(|_| *x).collect::<Vec<_>>())
                .flatten()
                .collect(),
        );
    }
    let mut parties = vec![];
    for i in 0..(1 << (log_n * 2)) {
        let mut open_point = vec![];
        for j in 0..log_d {
            open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
        }
        parties.push(AvssParty::new(
            log_d - terminate_round,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_d - terminate_round,
        &polynomial,
        &interpolate_cosets,
        &oracle,
        &folding_parameter,
    );
    dealer.send_evaluations(&mut parties);
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    oracle.borrow_mut().generate_queries(33);
    dealer.query();
}

fn bench_avss_deal(c: &mut Criterion) {
    let mut group = c.benchmark_group("create proof");
    group.sample_size(10);
    for i in 3..11 {
        let terminate_round = if i < 5 { i * 2 - 5 } else { 5 };
        group.bench_function(format!("vss prove {}", i), move |b| {
            b.iter(|| {
                avss_deal(i, 3, terminate_round);
            })
        });
    }
}

criterion_group!(
    benches,
    /*bench_vss_deal, bench_vss_verify,*/ bench_avss_deal
);
criterion_main!(benches);
