extern crate criterion;
use criterion::*;

use pcs::{prover::One2ManyProver, verifier::One2ManyVerifier};
use util::{
    algebra::{
        coset::Coset,
        field::{mersenne61_ext::Mersenne61Ext, Field},
        polynomial::MultilinearPolynomial,
    },
    random_oracle::RandomOracle,
};

fn commit(
    criterion: &mut Criterion,
    variable_num: usize,
    code_rate: usize,
    terminate_round: usize,
) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + code_rate),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, 33);

    criterion.bench_function(&format!("commit {}", variable_num), move |b| {
        b.iter_batched(
            || polynomial.clone(),
            |p| {
                let prover = One2ManyProver::new(
                    variable_num - terminate_round,
                    &interpolate_cosets,
                    p,
                    &oracle,
                );
                prover.commit_polynomial();
            },
            BatchSize::SmallInput,
        )
    });
}

fn bench_commit(c: &mut Criterion) {
    for i in 5..25 {
        commit(c, i, 4, if i < 9 { i - 4 } else { 5 });
    }
}

fn open(criterion: &mut Criterion, variable_num: usize, code_rate: usize, terminate_round: usize) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + code_rate),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, 33);
    let prover = One2ManyProver::new(
        variable_num - terminate_round,
        &interpolate_cosets,
        polynomial,
        &oracle,
    );
    let commit = prover.commit_polynomial();
    let mut verifier = One2ManyVerifier::new(
        variable_num - terminate_round,
        variable_num,
        &interpolate_cosets,
        commit,
        &oracle,
    );
    let open_point = verifier.get_open_point();

    criterion.bench_function(&format!("open {}", variable_num), move |b| {
        b.iter_batched(
            || prover.clone(),
            |mut p| {
                p.commit_functions(&open_point, &mut verifier);
                p.prove();
                p.commit_foldings(&mut verifier);
                p.query();
            },
            BatchSize::SmallInput,
        )
    });
}

fn bench_open(c: &mut Criterion) {
    for i in 5..25 {
        open(c, i, 4, if i < 9 { i - 4 } else { 5 });
    }
}

fn verify(
    criterion: &mut Criterion,
    variable_num: usize,
    code_rate: usize,
    terminate_round: usize,
) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + code_rate),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, 33);
    let mut prover = One2ManyProver::new(
        variable_num - terminate_round,
        &interpolate_cosets,
        polynomial,
        &oracle,
    );
    let commit = prover.commit_polynomial();
    let mut verifier = One2ManyVerifier::new(
        variable_num - terminate_round,
        variable_num,
        &interpolate_cosets,
        commit,
        &oracle,
    );
    let open_point = verifier.get_open_point();

    prover.commit_functions(&open_point, &mut verifier);
    prover.prove();
    prover.commit_foldings(&mut verifier);
    let (folding_proof, function_proof) = prover.query();
    criterion.bench_function(&format!("verify {}", variable_num), move |b| {
        b.iter(|| {
            verifier.verify(&folding_proof, &function_proof);
        })
    });
}

fn bench_verify(c: &mut Criterion) {
    for i in 5..25 {
        verify(c, i, 4, if i < 9 { i - 4 } else { 5 });
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_commit, bench_open, bench_verify
}

criterion_main!(benches);
