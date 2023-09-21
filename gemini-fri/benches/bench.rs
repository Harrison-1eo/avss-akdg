extern crate criterion;

use criterion::*;

use gemini_fri::{prover::FriProver, verifier::FriVerifier};
use util::{
    algebra::{
        coset::Coset, field::mersenne61_ext::Mersenne61Ext, field::Field,
        polynomial::MultilinearPolynomial,
    },
    random_oracle::RandomOracle,
};

use util::{CODE_RATE, SECURITY_BITS};

fn commit(variable_num: usize) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + CODE_RATE),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, SECURITY_BITS / CODE_RATE);
    let prover = FriProver::new(variable_num, &interpolate_cosets, polynomial, &oracle);
    prover.commit_first_polynomial();
}

fn bench_commit(c: &mut Criterion) {
    for i in 5..20 {
        c.bench_function(&format!("bench gemini commit {}", i), |b| {
            b.iter(|| commit(i))
        });
    }
}

fn open(criterion: &mut Criterion, variable_num: usize) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + CODE_RATE),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, SECURITY_BITS / CODE_RATE);
    let prover = FriProver::new(variable_num, &interpolate_cosets, polynomial, &oracle);
    let commitment = prover.commit_first_polynomial();
    let mut verifier = FriVerifier::new(variable_num, &interpolate_cosets, commitment, &oracle);
    let open_point = verifier.get_open_point();
    criterion.bench_function(&format!("gemini prove {}", variable_num), |b| {
        b.iter_batched(
            || (prover.clone(), verifier.clone()),
            |(mut p, mut v)| {
                p.commit_functions(&mut v, &open_point);
                p.compute_tuples();
                p.prove();
                p.commit_foldings(&mut v);
                p.query();
            },
            BatchSize::SmallInput,
        )
    });
}

fn bench_open(c: &mut Criterion) {
    for i in 5..20 {
        open(c, i);
    }
}

fn verify(criterion: &mut Criterion, variable_num: usize) {
    let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (variable_num + CODE_RATE),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..variable_num {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let oracle = RandomOracle::new(variable_num, SECURITY_BITS / CODE_RATE);
    let mut prover = FriProver::new(variable_num, &interpolate_cosets, polynomial, &oracle);
    let commitment = prover.commit_first_polynomial();
    let mut verifier = FriVerifier::new(variable_num, &interpolate_cosets, commitment, &oracle);
    let open_point = verifier.get_open_point();
    prover.commit_functions(&mut verifier, &open_point);
    let tuples = prover.compute_tuples();
    prover.prove();
    prover.commit_foldings(&mut verifier);
    let (folding_proofs, function_proofs) = prover.query();
    criterion.bench_function(&format!("gemini verify {}", variable_num), |b| {
        b.iter_batched(
            || verifier.clone(),
            |mut v| {
                v.set_tuples(&tuples);
                assert!(v.verify(&folding_proofs, &function_proofs));
            },
            BatchSize::SmallInput,
        )
    });
}

fn bench_verify(c: &mut Criterion) {
    for i in 5..20 {
        verify(c, i);
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_commit, bench_open, bench_verify
}

criterion_main!(benches);
