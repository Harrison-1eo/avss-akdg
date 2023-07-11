pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use crate::algebra::{
        coset::Coset, field::mersenne61_ext::Mersenne61Ext, field::Field, polynomial::Polynomial,
    };
    use crate::random_oracle::RandomOracle;
    use std::{cell::RefCell, rc::Rc};

    use super::{prover::One2ManyProver, verifier::One2ManyVerifier};

    #[test]
    fn test_one_to_many_rolling_fri() {
        let mut functions_value = vec![];
        let size_each_round = vec![1, 8, 16, 32, 64];
        let mut interpolate_cosets = vec![Coset::new(1 << 11, Mersenne61Ext::random_element())];
        for i in 1..8 {
            interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
        }
        for i in size_each_round.iter().enumerate() {
            let current_domain = &interpolate_cosets[i.0];
            let functions = (0..(*i.1))
                .into_iter()
                .map(|_| {
                    (
                        current_domain.fft(
                            Polynomial::random_polynomial(current_domain.size() / 8).coefficients(),
                        ),
                        Box::new(|v, x, c| v + c * x)
                            as Box<
                                dyn Fn(
                                    Mersenne61Ext,
                                    Mersenne61Ext,
                                    Mersenne61Ext,
                                ) -> Mersenne61Ext,
                            >,
                    )
                })
                .collect();
            functions_value.push(functions);
        }
        let oracle = Rc::new(RefCell::new(RandomOracle::new()));
        let verifiers = (0..4096)
            .into_iter()
            .map(|_| {
                Rc::new(RefCell::new(One2ManyVerifier::new(
                    5,
                    8,
                    &interpolate_cosets,
                    &oracle,
                )))
            })
            .collect::<Vec<_>>();
        verifiers.iter().for_each(|x| {
            for _ in 0..8 {
                x.borrow_mut().set_map(Rc::new(|v, x, c| v + c * x));
            }
        });
        let mut prover = One2ManyProver::new(5, &interpolate_cosets, functions_value, &oracle);
        prover.commit_functions(&verifiers);
        prover.prove();
        prover.commit_foldings(&verifiers);
        oracle.borrow_mut().generate_queries(10);
        let (folding, function) = prover.query();
        let mut folding715 = vec![];
        let mut function715 = vec![];
        for i in 0..5 {
            if i < 4 {
                folding715.push(folding[i][715 % folding[i].len()].clone());
            }
            function715.push(function[i][715 % function[i].len()].clone());
        }
        assert!(verifiers[715].borrow().verify(folding715, function715));
    }
}
