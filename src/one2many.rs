use rand::Rng;

use crate::algebra::field::Field;

mod prover;
mod verifier;

pub struct RandomOracle<T: Field> {
    folding_challenges: Vec<T>,
    usize_elements: Option<Vec<usize>>,
}

impl<T: Field> RandomOracle<T> {
    fn new() -> Self {
        RandomOracle {
            folding_challenges: vec![],
            usize_elements: None,
        }
    }

    fn generate_queries(&mut self, len: usize) {
        self.usize_elements = Some(
            (0..len)
                .into_iter()
                .map(|_| rand::thread_rng().gen())
                .collect(),
        )
    }

    fn get_challenge(&self, index: usize) -> T {
        self.folding_challenges[index]
    }

    fn generate_challenge(&mut self) -> T {
        let challenge = T::random_element();
        self.folding_challenges.push(challenge);
        challenge
    }
}

#[cfg(test)]
mod tests {
    use crate::algebra::{
        coset::Coset, field::mersenne61_ext::Mersenne61Ext, polynomial::Polynomial,
    };
    use std::{cell::RefCell, rc::Rc};

    use super::{prover::One2ManyProver, verifier::One2ManyVerifier, *};

    #[test]
    fn test_one_to_many_rolling_fri() {
        let shift = Mersenne61Ext::random_element();
        let interpolate_coset = Coset::new(1 << 11, shift);
        let mut functions_value = vec![];
        let size_each_round = vec![8, 16, 32, 64, 512, 1024, 2048];
        let mut domain_size = interpolate_coset.size() / 2;
        let mut shift = interpolate_coset.shift().pow(2);
        for size in size_each_round {
            let current_domain = Coset::new(domain_size, shift);
            let polies = (0..size)
                .into_iter()
                .map(|_| {
                    (
                        current_domain
                            .fft(Polynomial::random_polynomial(domain_size / 8).coefficients()),
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
            functions_value.push(polies);
            domain_size >>= 1;
            shift *= shift;
        }
        let oracle = Rc::new(RefCell::new(RandomOracle::new()));
        let mut verifiers: Vec<One2ManyVerifier<Mersenne61Ext, 8>> = (0..4096)
            .into_iter()
            .map(|_| One2ManyVerifier::new(&interpolate_coset, &oracle))
            .collect();
        let commited_polynomial = Polynomial::random_polynomial(1 << 8);
        let mut prover: One2ManyProver<Mersenne61Ext, 8> = One2ManyProver::new(
            &interpolate_coset,
            &commited_polynomial,
            functions_value,
            &oracle,
        );
        prover.commit_functions(&mut verifiers);
        prover.prove();
        prover.commit_foldings(&mut verifiers);
        oracle.borrow_mut().generate_queries(10);
        let (committed, folding, function) = prover.query();
        let mut folding715 = vec![];
        let mut function715 = vec![];
        for i in 0..7 {
            folding715.push(folding[i][715 % folding[i].len()].clone());
            function715.push(function[i][715 % function[i].len()].clone());
        }
        assert!(verifiers[715].verify(committed, folding715, function715));
    }
}
