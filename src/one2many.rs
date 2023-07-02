use crate::algebra::field::Field;

mod prover;
mod verifier;

pub struct RandomOracle<T: Field> {
    folding_challenges: Vec<T>,
    usize_elements: Vec<usize>,
}

impl<T: Field> RandomOracle<T> {
    fn new() -> Self {
        RandomOracle {
            folding_challenges: vec![],
            usize_elements: vec![],
        }
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
    use rand::Rng;
    use std::{cell::RefCell, rc::Rc};

    use super::{prover::One2ManyProver, verifier::One2ManyVerifier, *};

    #[test]
    fn test_one_to_many_rolling_fri() {
        let shift = Mersenne61Ext::random_element();
        let interpolate_coset = Coset::new(1 << 11, shift);
        let mut functions_value = vec![];
        let size_each_round = vec![1, 8, 16, 32, 64, 512, 1024, 2048];
        let mut domain_size = interpolate_coset.size();
        let mut shift = interpolate_coset.shift();
        for size in size_each_round {
            let current_domain = Coset::new(domain_size, shift);
            let mut polies: Vec<(
                Vec<Mersenne61Ext>,
                fn(Mersenne61Ext, Mersenne61Ext, Mersenne61Ext) -> Mersenne61Ext,
            )> = vec![];
            for _j in 0..size {
                polies.push((
                    current_domain
                        .fft(Polynomial::random_polynomial(domain_size / 8).coefficients()),
                    |v: Mersenne61Ext, x: Mersenne61Ext, poly: Mersenne61Ext| v,
                ))
            }
            functions_value.push(polies);
            domain_size >>= 1;
            shift *= shift;
        }
        let oracle = Rc::new(RefCell::new(RandomOracle::new()));
        let mut verifiers: Vec<One2ManyVerifier<Mersenne61Ext, 8>> = (0..2048)
            .into_iter()
            .map(|_| One2ManyVerifier::new(&interpolate_coset, &oracle))
            .collect();
        let mut prover: One2ManyProver<Mersenne61Ext, 8> =
            One2ManyProver::new(&interpolate_coset, functions_value, &oracle);
        prover.commit_functions(&mut verifiers);
        prover.prove();
        prover.commit_foldings(&mut verifiers);
        let mut points = vec![];
        for _i in 0..10 {
            points.push(rand::thread_rng().gen_range(0..interpolate_coset.size()));
        }
        let (folding, function) = prover.query(&points);
        let mut folding715 = vec![];
        let mut function715 = vec![];
        for i in 0..8 {
            if i < 7 {
                folding715.push(folding[i][715 % folding[i].len()].clone());
            }
            function715.push(function[i][715 % function[i].len()].clone());
        }
        assert!(verifiers[715].verify(points, folding715, function715));
    }
}
