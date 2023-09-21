pub mod prover;
pub mod verifier;

use util::algebra::field::Field;

#[derive(Debug, Clone)]
pub struct Tuple<T: Field> {
    pub a: T,
    pub b: T,
    pub c: T,
}

impl<T: Field> Tuple<T> {
    pub fn verify(&self, beta: T, folding_param: T) -> bool {
        let v = beta * (self.a + self.b) + folding_param * (self.a - self.b);
        v * T::INVERSE_2 == self.c * beta
    }
}

#[cfg(test)]
mod tests {
    use std::mem::size_of;

    use crate::{prover::FriProver, verifier::FriVerifier, Tuple};
    use util::{
        algebra::{
            coset::Coset, field::mersenne61_ext::Mersenne61Ext, field::Field,
            polynomial::MultilinearPolynomial,
        },
        merkle_tree::MERKLE_ROOT_SIZE,
        random_oracle::RandomOracle,
        CODE_RATE, SECURITY_BITS,
    };

    fn output_proof_size(variable_num: usize) -> usize {
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
        verifier.set_tuples(&tuples);
        assert!(verifier.verify(&folding_proofs, &function_proofs));
        tuples.len() * size_of::<Tuple<Mersenne61Ext>>()
            + variable_num * MERKLE_ROOT_SIZE * 2
            + 2 * size_of::<Mersenne61Ext>()
            + folding_proofs.iter().map(|x| x.proof_size()).sum::<usize>()
            + function_proofs
                .iter()
                .map(|x| x.proof_size())
                .sum::<usize>()
    }

    #[test]
    fn test_proof_size() {
        for i in 5..20 {
            let proof_size = output_proof_size(i);
            println!(
                "gemini proof size of {} variables is {} bytes",
                i, proof_size
            );
        }
    }
}
