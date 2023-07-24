pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use std::mem::size_of;

    use crate::{prover::One2ManyProver, verifier::One2ManyVerifier};
    use util::{
        algebra::{
            coset::Coset,
            field::{mersenne61_ext::Mersenne61Ext, Field},
            polynomial::MultilinearPolynomial,
        },
        merkle_tree::MERKLE_ROOT_SIZE,
        random_oracle::RandomOracle,
    };

    const SECURITY_BITS: usize = 100;

    fn output_proof_size(variable_num: usize, code_rate: usize, terminate_round: usize) -> usize {
        let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
        let mut interpolate_cosets = vec![Coset::new(
            1 << (variable_num + code_rate),
            Mersenne61Ext::random_element(),
        )];
        for i in 1..variable_num {
            interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
        }
        let oracle = RandomOracle::new(variable_num, SECURITY_BITS / code_rate);
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
        verifier.verify(&folding_proof, &function_proof);
        folding_proof.iter().map(|x| x.proof_size()).sum::<usize>()
            + function_proof.iter().map(|x| x.proof_size()).sum::<usize>()
            + (variable_num - terminate_round) * MERKLE_ROOT_SIZE * 2
            + ((1 << terminate_round) + 1) * size_of::<Mersenne61Ext>() * 2
    }

    #[test]
    fn test_proof_size() {
        for i in 5..21 {
            let proof_size = output_proof_size(i, 4, if i < 9 { i - 4 } else { 5 });
            println!(
                "frolling pcs proof size of {} variables is {} bytes",
                i, proof_size
            );
        }
    }
}
