use std::collections::HashMap;

use crate::{
    algebra::field::{as_bytes_vec, Field},
    merkle_tree::MerkleTreeVerifier,
};

mod prover;
mod verifier;

pub struct QueryResult<T: Field> {
    proof_bytes: Vec<u8>,
    proof_values: HashMap<usize, T>,
}

impl<T: Field> QueryResult<T> {
    fn verify_merkle_tree(
        &self,
        leaf_indices: &Vec<usize>,
        merkle_verifier: &MerkleTreeVerifier,
        half_len: Option<usize>,
    ) -> bool {
        let leaves: Vec<Vec<u8>> = if let Some(len) = half_len {
            leaf_indices
                .iter()
                .map(|x| {
                    as_bytes_vec(&[
                        self.proof_values.get(x).unwrap().clone(),
                        self.proof_values.get(&(x + len)).unwrap().clone(),
                    ])
                })
                .collect()
        } else {
            leaf_indices
                .iter()
                .map(|x| as_bytes_vec(&[self.proof_values.get(x).unwrap().clone()]))
                .collect()
        };
        let res = merkle_verifier.verify(self.proof_bytes.clone(), leaf_indices, &leaves);
        assert!(res);
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::{
        coset::Coset, field::mersenne61_ext::Mersenne61Ext, polynomial::Polynomial,
    };
    use prover::*;
    use rand::Rng;
    use std::{cell::RefCell, rc::Rc};
    use verifier::*;

    #[test]
    fn rolling_fri_test() {
        let shift = Mersenne61Ext::random_element();
        let domain = Coset::new(1 << 10, shift);
        let poly_degree_bound = 1 << 8;
        let mut functions = vec![];
        let mut function_maps: Vec<fn(Vec<Mersenne61Ext>) -> Mersenne61Ext> = vec![];
        let mut shift = domain.shift();
        let domain_size = domain.size();
        for i in 0..8 {
            let coset = Coset::new(domain_size >> i, shift);
            let poly_eval: Vec<_> = (0..2)
                .into_iter()
                .map(|_| {
                    let poly = Polynomial::random_polynomial(poly_degree_bound >> i);
                    coset.fft(poly.coefficients())
                })
                .collect();
            assert_eq!(poly_eval.len(), 2);
            functions.push(poly_eval);
            function_maps
                .push(|x: Vec<Mersenne61Ext>| x[0] + Mersenne61Ext::from_int(19260817) * x[1]);
            shift *= shift;
        }
        let verifier = Rc::new(RefCell::new(RollingFriVerifier::new(&domain, 8)));
        let prover = Rc::new(RefCell::new(RollingFriProver::new(
            8,
            functions,
            function_maps,
            &domain,
        )));
        verifier.borrow_mut().set_prover(&prover);
        prover.borrow_mut().set_verifier(&verifier);

        prover.borrow_mut().commit_functions();
        prover.borrow_mut().prove();
        let mut points = vec![];
        for _i in 0..1 {
            points.push(rand::thread_rng().gen_range(0..domain.size()));
        }
        let (folding_query, function_query) = prover.borrow().query(&points);
        assert!(verifier
            .borrow()
            .verify(points, folding_query, function_query));
    }
}
