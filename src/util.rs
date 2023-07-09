use crate::algebra::field::{as_bytes_vec, Field};
use crate::merkle_tree::MerkleTreeVerifier;
use std::collections::HashMap;

#[derive(Clone)]
pub struct QueryResult<T: Field> {
    pub proof_bytes: Vec<u8>,
    pub proof_values: HashMap<usize, T>,
}

impl<T: Field> QueryResult<T> {
    pub fn verify_merkle_tree(
        &self,
        leaf_indices: &Vec<usize>,
        merkle_verifier: &MerkleTreeVerifier,
    ) -> bool {
        let leaves: Vec<Vec<u8>> =// if let Some(len) = half_len {
            leaf_indices
                .iter()
                .map(|x| {
                    as_bytes_vec(&[
                        self.proof_values.get(x).unwrap().clone(),
                        self.proof_values.get(&(x + merkle_verifier.leave_number)).unwrap().clone(),
                    ])
                })
                .collect();
        let res = merkle_verifier.verify(self.proof_bytes.clone(), leaf_indices, &leaves);
        assert!(res);
        res
    }
}

pub fn split_n(mut n: usize) -> Vec<usize> {
    let mut res = vec![];
    let mut i = 1;
    while i < n {
        res.push(i);
        n -= i;
        i <<= 1;
    }
    if n > 0 {
        res.push(n);
    }
    res.sort_by(|x, y| y.trailing_zeros().cmp(&x.trailing_zeros()));
    res
}
