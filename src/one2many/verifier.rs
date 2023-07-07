use crate::rolling_fri::QueryResult;
use crate::{
    algebra::{coset::Coset, field::Field},
    merkle_tree::MerkleTreeVerifier,
};
use std::{cell::RefCell, rc::Rc};

use super::RandomOracle;

pub struct One2ManyVerifier<T: Field, const N: usize> {
    interpolate_cosets: Vec<Coset<T>>,
    // committed_polynomial: Option<MerkleTreeVerifier>,
    function_root: Vec<MerkleTreeVerifier>,
    function_maps: Vec<Box<dyn Fn(T, T, T) -> T>>,
    folding_root: Vec<MerkleTreeVerifier>,
    oracle: Rc<RefCell<RandomOracle<T>>>,
    final_value: Option<T>,
}

impl<T: Field, const N: usize> One2ManyVerifier<T, N> {
    pub fn new(coset: &Coset<T>, oracle: &Rc<RefCell<RandomOracle<T>>>) -> Self {
        let mut cosets = vec![coset.clone()];
        for _ in 1..N {
            cosets.push(cosets.last().as_ref().unwrap().pow(2));
        }
        One2ManyVerifier {
            interpolate_cosets: cosets,
            function_root: vec![],
            function_maps: vec![],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    // pub fn commit_polynomial(&mut self, leave_number: usize, function_root: &[u8; 32]) {
    //     self.committed_polynomial = Some(MerkleTreeVerifier {
    //         merkle_root: function_root.clone(),
    //         leave_number,
    //     });
    // }

    pub fn set_map(&mut self, function_map: Box<dyn Fn(T, T, T) -> T>) {
        self.function_maps.push(function_map);
    }

    pub fn set_function(&mut self, leave_number: usize, function_root: &[u8; 32]) {
        self.function_root.push(MerkleTreeVerifier {
            merkle_root: function_root.clone(),
            leave_number,
        });
    }

    pub fn receive_folding_root(&mut self, leave_number: usize, folding_root: [u8; 32]) {
        self.folding_root.push(MerkleTreeVerifier {
            leave_number,
            merkle_root: folding_root,
        });
    }

    pub fn set_final_value(&mut self, value: T) {
        self.final_value = Some(value);
    }

    pub fn verify(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
    ) -> bool {
        let mut leaf_indices = self.oracle.borrow().usize_elements.clone().unwrap();
        let mut shift_inv = self.interpolate_cosets[0].shift().inverse();
        let mut generator_inv = self.interpolate_cosets[0].generator().inverse();
        let mut domain_size = self.interpolate_cosets[0].size();
        for i in 0..N {
            leaf_indices = leaf_indices
                .iter_mut()
                .map(|v| *v % (domain_size >> 1))
                .collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                function_proofs[i].verify_merkle_tree(
                    &leaf_indices,
                    &self.function_root[0],
                    Some(domain_size / 2),
                );
            } else {
                folding_proofs[i - 1].verify_merkle_tree(
                    &leaf_indices,
                    &self.folding_root[i - 1],
                    Some(domain_size / 2),
                );
            }

            let challenge = self.oracle.borrow().get_challenge(i);
            let get_folding_value = |index: &usize| {
                if i == 0 {
                    self.function_maps[i](
                        function_proofs[i].proof_values[index],
                        self.interpolate_cosets[i].all_elements()[*index],
                        challenge,
                    )
                } else {
                    folding_proofs[i - 1].proof_values[index]
                }
            };

            let function_values = if i < N - 1 {
                let function_query_result = &function_proofs[i + 1];
                function_query_result.verify_merkle_tree(
                    &leaf_indices,
                    &self.function_root[i + 1],
                    None,
                );
                Some(&function_query_result.proof_values)
            } else {
                None
            };
            for j in &leaf_indices {
                let x = get_folding_value(j);
                let nx = get_folding_value(&(j + domain_size / 2));
                let v = x + nx + challenge * (x - nx) * shift_inv * generator_inv.pow(*j);
                let v = v * T::from_int(2).inverse();
                if i == N - 1 {
                    if v != self.final_value.unwrap() {
                        return false;
                    }
                } else {
                    let function_value = self.function_maps[i + 1](
                        function_values.as_ref().unwrap()[j],
                        self.interpolate_cosets[i + 1].all_elements()[*j],
                        challenge,
                    );
                    let v = v + challenge.pow(2) * function_value;
                    if i < N - 1 && v != folding_proofs[i].proof_values[j] {
                        return false;
                    }
                }
            }

            shift_inv *= shift_inv;
            generator_inv *= generator_inv;
            domain_size >>= 1;
        }
        true
    }
}
