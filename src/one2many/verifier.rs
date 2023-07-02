use crate::rolling_fri::QueryResult;
use crate::{
    algebra::{coset::Coset, field::Field},
    merkle_tree::MerkleTreeVerifier,
};
use std::{cell::RefCell, collections::HashMap, rc::Rc};

use super::RandomOracle;

pub struct One2ManyVerifier<T: Field, const N: usize> {
    coset: Coset<T>,
    function_root: Vec<MerkleTreeVerifier>,
    function_maps: Vec<fn(values: T, x: T, challenge: T) -> T>,
    folding_root: Vec<MerkleTreeVerifier>,
    oracle: Rc<RefCell<RandomOracle<T>>>,
    final_value: Option<T>,
}

impl<T: Field, const N: usize> One2ManyVerifier<T, N> {
    pub fn new(coset: &Coset<T>, oracle: &Rc<RefCell<RandomOracle<T>>>) -> Self {
        One2ManyVerifier {
            coset: coset.clone(),
            function_root: vec![],
            function_maps: vec![],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn set_function(
        &mut self,
        leave_number: usize,
        function_root: &[u8; 32],
        function_map: fn(T, T, T) -> T,
    ) {
        self.function_root.push(MerkleTreeVerifier {
            merkle_root: function_root.clone(),
            leave_number,
        });
        self.function_maps.push(function_map);
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
        mut leaf_indices: Vec<usize>,
        mut folding_proofs: Vec<QueryResult<T>>,
        mut function_proofs: Vec<QueryResult<T>>,
    ) -> bool {
        let mut shift_inv = self.coset.shift().inverse();
        let mut generator_inv = self.coset.generator().inverse();
        let mut domain_size = self.coset.size();
        for i in 0..N {
            leaf_indices = leaf_indices
                .iter_mut()
                .map(|v| *v % (domain_size >> 1))
                .collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let folding_values = if i == 0 {
                let query_results = function_proofs.remove(0);
                query_results.verify_merkle_tree(
                    &leaf_indices,
                    &self.function_root[0],
                    Some(domain_size / 2),
                );
                query_results.proof_values
            } else {
                let query_result = folding_proofs.remove(0);
                query_result.verify_merkle_tree(
                    &leaf_indices,
                    &self.folding_root[i - 1],
                    Some(domain_size / 2),
                );
                let mut interpolation_value = HashMap::new();
                for j in &leaf_indices {
                    interpolation_value
                        .insert(*j, query_result.proof_values.get(j).unwrap().clone());
                    interpolation_value.insert(
                        j + domain_size / 2,
                        query_result
                            .proof_values
                            .get(&(j + domain_size / 2))
                            .unwrap()
                            .clone(),
                    );
                }
                interpolation_value
            };

            if i < N - 1 {
                let function_query_result = function_proofs.remove(0);
                function_query_result.verify_merkle_tree(
                    &leaf_indices,
                    &self.function_root[i + 1],
                    None,
                );
                for j in &leaf_indices {
                    let x = folding_values.get(j).unwrap().clone();
                    let nx = folding_values.get(&(j + domain_size / 2)).unwrap().clone();
                    let v = x
                        + nx
                        + self.oracle.borrow().get_challenge(i)
                            * (x - nx)
                            * shift_inv
                            * generator_inv.pow(*j as u64);
                    let v = v * T::from_int(2).inverse();
                    let function_value = self.function_maps[i + 1](
                        function_query_result.proof_values.get(j).unwrap().clone(),
                        T::random_element(),
                        T::random_element(),
                    );
                    let v = v + self.oracle.borrow().get_challenge(i).pow(2) * function_value;
                    if v != *folding_proofs[0]
                        .proof_values
                        .get(&j)
                        .expect("query missing")
                    {
                        return false;
                    }
                }
            } else {
                for j in &leaf_indices {
                    let x = folding_values.get(j).unwrap().clone();
                    let nx = folding_values.get(&(j + domain_size / 2)).unwrap().clone();
                    let v = x
                        + nx
                        + self.oracle.borrow().get_challenge(i)
                            * (x - nx)
                            * shift_inv
                            * generator_inv.pow(*j as u64);
                    let v = v * T::from_int(2).inverse();
                    if v != self.final_value.unwrap() {
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
