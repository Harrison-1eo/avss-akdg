use super::prover::RollingFriProver;
use super::QueryResult;
use crate::algebra::{coset::Coset, field::Field};
use crate::merkle_tree::MerkleTreeVerifier;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct RollingFriVerifier<T: Field> {
    total_round: usize,
    coset: Coset<T>,
    function_root: Vec<Vec<MerkleTreeVerifier>>,
    function_maps: Vec<fn(Vec<T>) -> T>,
    challenges: Vec<T>,
    folding_root: Vec<MerkleTreeVerifier>,
    prover: Option<Rc<RefCell<RollingFriProver<T>>>>,
    final_value: Option<T>,
}

impl<T: Field> RollingFriVerifier<T> {
    pub fn new(coset: &Coset<T>, total_round: usize) -> RollingFriVerifier<T> {
        RollingFriVerifier {
            coset: coset.clone(),
            total_round,
            function_root: vec![],
            function_maps: vec![],
            challenges: vec![],
            folding_root: vec![MerkleTreeVerifier {
                merkle_root: [0; 32],
                leave_number: 0,
            }],
            prover: None,
            final_value: None,
        }
    }

    pub fn set_function_root(
        &mut self,
        leave_number: usize,
        function_root: Vec<[u8; 32]>,
        function_map: fn(Vec<T>) -> T,
    ) {
        let v = function_root
            .into_iter()
            .map(|x| MerkleTreeVerifier {
                merkle_root: x,
                leave_number,
            })
            .collect();
        self.function_root.push(v);
        self.function_maps.push(function_map);
    }

    pub fn receive_root(&mut self, leave_number: usize, folding_root: &[u8; 32]) {
        self.folding_root.push(MerkleTreeVerifier {
            leave_number,
            merkle_root: folding_root.clone(),
        });
    }

    pub fn get_challenge(&mut self) -> T {
        let challenge = T::random_element();
        self.challenges.push(challenge);
        challenge
    }

    pub fn set_prover(&mut self, prover: &Rc<RefCell<RollingFriProver<T>>>) {
        self.prover = Some(prover.clone());
    }

    pub fn verify(
        &self,
        mut leaf_indices: Vec<usize>,
        mut folding_proofs: Vec<QueryResult<T>>,
        mut function_proofs: Vec<Vec<QueryResult<T>>>,
    ) -> bool {
        let mut shift_inv = self.coset.shift().inverse();
        let mut generator_inv = self.coset.generator().inverse();
        let mut domain_size = self.coset.size();
        for i in 0..self.total_round {
            leaf_indices = leaf_indices
                .iter_mut()
                .map(|v| *v % (domain_size >> 1))
                .collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let folding_values = if i == 0 {
                let query_results = function_proofs.remove(0);
                for j in 0..query_results.len() {
                    query_results[j].verify_merkle_tree(
                        &leaf_indices,
                        &self.function_root[0][j],
                        Some(domain_size / 2),
                    );
                }
                leaf_indices
                    .iter()
                    .flat_map(|j| {
                        let values: Vec<T> = query_results
                            .iter()
                            .map(|x| x.proof_values.get(j).unwrap().clone())
                            .collect();
                        let value1 = self.function_maps[i](values);
                        let values: Vec<T> = query_results
                            .iter()
                            .map(|x| x.proof_values.get(&(j + domain_size / 2)).unwrap().clone())
                            .collect();
                        let value2 = self.function_maps[i](values);
                        [(*j, value1), (*j + domain_size / 2, value2)]
                    })
                    .collect()
            } else {
                let query_result = folding_proofs.remove(0);
                query_result.verify_merkle_tree(
                    &leaf_indices,
                    &self.folding_root[i],
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

            if i < self.total_round - 1 {
                let function_query_result = function_proofs.remove(0);
                for j in 0..function_query_result.len() {
                    function_query_result[j].verify_merkle_tree(
                        &leaf_indices,
                        &self.function_root[i + 1][j],
                        None,
                    );
                }
                for j in &leaf_indices {
                    let x = folding_values.get(j).unwrap().clone();
                    let nx = folding_values.get(&(j + domain_size / 2)).unwrap().clone();
                    let v =
                        x + nx + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(*j);
                    let v = v * T::from_int(2).inverse();
                    let function_value = self.function_maps[i + 1](
                        function_query_result
                            .iter()
                            .map(|x| x.proof_values.get(j).unwrap().clone())
                            .collect(),
                    );
                    let v = v + self.challenges[i].pow(2) * function_value;
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
                    let v =
                        x + nx + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(*j);
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

    pub fn set_final_value(&mut self, final_value: T) {
        self.final_value = Some(final_value);
    }
}
