use super::verifier::RollingFriVerifier;
use super::QueryResult;
use crate::algebra::{
    coset::Coset,
    field::{as_bytes_vec, Field},
};
use crate::merkle_tree::MerkleTreeProver;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

trait InterpolateValue<T: Field> {
    fn get_value(&self, index: usize) -> T;
}

struct FunctionValues<T: Field> {
    functions: Vec<Vec<T>>,
    map: fn(Vec<T>) -> T,
}

impl<T: Field> FunctionValues<T> {
    fn get_value(&self, index: usize) -> T {
        let v = self.functions.iter().map(|x| x[index]).collect();
        let map = self.map;
        let res = map(v);
        res
    }

    fn to_interpolate_values(&self) -> Vec<T> {
        (0..self.leave_number())
            .map(|i| {
                let v = self.functions.iter().map(|x| x[i]).collect();
                let map = self.map;
                map(v)
            })
            .collect()
    }

    fn commit(&self, half: bool) -> Vec<MerkleTreeProver> {
        self.functions
            .iter()
            .map(|x| {
                MerkleTreeProver::new(if half {
                    let len = x.len() / 2;
                    (0..len)
                        .map(|i| as_bytes_vec(&[x[i], x[i + len]]))
                        .collect()
                } else {
                    (0..x.len()).map(|i| as_bytes_vec(&[x[i]])).collect()
                })
            })
            .collect()
    }

    fn leave_number(&self) -> usize {
        self.functions[0].len()
    }
}

pub struct RollingFriProver<T: Field> {
    total_round: usize,
    coset: Coset<T>,
    functions: Vec<FunctionValues<T>>,
    folding_values: Vec<Vec<T>>,
    functions_tree: Vec<Vec<MerkleTreeProver>>,
    folding_tree: Vec<MerkleTreeProver>,
    verifier: Option<Rc<RefCell<RollingFriVerifier<T>>>>,
}

impl<T: Field> RollingFriProver<T> {
    pub fn new(
        total_round: usize,
        function_values: Vec<Vec<Vec<T>>>,
        function_maps: Vec<fn(Vec<T>) -> T>,
        coset: &Coset<T>,
    ) -> RollingFriProver<T> {
        assert_eq!(function_values.len(), function_maps.len());
        RollingFriProver {
            total_round,
            coset: coset.clone(),
            functions: function_values
                .into_iter()
                .zip(function_maps.into_iter())
                .map(|(x, y)| FunctionValues {
                    functions: x,
                    map: y,
                })
                .collect(),
            folding_values: vec![vec![]],
            functions_tree: vec![],
            folding_tree: vec![MerkleTreeProver::new(vec![])],
            verifier: None,
        }
    }

    pub fn set_verifier(&mut self, verifier: &Rc<RefCell<RollingFriVerifier<T>>>) {
        self.verifier = Some(verifier.clone());
    }

    pub fn commit_functions(&mut self) {
        let verifier = self.verifier.clone().unwrap();
        for (idx, fv) in self.functions.iter().enumerate() {
            let merkle_prover = fv.commit(idx == 0);
            let commit = merkle_prover.iter().map(|x| x.commit()).collect();

            verifier.borrow_mut().set_function_root(
                if idx == 0 {
                    fv.leave_number() / 2
                } else {
                    fv.leave_number()
                },
                commit,
                fv.map,
            );
            self.functions_tree.push(merkle_prover);
        }
    }

    fn evaluation_next_domain(
        &self,
        round: usize,
        current_domain: &Coset<T>,
        challenge: T,
    ) -> Vec<T> {
        let mut res = vec![];
        let len = if round == 0 {
            self.functions[round].leave_number()
        } else {
            self.folding_values[round].len()
        };
        let get_folding_value = |i: usize| {
            if round == 0 {
                let fv = &self.functions[round];
                let v = fv.functions.iter().map(|x| x[i]).collect();
                let map = fv.map;
                map(v)
            } else {
                self.folding_values[round][i]
            }
        };
        assert_eq!(len, current_domain.size());
        let inv_2 = T::from_int(2).inverse();
        let mut shift_inv = current_domain.shift().inverse();
        let generator_inv = current_domain.generator().inverse();
        for i in 0..(len / 2) {
            let x = get_folding_value(i);
            let nx = get_folding_value(i + len / 2);
            let new_v = ((x + nx) + challenge * (x - nx) * shift_inv) * inv_2;
            if round < self.total_round - 1 {
                res.push(new_v + challenge.pow(2) * self.functions[round + 1].get_value(i));
            } else {
                res.push(new_v);
            }
            shift_inv *= generator_inv;
        }
        res
    }

    fn merkle_tree_commit(value: &Vec<T>) -> MerkleTreeProver {
        let mut leaf_values = vec![];
        for i in 0..(value.len() / 2) {
            leaf_values.push(as_bytes_vec(&[value[i], value[i + value.len() / 2]]));
        }
        MerkleTreeProver::new(leaf_values)
    }

    pub fn prove(&mut self) {
        let mut domain_size = self.coset.size();
        let mut domain = self.coset.clone();
        let mut shift = domain.shift();
        let verifier = self.verifier.clone().unwrap();
        for i in 0..self.total_round {
            let challenge = verifier.borrow_mut().get_challenge();
            let next_evalutation = self.evaluation_next_domain(i, &domain, challenge);
            self.folding_values.push(next_evalutation);

            shift *= shift;
            domain_size >>= 1;
            domain = Coset::new(domain_size, shift);

            if i < self.total_round - 1 {
                let merkle_tree_prover =
                    Self::merkle_tree_commit(self.folding_values.last().unwrap());
                let commit = merkle_tree_prover.commit();
                verifier.borrow_mut().receive_root(domain_size / 2, &commit);
                self.folding_tree.push(merkle_tree_prover);
            }
        }

        verifier
            .borrow_mut()
            .set_final_value(self.folding_values.last().unwrap()[0]);
    }

    pub fn query(&self, points: &Vec<usize>) -> (Vec<QueryResult<T>>, Vec<Vec<QueryResult<T>>>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = points.clone();

        for i in 0..self.total_round {
            let len = self.functions[i].leave_number();

            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                let query_result = self.functions[i]
                    .functions
                    .iter()
                    .zip(self.functions_tree[i].iter())
                    .map(|(v, p)| {
                        let mut proof_values = HashMap::new();
                        for k in &leaf_indices {
                            proof_values.insert(*k, v[*k]);
                            proof_values.insert(*k + len / 2, v[*k + len / 2]);
                        }
                        let proof_bytes = p.open(&leaf_indices);
                        QueryResult {
                            proof_bytes,
                            proof_values,
                        }
                    })
                    .collect();
                functions_res.push(query_result);
            }

            if i < self.total_round - 1 {
                let query_result = self.functions[i + 1]
                    .functions
                    .iter()
                    .zip(self.functions_tree[i + 1].iter())
                    .map(|(v, p)| {
                        let mut proof_values = HashMap::new();
                        for k in &leaf_indices {
                            proof_values.insert(*k, v[*k]);
                        }
                        let proof_bytes = p.open(&leaf_indices);
                        QueryResult {
                            proof_bytes,
                            proof_values,
                        }
                    })
                    .collect();
                functions_res.push(query_result);
            }

            if i > 0 {
                let mut values = HashMap::new();
                for j in &leaf_indices {
                    values.insert(*j, self.folding_values[i][*j]);
                    values.insert(j + len / 2, self.folding_values[i][*j + len / 2]);
                }
                let proof_bytes = self.folding_tree[i].open(&leaf_indices);
                folding_res.push(QueryResult {
                    proof_bytes,
                    proof_values: values,
                });
            }
        }
        (folding_res, functions_res)
    }
}
