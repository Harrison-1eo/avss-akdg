use crate::algebra::{
    coset::Coset,
    field::{as_bytes_vec, Field},
    polynomial::*,
};
use crate::merkle_tree::{MerkleTreeProver, MerkleTreeVerifier};
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct RollingFriVerifier<T: Field> {
    total_round: usize,
    coset: Coset<T>,
    function_root: Vec<MerkleTreeVerifier>,
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
            challenges: vec![],
            folding_root: vec![MerkleTreeVerifier {merkle_root: [0; 32], leave_number: 0}],
            prover: None,
            final_value: None,
        }
    }

    fn set_function_root(&mut self, leave_number: usize, function_root: &[u8; 32]) {
        self.function_root.push(MerkleTreeVerifier { 
            merkle_root: function_root.clone(), 
            leave_number 
        });
    }

    fn receive_root(&mut self, leave_number: usize, folding_root: &[u8; 32]) {
        self.folding_root.push(MerkleTreeVerifier {
            leave_number,
            merkle_root: folding_root.clone(),
        });
    }

    fn get_challenge(&mut self) -> T {
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
        mut folding_proofs: Vec<(Vec<u8>, HashMap<usize, T>)>,
        mut function_proofs: Vec<(Vec<u8>, HashMap<usize, T>)>
    ) -> bool {
        let mut shift_inv = self.coset.shift().inverse();
        let mut generator_inv = self.coset.generator().inverse();
        let mut domain_size = self.coset.size();
        for i in 0..self.total_round {
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (domain_size >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();
            let (folding_proof_bytes, folding_values) = if i == 0 {
                function_proofs.remove(0)
            } else {
                folding_proofs.remove(0)
            };
            let open_values = leaf_indices
                .iter()
                .map(|v| {
                    as_bytes_vec(&[
                        folding_values.get(v).unwrap().clone(),
                        folding_values.get(&(v + domain_size / 2)).unwrap().clone(),
                    ])
                })
                .collect();
            if i == 0 {
                if !self.function_root[i].verify(folding_proof_bytes, &leaf_indices, &open_values) {
                    return false;
                }
            } else {
                if !self.folding_root[i].verify(folding_proof_bytes, &leaf_indices, &open_values) {
                    return false;
                }
            }
            
            if i < self.total_round - 1 {
                let (function_proof_bytes, function_values) = function_proofs.remove(0);
                let open_values = leaf_indices.iter().map(|v| {
                    as_bytes_vec(&[function_values.get(v).unwrap().clone()])
                })
                .collect();
                if !self.function_root[i + 1].verify(function_proof_bytes, &leaf_indices, &open_values) {
                    return false;
                }
                for j in &leaf_indices {
                    let x = folding_values.get(j).unwrap().clone();
                    let nx = folding_values.get(&(j + domain_size / 2)).unwrap().clone();
                    let v = x
                        + nx
                        + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(*j as u64);
                    let v = v * T::from_int(2).inverse();
                    let v = v + self.challenges[i].pow(2) * (*function_values.get(j).unwrap());
                    if v != *folding_proofs[0].1.get(&j).expect("query missing") {
                        return false;
                    }
                }
            } else {
                for j in &leaf_indices {
                    let x = folding_values.get(j).unwrap().clone();
                    let nx = folding_values.get(&(j + domain_size / 2)).unwrap().clone();
                    let v = x
                        + nx
                        + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(*j as u64);
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

    fn set_final_value(&mut self, final_value: T) {
        self.final_value = Some(final_value);
    }
}

pub struct RollingFriProver<T: Field> {
    total_round: usize,
    coset: Coset<T>,
    function_values: Vec<Vec<T>>,
    folding_values: Vec<Vec<T>>,
    functions_tree: Vec<MerkleTreeProver>,
    folding_tree: Vec<MerkleTreeProver>,
    verifier: Option<Rc<RefCell<RollingFriVerifier<T>>>>,
}

impl<T: Field> RollingFriProver<T> {
    pub fn new(
        total_round: usize,
        function_values: Vec<Vec<T>>,
        coset: &Coset<T>,
    ) -> RollingFriProver<T> {
        RollingFriProver {
            total_round,
            coset: coset.clone(),
            function_values,
            folding_values: vec![vec![]],
            functions_tree: vec![],
            folding_tree: vec![MerkleTreeProver::new(vec![])],
            verifier: None,
        }
    }

    pub fn set_verifier(&mut self, verifier: &Rc<RefCell<RollingFriVerifier<T>>>) {
        self.verifier = Some(verifier.clone());
    }

    fn commit_functions(&mut self) {
        let verifier = self.verifier.clone().unwrap();
        for i in 0..self.function_values.len() {
            let len = self.function_values[i].len();
            let leaf_values: Vec<Vec<u8>> = if i > 0 {
                (0..len)
                .map(|j| as_bytes_vec(&[self.function_values[i][j]]))
                .collect()
            } else {
                (0..len / 2)
                .map(|j| as_bytes_vec(&[self.function_values[i][j], self.function_values[i][j + len / 2]]))
                .collect()
            };
            let leave_number = leaf_values.len();
            let merkle_tree_prover = MerkleTreeProver::new(leaf_values);
            let commit = merkle_tree_prover.commit();
            
            verifier.borrow_mut().set_function_root(leave_number, &commit);
            self.functions_tree.push(merkle_tree_prover);
        }
    }

    fn evaluation_next_domain(
        &self,
        round: usize,
        current_domain: &Coset<T>,
        challenge: T,
    ) -> Vec<T> {
        let mut res = vec![];
        let last_folding_values = if round == 0 {
            &self.function_values[0]
        } else {
            &self.folding_values[round]
        };
        assert_eq!(last_folding_values.len(), current_domain.size());
        let inv_2 = T::from_int(2).inverse();
        let mut shift_inv = current_domain.shift().inverse();
        let generator_inv = current_domain.generator().inverse();
        for i in 0..(last_folding_values.len() / 2) {
            let x = last_folding_values[i];
            let nx = last_folding_values[i + last_folding_values.len() / 2];
            let new_v = ((x + nx) + challenge * (x - nx) * shift_inv) * inv_2;
            if round < self.total_round - 1 {
                res.push(new_v + challenge.pow(2) * self.function_values[round + 1][i]);
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
            let next_evalutation = self.evaluation_next_domain(
                i,
                &domain,
                challenge,
            );
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

    pub fn query(&self, points: &Vec<usize>) -> (Vec<(Vec<u8>, HashMap<usize, T>)>, Vec<(Vec<u8>, HashMap<usize, T>)>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = points.clone();

        for i in 0..self.total_round {
            let len = self.function_values[i].len();

            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                let mut values = HashMap::new();
                for j in &leaf_indices {
                    values.insert(*j, self.function_values[i][*j]);
                    values.insert(j + len / 2, self.function_values[i][*j + len / 2]);
                }
                let proof_bytes = self.functions_tree[i].open(&leaf_indices);
                functions_res.push((proof_bytes, values));
            }

            if i < self.total_round - 1 {
                let mut values = HashMap::new();
                for j in &leaf_indices {
                    values.insert(*j, self.function_values[i + 1][*j]);
                }
                let proof_bytes = self.functions_tree[i + 1].open(&leaf_indices);
                functions_res.push((proof_bytes, values));
            }

            if i > 0 {
                let mut values = HashMap::new();
                for j in &leaf_indices {
                    values.insert(*j, self.folding_values[i][*j]);
                    values.insert(j + len / 2, self.folding_values[i][*j + len / 2]);
                }
                let proof_bytes = self.folding_tree[i].open(&leaf_indices);
                folding_res.push((proof_bytes, values));
            }
        }
        (folding_res, functions_res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::mersenne61_ext::Mersenne61Ext;
    use rand::Rng;

    #[test]
    fn rolling_fri_test() {
        let shift = Mersenne61Ext::random_element();
        let domain = Coset::new(1 << 10, shift);
        let poly_degree_bound = 1 << 8;
        let mut functions = vec![];
        let mut shift = domain.shift();
        let domain_size = domain.size();
        for i in 0..8 {
            let poly = Polynomial::random_polynomial(poly_degree_bound >> i);
            let coset = Coset::new(domain_size >> i, shift);
            functions.push(coset.fft(poly.coefficients()));
            shift *= shift;
        }
        let verifier = Rc::new(RefCell::new(RollingFriVerifier::new(&domain, 8)));
        let prover = Rc::new(RefCell::new(RollingFriProver::new(
            8,
            functions, 
            &domain
        )));
        verifier.borrow_mut().set_prover(&prover);
        prover.borrow_mut().set_verifier(&verifier);

        prover.borrow_mut().commit_functions();
        prover.borrow_mut().prove();
        let mut points = vec![];
        for _i in 0..10 {
            points.push(rand::thread_rng().gen_range(0..domain.size()));
        }
        let (folding_query, function_query) = prover.borrow().query(&points);
        assert!(verifier.borrow().verify(points, folding_query, function_query));
    }
}
