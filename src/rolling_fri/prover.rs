use super::verifier::RollingFriVerifier;
use super::QueryResult;
use crate::algebra::{
    coset::Coset,
    field::{as_bytes_vec, Field},
};
use crate::merkle_tree::MerkleTreeProver;
use std::{cell::RefCell, rc::Rc};

struct InterpolateValue<T: Field> {
    value: Vec<T>,
    merkle_tree: MerkleTreeProver,
}

impl<T: Field> InterpolateValue<T> {
    fn new(value: Vec<T>, half: bool) -> Self {
        let merkle_tree = MerkleTreeProver::new(if half {
            let len = value.len() / 2;
            (0..len)
                .map(|i| as_bytes_vec(&[value[i], value[i + len]]))
                .collect()
        } else {
            let len = value.len();
            (0..len).map(|i| as_bytes_vec(&[value[i]])).collect()
        });
        Self { value, merkle_tree }
    }

    fn query(&self, leaf_indices: &Vec<usize>, half: bool) -> QueryResult<T> {
        if half {
            let len = self.value.len() / 2;
            let proof_values = leaf_indices
                .iter()
                .flat_map(|j| [(*j, self.value[*j]), (*j + len, self.value[*j + len])])
                .collect();
            let proof_bytes = self.merkle_tree.open(&leaf_indices);
            QueryResult {
                proof_bytes,
                proof_values,
            }
        } else {
            let proof_values = leaf_indices.iter().map(|j| (*j, self.value[*j])).collect();
            let proof_bytes = self.merkle_tree.open(&leaf_indices);
            QueryResult {
                proof_bytes,
                proof_values,
            }
        }
    }
}

struct FunctionValue<T: Field> {
    interpolates: Vec<InterpolateValue<T>>,
    map: fn(Vec<T>) -> T,
}

impl<T: Field> FunctionValue<T> {
    fn new(values: Vec<Vec<T>>, map: fn(Vec<T>) -> T, half: bool) -> Self {
        Self {
            interpolates: values
                .into_iter()
                .map(|x| InterpolateValue::new(x, half))
                .collect(),
            map,
        }
    }

    fn get_value(&self, index: usize) -> T {
        let v = self.interpolates.iter().map(|x| x.value[index]).collect();
        let map = self.map;
        map(v)
    }

    fn field_size(&self) -> usize {
        self.interpolates[0].value.len()
    }
}

pub struct RollingFriProver<T: Field> {
    total_round: usize,
    coset: Coset<T>,
    functions: Vec<FunctionValue<T>>,
    foldings: Vec<InterpolateValue<T>>,
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
                .enumerate()
                .map(|(i, (x, y))| FunctionValue::new(x, y, i == 0))
                .collect(),
            foldings: vec![],
            verifier: None,
        }
    }

    pub fn set_verifier(&mut self, verifier: &Rc<RefCell<RollingFriVerifier<T>>>) {
        self.verifier = Some(verifier.clone());
    }

    pub fn commit_functions(&mut self) {
        let verifier = self.verifier.clone().unwrap();
        for (_idx, fv) in self.functions.iter().enumerate() {
            let commit = fv
                .interpolates
                .iter()
                .map(|x| x.merkle_tree.commit())
                .collect();

            verifier.borrow_mut().set_function_root(
                fv.interpolates[0].merkle_tree.leave_num(),
                commit,
                fv.map,
            );
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
            self.functions[round].field_size()
        } else {
            self.foldings[round - 1].value.len()
        };
        let get_folding_value = |i: usize| {
            if round == 0 {
                let fv = &self.functions[round];
                let v = fv.interpolates.iter().map(|x| x.value[i]).collect();
                let map = fv.map;
                map(v)
            } else {
                self.foldings[round - 1].value[i]
            }
        };
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

    pub fn prove(&mut self) {
        let mut domain_size = self.coset.size();
        let mut domain = self.coset.clone();
        let mut shift = domain.shift();
        let verifier = self.verifier.clone().unwrap();
        for i in 0..self.total_round {
            let challenge = verifier.borrow_mut().get_challenge();
            let next_evalutation = self.evaluation_next_domain(i, &domain, challenge);

            shift *= shift;
            domain_size >>= 1;
            domain = Coset::new(domain_size, shift);

            if i < self.total_round - 1 {
                let interpolate_value = InterpolateValue::new(next_evalutation, true);
                let commit = interpolate_value.merkle_tree.commit();
                self.foldings.push(interpolate_value);
                verifier.borrow_mut().receive_root(domain_size / 2, &commit);
            } else {
                verifier.borrow_mut().set_final_value(next_evalutation[0]);
            }
        }
    }

    pub fn query(&self, points: &Vec<usize>) -> (Vec<QueryResult<T>>, Vec<Vec<QueryResult<T>>>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = points.clone();

        for i in 0..self.total_round {
            let len = self.functions[i].field_size();

            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let f = |x: &FunctionValue<T>, y| {
                x.interpolates
                    .iter()
                    .map(|x| x.query(&leaf_indices, y))
                    .collect::<Vec<QueryResult<T>>>()
            };
            if i == 0 {
                let query_result = f(&self.functions[i], true);
                functions_res.push(query_result);
            }

            if i < self.total_round - 1 {
                let query_result = f(&self.functions[i + 1], false);
                functions_res.push(query_result);
            }

            if i > 0 {
                folding_res.push(self.foldings[i - 1].query(&leaf_indices, true));
            }
        }
        (folding_res, functions_res)
    }
}
