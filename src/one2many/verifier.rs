use crate::random_oracle::RandomOracle;
use crate::util::QueryResult;
use crate::{
    algebra::{coset::Coset, field::Field},
    merkle_tree::MerkleTreeVerifier,
};
use std::{cell::RefCell, rc::Rc};

#[derive(Clone)]
pub struct One2ManyVerifier<T: Field> {
    total_round: usize,
    interpolate_cosets: Vec<Coset<T>>,
    function_root: Vec<MerkleTreeVerifier>,
    function_maps: Vec<Rc<dyn Fn(T, T, T) -> T>>,
    folding_root: Vec<MerkleTreeVerifier>,
    oracle: Rc<RefCell<RandomOracle<T>>>,
    final_value: Option<T>,
}

impl<T: Field> One2ManyVerifier<T> {
    pub fn new_with_default_map(
        total_round: usize,
        coset: &Coset<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
    ) -> Self {
        let mut cosets = vec![coset.clone()];
        for _ in 1..total_round {
            cosets.push(cosets.last().as_ref().unwrap().pow(2));
        }
        One2ManyVerifier {
            total_round,
            interpolate_cosets: cosets,
            function_root: vec![],
            function_maps: (0..total_round)
                .map(|_| Rc::new(move |v: T, _: T, _: T| v) as Rc<dyn Fn(T, T, T) -> T>)
                .collect(),
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn new(
        total_round: usize,
        coset: &Coset<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
    ) -> Self {
        let mut cosets = vec![coset.clone()];
        for _ in 1..total_round {
            cosets.push(cosets.last().as_ref().unwrap().pow(2));
        }
        One2ManyVerifier {
            total_round,
            interpolate_cosets: cosets,
            function_root: vec![],
            function_maps: vec![],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn set_map(&mut self, function_map: Rc<dyn Fn(T, T, T) -> T>) {
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

    fn verify_both_condition(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
        extra_folding_param: Option<&Vec<T>>,
        extra_final_value: Option<T>,
    ) -> bool {
        let has_extra = match extra_final_value {
            Some(_) => true,
            None => false,
        };
        let mut leaf_indices = self.oracle.borrow().query_list();
        let mut shift_inv = self.interpolate_cosets[0].shift().inverse();
        let mut generator_inv = self.interpolate_cosets[0].generator().inverse();
        let mut domain_size = self.interpolate_cosets[0].size();
        for i in 0..self.total_round {
            leaf_indices = leaf_indices
                .iter_mut()
                .map(|v| *v % (domain_size >> 1))
                .collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                function_proofs[i].verify_merkle_tree(&leaf_indices, &self.function_root[0]);
            } else {
                folding_proofs[i - 1].verify_merkle_tree(&leaf_indices, &self.folding_root[i - 1]);
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

            let function_values = if i != 0 {
                let function_query_result = &function_proofs[i];
                function_query_result.verify_merkle_tree(&leaf_indices, &self.function_root[i]);
                Some(&function_query_result.proof_values)
            } else {
                None
            };
            for j in &leaf_indices {
                let x = get_folding_value(j);
                let nx = get_folding_value(&(j + domain_size / 2));
                let v = x + nx + challenge * (x - nx) * shift_inv * generator_inv.pow(*j);
                if i == self.total_round - 1 {
                    if v != self.final_value.unwrap() {
                        return false;
                    }
                } else if i != 0 {
                    let x = self.function_maps[i](
                        function_values.as_ref().unwrap()[j],
                        self.interpolate_cosets[i].all_elements()[*j],
                        challenge,
                    );
                    let nx = self.function_maps[i](
                        function_values.as_ref().unwrap()[&(j + domain_size / 2)],
                        self.interpolate_cosets[i].all_elements()[*j + domain_size / 2],
                        challenge,
                    );
                    let v = (v * challenge + (x + nx)) * challenge
                        + (x - nx) * shift_inv * generator_inv.pow(*j);
                    if v != folding_proofs[i].proof_values[j] {
                        return false;
                    }
                } else {
                    if v != folding_proofs[i].proof_values[j] {
                        return false;
                    }
                }
                if has_extra {
                    let x = function_proofs[i].proof_values[j];
                    let nx = function_proofs[i].proof_values[&(j + domain_size / 2)];
                    let v = x
                        + nx
                        + extra_folding_param.unwrap()[i]
                            * (x - nx)
                            * shift_inv
                            * generator_inv.pow(*j);
                    if i < self.total_round - 1 {
                        assert_eq!(v, function_proofs[i + 1].proof_values[j] * T::from_int(2));
                    } else {
                        assert_eq!(v, extra_final_value.unwrap() * T::from_int(2));
                    }
                }
            }

            shift_inv *= shift_inv;
            generator_inv *= generator_inv;
            domain_size >>= 1;
        }
        true
    }

    pub fn verify_with_extra_folding(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
        extra_folding_param: &Vec<T>,
        extra_final_value: T,
    ) -> bool {
        self.verify_both_condition(
            folding_proofs,
            function_proofs,
            Some(extra_folding_param),
            Some(extra_final_value),
        )
    }

    pub fn verify(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
    ) -> bool {
        self.verify_both_condition(folding_proofs, function_proofs, None, None)
    }
}
