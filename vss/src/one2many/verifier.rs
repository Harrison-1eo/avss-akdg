use util::algebra::polynomial::{MultilinearPolynomial, Polynomial};
use util::random_oracle::RandomOracle;
use util::{
    algebra::{coset::Coset, field::Field},
    merkle_tree::MerkleTreeVerifier,
    query_result::QueryResult,
};

#[derive(Clone)]
pub struct One2ManyVerifier<T: Field> {
    total_round: usize,
    log_max_degree: usize,
    interpolate_cosets: Vec<Coset<T>>,
    function_root: Vec<MerkleTreeVerifier>,
    folding_root: Vec<MerkleTreeVerifier>,
    oracle: RandomOracle<T>,
    final_value: Option<Polynomial<T>>,
}

impl<T: Field> One2ManyVerifier<T> {
    pub fn new_with_default_map(
        total_round: usize,
        log_max_degree: usize,
        coset: &Vec<Coset<T>>,
        oracle: &RandomOracle<T>,
    ) -> Self {
        One2ManyVerifier {
            total_round,
            log_max_degree,
            interpolate_cosets: coset.clone(),
            function_root: vec![],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn new(
        total_round: usize,
        log_max_degree: usize,
        coset: &Vec<Coset<T>>,
        oracle: &RandomOracle<T>,
    ) -> Self {
        One2ManyVerifier {
            total_round,
            log_max_degree,
            interpolate_cosets: coset.clone(),
            function_root: vec![],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
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

    pub fn set_final_value(&mut self, value: &Polynomial<T>) {
        assert!(value.degree() <= 1 << (self.log_max_degree - self.total_round));
        self.final_value = Some(value.clone());
    }

    fn verify_both_condition(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
        extra_folding_param: Option<&Vec<T>>,
        extra_final_poly: Option<&MultilinearPolynomial<T>>,
    ) -> bool {
        let has_extra = match extra_final_poly {
            Some(_) => true,
            None => false,
        };
        let mut leaf_indices = self.oracle.query_list.clone();
        for i in 0..self.total_round {
            let domain_size = self.interpolate_cosets[i].size();
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

            let challenge = self.oracle.folding_challenges[i];
            let get_folding_value = if i == 0 {
                &function_proofs[i].proof_values
            } else {
                &folding_proofs[i - 1].proof_values
            };

            let function_values = if i != 0 {
                let function_query_result = &function_proofs[i];
                function_query_result.verify_merkle_tree(&leaf_indices, &self.function_root[i]);
                Some(&function_query_result.proof_values)
            } else {
                None
            };
            for j in &leaf_indices {
                let x = get_folding_value[j];
                let nx = get_folding_value[&(j + domain_size / 2)];
                let v =
                    x + nx + challenge * (x - nx) * self.interpolate_cosets[i].element_inv_at(*j);
                if i != 0 {
                    let x = function_values.as_ref().unwrap()[j];
                    let nx = function_values.as_ref().unwrap()[&(j + domain_size / 2)];
                    let v = (v * challenge + (x + nx)) * challenge
                        + (x - nx) * self.interpolate_cosets[i].element_inv_at(*j);
                    if i == self.total_round - 1 {
                        let x = self.interpolate_cosets[i + 1].element_at(*j);
                        if v != self.final_value.as_ref().unwrap().evaluation_at(x) {
                            return false;
                        }
                    } else if v != folding_proofs[i].proof_values[j] {
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
                            * self.interpolate_cosets[i].element_inv_at(*j);
                    if i < self.total_round - 1 {
                        assert_eq!(v, function_proofs[i + 1].proof_values[j] * T::from_int(2));
                    } else {
                        let x = self.interpolate_cosets[i + 1].element_at(*j);
                        let poly_v = extra_final_poly.unwrap().evaluate_as_polynomial(x);
                        assert_eq!(v, poly_v * T::from_int(2));
                    }
                }
            }
        }
        true
    }

    pub fn verify_with_extra_folding(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
        extra_folding_param: &Vec<T>,
        extra_final_poly: &MultilinearPolynomial<T>,
    ) -> bool {
        self.verify_both_condition(
            folding_proofs,
            function_proofs,
            Some(extra_folding_param),
            Some(extra_final_poly),
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
