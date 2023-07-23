use super::Tuple;
use util::query_result::QueryResult;
use util::random_oracle::RandomOracle;
use util::{
    algebra::{coset::Coset, field::Field},
    merkle_tree::MerkleTreeVerifier,
};

#[derive(Clone)]
pub struct FriVerifier<T: Field> {
    total_round: usize,
    interpolate_cosets: Vec<Coset<T>>,
    function_root: Vec<(MerkleTreeVerifier, Vec<(T, T)>)>,
    folding_root: Vec<MerkleTreeVerifier>,
    oracle: RandomOracle<T>,
    final_value: Option<T>,
    open_point: Option<Vec<T>>,
}

impl<T: Field> FriVerifier<T> {
    pub fn new(
        total_round: usize,
        coset: &Vec<Coset<T>>,
        polynomial_commitment: [u8; 32],
        oracle: &RandomOracle<T>,
    ) -> Self {
        FriVerifier {
            total_round,
            interpolate_cosets: coset.clone(),
            function_root: vec![(
                MerkleTreeVerifier {
                    leave_number: coset[0].size() / 2,
                    merkle_root: polynomial_commitment,
                },
                vec![],
            )],
            folding_root: vec![],
            oracle: oracle.clone(),
            final_value: None,
            open_point: None,
        }
    }

    pub fn get_open_point(&mut self) -> Vec<T> {
        let point = (0..self.total_round)
            .map(|_| T::random_element())
            .collect::<Vec<T>>();
        self.open_point = Some(point.clone());
        point
    }

    pub fn append_function(&mut self, function_root: [u8; 32]) {
        self.function_root.push((
            MerkleTreeVerifier {
                merkle_root: function_root,
                leave_number: self.interpolate_cosets[0].size() / 2,
            },
            vec![],
        ));
    }

    pub fn set_tuples(&mut self, tuples: &Vec<Tuple<T>>) {
        let beta = self.oracle.beta;
        for i in 0..tuples.len() {
            assert!(tuples[i].verify(beta, self.open_point.as_ref().unwrap()[i]));
            self.function_root[i].1.push((beta, tuples[i].a));
            self.function_root[i].1.push((-beta, tuples[i].b));
            if i < tuples.len() - 1 {
                self.function_root[i + 1].1.push((beta * beta, tuples[i].c))
            }
        }
    }

    pub fn receive_folding_root(&mut self, leave_number: usize, folding_root: [u8; 32]) {
        self.folding_root.push(MerkleTreeVerifier {
            leave_number,
            merkle_root: folding_root,
        });
    }

    pub fn set_final_value(&mut self, value: T) {
        assert_ne!(value, T::from_int(0));
        self.final_value = Some(value);
    }

    pub fn verify(
        &self,
        folding_proofs: &Vec<QueryResult<T>>,
        function_proofs: &Vec<QueryResult<T>>,
    ) -> bool {
        let mut leaf_indices = self.oracle.query_list.clone();
        let rlc = self.oracle.rlc;
        for i in 0..self.total_round {
            let domain_size = self.interpolate_cosets[i].size();
            leaf_indices = leaf_indices
                .iter_mut()
                .map(|v| *v % (domain_size >> 1))
                .collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                for j in 0..function_proofs.len() {
                    assert!(function_proofs[j]
                        .verify_merkle_tree(&leaf_indices, &self.function_root[j].0));
                }
            } else {
                folding_proofs[i - 1].verify_merkle_tree(&leaf_indices, &self.folding_root[i - 1]);
            }

            let challenge = self.oracle.folding_challenges[i];
            let get_folding_value = |index: &usize| {
                if i == 0 {
                    let mut tmp_rlc = T::from_int(1);
                    let mut res = T::from_int(0); //function_proofs[0].proof_values[index];
                    for f in 0..self.function_root.len() {
                        let this_v = function_proofs[f].proof_values[index];
                        res += this_v * tmp_rlc;
                        tmp_rlc *= rlc;
                        for (x, y) in &self.function_root[f].1 {
                            res += tmp_rlc
                                * (this_v - *y)
                                * (self.interpolate_cosets[0].element_at(*index) - *x).inverse();
                            tmp_rlc *= rlc
                        }
                    }
                    res
                } else {
                    folding_proofs[i - 1].proof_values[index]
                }
            };

            for j in &leaf_indices {
                let x = get_folding_value(j);
                let nx = get_folding_value(&(j + domain_size / 2));
                let v =
                    x + nx + challenge * (x - nx) * self.interpolate_cosets[i].element_inv_at(*j);
                if i < self.total_round - 1 {
                    if v != folding_proofs[i].proof_values[j] {
                        return false;
                    }
                } else {
                    if v != self.final_value.unwrap() {
                        return false;
                    }
                }
            }
        }
        true
    }
}
