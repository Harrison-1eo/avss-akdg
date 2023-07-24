use std::{cell::RefCell, rc::Rc};

use super::verifier::One2ManyVerifier;
use util::algebra::polynomial::Polynomial;

use util::merkle_tree::MERKLE_ROOT_SIZE;
use util::query_result::QueryResult;
use util::{
    algebra::{
        coset::Coset,
        field::{as_bytes_vec, Field},
    },
    merkle_tree::MerkleTreeProver,
    random_oracle::RandomOracle,
};

struct InterpolateValue<T: Field> {
    value: Vec<T>,
    merkle_tree: MerkleTreeProver,
}

impl<T: Field> InterpolateValue<T> {
    fn new(value: Vec<T>) -> Self {
        let len = value.len() / 2;
        let merkle_tree = MerkleTreeProver::new(
            (0..len)
                .map(|i| as_bytes_vec(&[value[i], value[i + len]]))
                .collect(),
        );
        Self { value, merkle_tree }
    }

    fn leave_num(&self) -> usize {
        self.merkle_tree.leave_num()
    }

    fn commit(&self) -> [u8; MERKLE_ROOT_SIZE] {
        self.merkle_tree.commit()
    }

    fn query(&self, leaf_indices: &Vec<usize>) -> QueryResult<T> {
        let len = self.merkle_tree.leave_num();
        let proof_values = leaf_indices
            .iter()
            .flat_map(|j| [(*j, self.value[*j]), (*j + len, self.value[*j + len])])
            .collect();
        let proof_bytes = self.merkle_tree.open(&leaf_indices);
        QueryResult {
            proof_bytes,
            proof_values,
        }
    }
}

struct CosetInterpolate<T: Field> {
    interpolates: Vec<InterpolateValue<T>>,
}

impl<T: Field> CosetInterpolate<T> {
    fn len(&self) -> usize {
        self.interpolates.len()
    }
    fn new(functions: Vec<Vec<T>>) -> Self {
        CosetInterpolate {
            interpolates: functions
                .into_iter()
                .map(|values| InterpolateValue::new(values))
                .collect(),
        }
    }
    fn field_size(&self) -> usize {
        self.interpolates[0].value.len()
    }
    fn from_interpolates(interpolates: Vec<InterpolateValue<T>>) -> Self {
        CosetInterpolate { interpolates }
    }

    fn get_interpolation(&self, index: usize) -> &InterpolateValue<T> {
        let len = self.interpolates.len();
        assert!((len & (len - 1)) == 0);
        &self.interpolates[index & (len - 1)]
    }
}

pub struct One2ManyProver<T: Field> {
    total_round: usize,
    interpolate_cosets: Vec<Coset<T>>,
    functions: Vec<CosetInterpolate<T>>,
    foldings: Vec<CosetInterpolate<T>>,
    oracle: RandomOracle<T>,
    final_value: Vec<Polynomial<T>>,
}

impl<T: Field> One2ManyProver<T> {
    pub fn new(
        total_round: usize,
        interpolate_coset: &Vec<Coset<T>>,
        functions: Vec<Vec<Vec<T>>>,
        oracle: &RandomOracle<T>,
    ) -> One2ManyProver<T> {
        assert_eq!(total_round, functions.len());
        let functions: Vec<CosetInterpolate<T>> = functions
            .into_iter()
            .map(|x| CosetInterpolate::new(x))
            .collect();

        One2ManyProver {
            total_round,
            interpolate_cosets: interpolate_coset.clone(),
            functions,
            foldings: vec![],
            oracle: oracle.clone(),
            final_value: vec![],
        }
    }

    pub fn commit_functions(&self, verifiers: &Vec<Rc<RefCell<One2ManyVerifier<T>>>>) {
        for i in 0..self.total_round {
            for (idx, j) in verifiers.into_iter().enumerate() {
                let function = self.functions[i].get_interpolation(idx);
                j.borrow_mut()
                    .set_function(function.leave_num(), &function.commit());
            }
        }
    }

    pub fn commit_foldings(&self, verifiers: &Vec<Rc<RefCell<One2ManyVerifier<T>>>>) {
        for i in 0..(self.total_round - 1) {
            for (idx, j) in verifiers.into_iter().enumerate() {
                let interpolation = self.foldings[i].get_interpolation(idx);
                j.borrow_mut()
                    .receive_folding_root(interpolation.leave_num(), interpolation.commit());
            }
        }
        for i in 0..verifiers.len() {
            verifiers[i]
                .borrow_mut()
                .set_final_value(&self.final_value[i % self.final_value.len()]);
        }
    }

    fn evaluation_next_domain(
        &self,
        round: usize,
        rolling_function_index: usize,
        challenge: T,
    ) -> Vec<T> {
        let mut res = vec![];
        let len = self.functions[round].field_size();
        let get_folding_value = if round == 0 {
            self.functions[round].get_interpolation(rolling_function_index)
        } else {
            self.foldings[round - 1].get_interpolation(rolling_function_index)
        };
        let coset = &self.interpolate_cosets[round];
        for i in 0..(len / 2) {
            let x = get_folding_value.value[i];
            let nx = get_folding_value.value[i + len / 2];
            let new_v = (x + nx) + challenge * (x - nx) * coset.element_inv_at(i);
            if round == 0 {
                res.push(new_v);
            } else {
                let fv = &self.functions[round].interpolates[rolling_function_index];
                let x = fv.value[i];
                let nx = fv.value[i + len / 2];
                let new_v =
                    (new_v * challenge + (x + nx)) * challenge + (x - nx) * coset.element_inv_at(i);
                res.push(new_v);
            }
        }
        res
    }

    pub fn prove(&mut self) {
        for i in 0..self.total_round {
            let challenge = self.oracle.folding_challenges[i];
            if i < self.total_round - 1 {
                let mut interpolates = vec![];
                for j in 0..self.functions[i].len() {
                    let next_evalutation = self.evaluation_next_domain(i, j, challenge);
                    let interpolate_value = InterpolateValue::new(next_evalutation);
                    interpolates.push(interpolate_value);
                }
                self.foldings
                    .push(CosetInterpolate::from_interpolates(interpolates));
            } else {
                for j in 0..self.functions[i].len() {
                    let next_evalutation = self.evaluation_next_domain(i, j, challenge);
                    let coefficients = self.interpolate_cosets[i + 1].ifft(next_evalutation);
                    self.final_value.push(Polynomial::new(coefficients));
                }
            }
        }
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = self.oracle.query_list.clone();

        for i in 0..self.total_round {
            let len = self.functions[i].field_size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                let query_result = self.functions[0].get_interpolation(0).query(&leaf_indices);
                functions_res.push(vec![query_result]);
            } else {
                let query_result = self.functions[i]
                    .interpolates
                    .iter()
                    .map(|x| x.query(&leaf_indices))
                    .collect();
                functions_res.push(query_result);
            }

            if i > 0 {
                folding_res.push(
                    self.foldings[i - 1]
                        .interpolates
                        .iter()
                        .map(|x| x.query(&leaf_indices))
                        .collect(),
                );
            }
        }
        (folding_res, functions_res)
    }
}
