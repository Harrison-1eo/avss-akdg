use super::verifier::One2ManyVerifier;
use util::algebra::polynomial::{MultilinearPolynomial, Polynomial};

use util::query_result::QueryResult;
use util::{
    algebra::{
        coset::Coset,
        field::{as_bytes_vec, Field},
    },
    merkle_tree::MerkleTreeProver,
    random_oracle::RandomOracle,
};

#[derive(Clone)]
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

    fn commit(&self) -> [u8; 32] {
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

#[derive(Clone)]
pub struct One2ManyProver<T: Field> {
    total_round: usize,
    variable_num: usize,
    interpolate_cosets: Vec<Coset<T>>,
    functions: Vec<InterpolateValue<T>>,
    foldings: Vec<InterpolateValue<T>>,
    oracle: RandomOracle<T>,
    final_value: Option<Polynomial<T>>,
}

impl<T: Field> One2ManyProver<T> {
    pub fn new(
        total_round: usize,
        interpolate_coset: &Vec<Coset<T>>,
        polynomial: MultilinearPolynomial<T>,
        oracle: &RandomOracle<T>,
    ) -> One2ManyProver<T> {
        let interpolation = interpolate_coset[0].fft(polynomial.coefficients());

        One2ManyProver {
            total_round,
            variable_num: polynomial.variable_num(),
            interpolate_cosets: interpolate_coset.clone(),
            functions: vec![InterpolateValue::new(interpolation)],
            foldings: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn commit_polynomial(&self) -> [u8; 32] {
        assert_eq!(self.functions.len(), 1);
        self.functions[0].commit()
    }

    fn fold(values: &Vec<T>, parameter: T, coset: &Coset<T>) -> Vec<T> {
        let len = values.len() / 2;
        let res = (0..len)
            .into_iter()
            .map(|i| {
                let x = values[i];
                let nx = values[i + len];
                let new_v = (x + nx) + parameter * (x - nx) * coset.element_inv_at(i);
                new_v * T::INVERSE_2
            })
            .collect();
        res
    }

    pub fn commit_functions(&mut self, open_point: &Vec<T>, verifier: &mut One2ManyVerifier<T>) {
        let mut evaluation = None;
        for round in 0..self.total_round {
            let next_evaluation = Self::fold(
                &self.functions[round].value,
                open_point[round],
                &self.interpolate_cosets[round],
            );
            if round < self.total_round - 1 {
                self.functions.push(InterpolateValue::new(next_evaluation));
            } else {
                let mut coefficients = self.interpolate_cosets[round + 1].ifft(&next_evaluation);
                coefficients.truncate(1 << (self.variable_num - self.total_round));
                evaluation = Some(MultilinearPolynomial::new(coefficients));
            }
        }
        for i in 1..self.total_round {
            let function = &self.functions[i];
            verifier.set_function(function.leave_num(), &function.commit());
        }
        verifier.set_evaluation(evaluation.unwrap());
    }

    pub fn commit_foldings(&self, verifier: &mut One2ManyVerifier<T>) {
        for i in 0..(self.total_round - 1) {
            let interpolation = &self.foldings[i];
            verifier.receive_folding_root(interpolation.leave_num(), interpolation.commit());
        }
        verifier.set_final_value(self.final_value.as_ref().unwrap());
    }

    fn evaluation_next_domain(&self, round: usize, challenge: T) -> Vec<T> {
        let mut res = vec![];
        let len = self.interpolate_cosets[round].size();
        let get_folding_value = if round == 0 {
            &self.functions[round]
        } else {
            &self.foldings[round - 1]
        };
        let coset = &self.interpolate_cosets[round];
        for i in 0..(len / 2) {
            let x = get_folding_value.value[i];
            let nx = get_folding_value.value[i + len / 2];
            let new_v = (x + nx) + challenge * (x - nx) * coset.element_inv_at(i);
            if round == 0 {
                res.push(new_v);
            } else {
                let fv = &self.functions[round];
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
                let next_evalutation = self.evaluation_next_domain(i, challenge);
                self.foldings.push(InterpolateValue::new(next_evalutation));
            } else {
                let next_evalutation = self.evaluation_next_domain(i, challenge);
                let coefficients = self.interpolate_cosets[i + 1].ifft(&next_evalutation);
                self.final_value = Some(Polynomial::new(coefficients));
            }
        }
    }

    pub fn query(&self) -> (Vec<QueryResult<T>>, Vec<QueryResult<T>>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = self.oracle.query_list.clone();

        for i in 0..self.total_round {
            let len = self.interpolate_cosets[i].size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let query_result = self.functions[i].query(&leaf_indices);
            functions_res.push(query_result);

            if i > 0 {
                folding_res.push(self.foldings[i - 1].query(&leaf_indices));
            }
        }
        (folding_res, functions_res)
    }
}
