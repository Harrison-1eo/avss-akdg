use super::verifier::FriVerifier;
use util::{
    algebra::{field, polynomial::MultilinearPolynomial},
    random_oracle::RandomOracle,
};

use util::query_result::QueryResult;
use util::{
    algebra::{
        coset::Coset,
        field::{as_bytes_vec, Field},
    },
    merkle_tree::MerkleTreeProver,
};

use super::Tuple;

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
pub struct Function<T: Field> {
    interpolation: InterpolateValue<T>,
    evaluations: Vec<(T, T)>,
}

impl<T: Field> Function<T> {
    pub fn new(value: Vec<T>, evaluations: Vec<(T, T)>) -> Self {
        Function {
            interpolation: InterpolateValue::new(value),
            evaluations,
        }
    }
}

#[derive(Clone)]
pub struct FriProver<T: Field> {
    total_round: usize,
    interpolate_cosets: Vec<Coset<T>>,
    functions: Vec<Function<T>>,
    polynomials: Vec<MultilinearPolynomial<T>>,
    foldings: Vec<InterpolateValue<T>>,
    oracle: RandomOracle<T>,
    final_value: Option<T>,
}

impl<T: Field> FriProver<T> {
    pub fn new(
        total_round: usize,
        interpolate_coset: &Vec<Coset<T>>,
        polynomial: MultilinearPolynomial<T>,
        oracle: &RandomOracle<T>,
    ) -> FriProver<T> {
        let interpolation = interpolate_coset[0].fft(polynomial.coefficients());
        FriProver {
            total_round,
            interpolate_cosets: interpolate_coset.clone(),
            functions: vec![Function::new(interpolation, vec![])],
            polynomials: vec![polynomial],
            foldings: vec![],
            oracle: oracle.clone(),
            final_value: None,
        }
    }

    pub fn compute_tuples(&mut self) -> Vec<Tuple<T>> {
        let beta = self.oracle.beta;
        let mut tuples = vec![];
        for i in 0..self.total_round {
            tuples.push(Tuple {
                a: self.polynomials[i].evaluate_as_polynomial(beta),
                b: self.polynomials[i].evaluate_as_polynomial(-beta),
                c: self.polynomials[i + 1].evaluate_as_polynomial(beta * beta),
            });
        }
        self.functions[0].evaluations = vec![(beta, tuples[0].a), (-beta, tuples[0].b)];
        for i in 1..tuples.len() {
            self.functions[i].evaluations = vec![
                (beta * beta, tuples[i - 1].c),
                (beta, tuples[i].a),
                (-beta, tuples[i].b),
            ];
        }
        tuples
    }

    pub fn commit_first_polynomial(&self) -> [u8; 32] {
        self.functions[0].interpolation.commit()
    }

    pub fn commit_functions(&mut self, verifier: &mut FriVerifier<T>, open_point: &Vec<T>) {
        assert_eq!(open_point.len(), self.total_round);
        for i in 0..self.total_round {
            let last = self.polynomials.last().unwrap();
            let next_polynomial = last.folding(open_point[i]);
            self.polynomials.push(next_polynomial);
        }
        for i in 1..self.total_round {
            self.functions.push(Function::new(
                self.interpolate_cosets[0].fft(self.polynomials[i].coefficients()),
                vec![],
            ));
            verifier.append_function(self.functions[i].interpolation.commit());
        }
    }

    pub fn commit_foldings(&self, verifier: &mut FriVerifier<T>) {
        for i in 0..(self.total_round - 1) {
            verifier.receive_folding_root(self.foldings[i].leave_num(), self.foldings[i].commit());
        }
        verifier.set_final_value(self.final_value.unwrap())
    }

    fn initial_interpolation(&self) -> Vec<T> {
        let rlc = self.oracle.rlc;
        let mut acc = rlc;
        let mut res = self.functions[0].interpolation.value.clone();
        for i in 0..self.functions.len() {
            let interpolation = &self.functions[i].interpolation.value;
            if i != 0 {
                for j in 0..interpolation.len() {
                    res[j] += interpolation[j] * acc;
                }
                acc *= rlc;
            }
            for (x, y) in &self.functions[i].evaluations {
                let v = self.interpolate_cosets[0]
                    .all_elements()
                    .into_iter()
                    .map(|i| i - *x)
                    .collect();

                let inv = field::batch_inverse(&v);
                for j in 0..interpolation.len() {
                    res[j] += inv[j] * (interpolation[j] - *y) * acc;
                }
                acc *= rlc
            }
        }
        res
    }

    fn evaluation_next_domain(&self, round: usize, challenge: T) -> Vec<T> {
        let mut res = vec![];
        let coset = &self.interpolate_cosets[round];
        let len = coset.size();
        if round == 0 {
            let function = self.initial_interpolation();
            for i in 0..(len / 2) {
                let x = function[i];
                let nx = function[i + len / 2];
                let new_v = (x + nx) + challenge * (x - nx) * coset.element_inv_at(i);
                res.push(new_v);
            }
        } else {
            let last_folding = &self.foldings.last().unwrap().value;
            for i in 0..(len / 2) {
                let x = last_folding[i];
                let nx = last_folding[i + len / 2];
                let new_v = (x + nx) + challenge * (x - nx) * coset.element_inv_at(i);
                res.push(new_v);
            }
        }
        res
    }

    pub fn prove(&mut self) {
        for i in 0..self.total_round {
            let challenge = self.oracle.folding_challenges[i];
            let next_evalutation = self.evaluation_next_domain(i, challenge);
            if i < self.total_round - 1 {
                let interpolate_value = InterpolateValue::new(next_evalutation);
                self.foldings.push(interpolate_value);
            } else {
                self.final_value = Some(next_evalutation[0]);
            }
        }
    }

    pub fn query(&self) -> (Vec<QueryResult<T>>, Vec<QueryResult<T>>) {
        let mut folding_res = vec![];
        let mut functions_res = None;
        let mut leaf_indices = self.oracle.query_list.clone();

        for i in 0..self.total_round {
            let len = self.interpolate_cosets[i].size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            if i == 0 {
                functions_res = Some(
                    self.functions
                        .iter()
                        .map(|x| x.interpolation.query(&leaf_indices))
                        .collect(),
                );
            } else {
                let query_result = self.foldings[i - 1].query(&leaf_indices);
                folding_res.push(query_result);
            }
        }
        (folding_res, functions_res.unwrap())
    }
}
