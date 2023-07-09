use std::{cell::RefCell, rc::Rc};

use super::verifier::One2ManyVerifier;
use crate::random_oracle::RandomOracle;

use crate::util::QueryResult;
use crate::{
    algebra::{
        coset::Coset,
        field::{as_bytes_vec, Field},
    },
    merkle_tree::MerkleTreeProver,
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

    fn field_size(&self) -> usize {
        self.value.len()
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

struct FunctionValue<T: Field> {
    interpolate: InterpolateValue<T>,
    map: Box<dyn Fn(T, T, T) -> T>, // value, x, challenge
}

impl<T: Field> FunctionValue<T> {
    fn new(value: Vec<T>, map: Box<dyn Fn(T, T, T) -> T>) -> Self {
        Self {
            interpolate: InterpolateValue::new(value),
            map,
        }
    }

    fn get_value(&self, index: usize, coset: &Coset<T>, challenge: T) -> T {
        let v = self.interpolate.value[index];
        (self.map)(v, coset.element_at(index), challenge)
    }

    fn commit(&self) -> [u8; 32] {
        self.interpolate.merkle_tree.commit()
    }

    fn leave_number(&self) -> usize {
        self.interpolate.merkle_tree.leave_num()
    }

    fn field_size(&self) -> usize {
        self.interpolate.value.len()
    }
}

struct CosetFunction<T: Field> {
    pub functions: Vec<FunctionValue<T>>,
}

impl<T: Field> CosetFunction<T> {
    fn new(functions: Vec<(Vec<T>, Box<dyn Fn(T, T, T) -> T>)>) -> Self {
        CosetFunction {
            functions: functions
                .into_iter()
                .map(|(values, map)| FunctionValue::new(values, map))
                .collect(),
        }
    }

    fn field_size(&self) -> usize {
        self.functions[0].field_size()
    }

    fn get_function(&self, index: usize) -> &FunctionValue<T> {
        let len = self.functions.len();
        assert_eq!(len & (len - 1), 0);
        &self.functions[index % len]
    }

    fn len(&self) -> usize {
        self.functions.len()
    }
}

struct CosetInterpolate<T: Field> {
    interpolates: Vec<InterpolateValue<T>>,
}

impl<T: Field> CosetInterpolate<T> {
    fn from_interpolates(interpolates: Vec<InterpolateValue<T>>) -> Self {
        let len = interpolates[0].field_size();
        let leave_num = interpolates[0].leave_num();
        for i in &interpolates {
            assert_eq!(i.value.len(), len);
            assert_eq!(i.leave_num(), leave_num);
        }
        CosetInterpolate { interpolates }
    }

    fn field_size(&self) -> usize {
        self.interpolates[0].field_size()
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
    functions: Vec<CosetFunction<T>>,
    foldings: Vec<CosetInterpolate<T>>,
    oracle: Rc<RefCell<RandomOracle<T>>>,
    final_value: Vec<T>,
}

impl<T: Field> One2ManyProver<T> {
    pub fn new(
        total_round: usize,
        interpolate_coset: &Coset<T>,
        functions: Vec<Vec<(Vec<T>, Box<dyn Fn(T, T, T) -> T>)>>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
    ) -> One2ManyProver<T> {
        assert_eq!(total_round, functions.len());
        let functions: Vec<CosetFunction<T>> = functions
            .into_iter()
            .map(|x| CosetFunction::new(x))
            .collect();
        let mut cosets = vec![interpolate_coset.clone()];
        for _ in 1..total_round {
            cosets.push(cosets.last().as_ref().unwrap().pow(2));
        }

        One2ManyProver {
            total_round,
            interpolate_cosets: cosets,
            functions,
            foldings: vec![],
            oracle: oracle.clone(),
            final_value: vec![],
        }
    }

    pub fn commit_functions(&self, verifiers: &Vec<Rc<RefCell<One2ManyVerifier<T>>>>) {
        for i in 0..self.total_round {
            for (idx, j) in verifiers.into_iter().enumerate() {
                let function = self.functions[i].get_function(idx);
                j.borrow_mut()
                    .set_function(function.leave_number(), &function.commit());
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
                .set_final_value(self.final_value[i % self.final_value.len()]);
        }
    }

    fn evaluation_next_domain(
        &self,
        round: usize,
        rolling_function_index: usize,
        challenge: T,
    ) -> Vec<T> {
        let mut res = vec![];
        let len = if round == 0 {
            self.functions[round].field_size()
        } else {
            self.foldings[round - 1].field_size()
        };
        let get_folding_value = |i: usize| {
            if round == 0 {
                self.functions[round].get_function(0).get_value(
                    i,
                    &self.interpolate_cosets[round],
                    challenge,
                )
            } else {
                self.foldings[round - 1]
                    .get_interpolation(rolling_function_index)
                    .value[i]
            }
        };
        let mut shift_inv = self.interpolate_cosets[round].shift().inverse();
        let generator_inv = self.interpolate_cosets[round].generator().inverse();
        for i in 0..(len / 2) {
            let x = get_folding_value(i);
            let nx = get_folding_value(i + len / 2);
            let new_v = (x + nx) + challenge * (x - nx) * shift_inv;
            if round == 0 || round == self.total_round - 1 {
                res.push(new_v);
            } else {
                let fv = &self.functions[round].functions[rolling_function_index];
                let x = fv.get_value(i, &self.interpolate_cosets[round], challenge);
                let nx = fv.get_value(i + len / 2, &self.interpolate_cosets[round], challenge);
                let new_v = (new_v * challenge + (x + nx)) * challenge + (x - nx) * shift_inv;
                res.push(new_v);
            }
            shift_inv *= generator_inv;
        }
        res
    }

    pub fn prove(&mut self) {
        for i in 0..self.total_round {
            let challenge = self.oracle.borrow_mut().generate_challenge();
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
                    self.final_value.push(next_evalutation[0]);
                }
            }
        }
    }

    pub fn query(&self) -> (Vec<Vec<QueryResult<T>>>, Vec<Vec<QueryResult<T>>>) {
        let mut folding_res = vec![];
        let mut functions_res = vec![];
        let mut leaf_indices = self.oracle.borrow().query_list();

        for i in 0..self.total_round {
            let len = self.functions[i].field_size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let f = |x: &FunctionValue<T>| x.interpolate.query(&leaf_indices);
            if i == 0 {
                let query_result = f(self.functions[0].get_function(0));
                functions_res.push(vec![query_result]);
            } else {
                let query_result = self.functions[i].functions.iter().map(|x| f(x)).collect();
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
