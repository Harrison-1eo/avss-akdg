use crate::algebra::{
    coset::Coset,
    field::{as_bytes_vec, Field},
    polynomial::*,
};
use crate::merkle_tree::{MerkleTreeProver, MerkleTreeVerifier};
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct FriProver<T: Field> {
    log_poly_degree: usize,
    coset: Coset<T>,
    interpolate_values: Vec<Vec<T>>,
    merkle_tree: Vec<MerkleTreeProver>,
    verifier: Option<Rc<RefCell<FriVerifier<T>>>>,
}

pub struct FriVerifier<T: Field> {
    coset: Coset<T>,
    log_poly_degree: usize,
    challenges: Vec<T>,
    merkle_root: Vec<MerkleTreeVerifier>,
    prover: Option<Rc<RefCell<FriProver<T>>>>,
    final_value: Option<T>,
}

impl<T: Field> FriVerifier<T> {
    pub fn new(coset: &Coset<T>, log_poly_degree: usize) -> FriVerifier<T> {
        FriVerifier {
            coset: coset.clone(),
            log_poly_degree,
            challenges: vec![],
            merkle_root: vec![],
            prover: None,
            final_value: None,
        }
    }

    fn receive_root(&mut self, leave_number: usize, merkle_root: &[u8; 32]) {
        self.merkle_root.push(MerkleTreeVerifier {
            leave_number,
            merkle_root: merkle_root.clone(),
        })
    }

    fn get_challenge(&mut self) -> T {
        let challenge = T::random_element();
        self.challenges.push(challenge);
        challenge
    }

    pub fn set_prover(&mut self, prover: &Rc<RefCell<FriProver<T>>>) {
        self.prover = Some(prover.clone());
    }

    pub fn verify(
        &self,
        mut leaf_indices: Vec<usize>,
        mut evaluation: Vec<(Vec<u8>, HashMap<usize, T>)>,
    ) -> bool {
        let mut shift_inv = self.coset.shift().inverse();
        let mut generator_inv = self.coset.generator().inverse();
        let mut len = self.coset.num_elements();
        for i in 0..self.log_poly_degree {
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();
            let (proof_bytes, values) = evaluation.remove(0);

            let open_values = leaf_indices
                .iter()
                .map(|v| {
                    as_bytes_vec(&[
                        values.get(v).unwrap().clone(),
                        values.get(&(v + len / 2)).unwrap().clone(),
                    ])
                })
                .collect();
            if !self.merkle_root[i].verify(proof_bytes, &leaf_indices, &open_values) {
                return false;
            }

            for j in &leaf_indices {
                let x = values.get(j).unwrap().clone();
                let nx = values.get(&(j + len / 2)).unwrap().clone();
                let v = x
                    + nx
                    + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(*j as u64);
                let v = v * T::from_int(2).inverse();
                if i < self.log_poly_degree - 1 {
                    if v != *evaluation[0].1.get(&j).expect("query missing") {
                        return false;
                    }
                } else {
                    if v != self.final_value.unwrap() {
                        return false;
                    }
                }
            }
            shift_inv *= shift_inv;
            generator_inv *= generator_inv;
            len >>= 1;
        }
        true
    }

    fn set_final_value(&mut self, final_value: T) {
        self.final_value = Some(final_value);
    }
}

impl<T: Field> FriProver<T> {
    pub fn new(
        log_poly_degree: usize,
        interpolate_value: Vec<T>,
        coset: &Coset<T>,
    ) -> FriProver<T> {
        FriProver {
            log_poly_degree,
            coset: coset.clone(),
            interpolate_values: vec![interpolate_value],
            merkle_tree: vec![],
            verifier: None,
        }
    }

    pub fn from_polynomial(
        log_poly_degree: usize,
        polynomial: Polynomial<T>,
        coset: &Coset<T>,
    ) -> FriProver<T> {
        assert_eq!(1 << log_poly_degree, polynomial.degree() + 1);
        let interpolate_values = vec![polynomial.evaluation_over_coset(coset)];
        FriProver {
            log_poly_degree,
            coset: coset.clone(),
            interpolate_values,
            merkle_tree: vec![],
            verifier: None,
        }
    }

    pub fn set_verifier(&mut self, verifier: &Rc<RefCell<FriVerifier<T>>>) {
        self.verifier = Some(verifier.clone());
    }

    fn batch_inverse_and_mul(vec: Vec<T>, k: T) -> Vec<T> {
        let mut res = Vec::with_capacity(vec.len());
        let mut c = vec[0];
        res.push(c);
        for i in 1..vec.len() {
            c *= vec[i];
            res.push(c);
        }
        let mut c_inv = c.inverse() * k;
        for i in (1..vec.len()).rev() {
            res[i] = res[i - 1] * c_inv;
            c_inv *= vec[i];
        }
        res[0] = c_inv;
        res
    }

    fn evaluation_next_domain(
        interpolate_value: &Vec<T>,
        current_domain: &Coset<T>,
        challenge: T,
    ) -> Vec<T> {
        let mut res = vec![];
        assert_eq!(interpolate_value.len(), current_domain.num_elements());
        let inv_2 = T::from_int(2).inverse();
        let mut shift_inv = current_domain.shift().inverse();
        let generator_inv = current_domain.generator().inverse();
        for i in 0..(interpolate_value.len() / 2) {
            let x = interpolate_value[i];
            let nx = interpolate_value[i + interpolate_value.len() / 2];
            let new_v = (x + nx) + challenge * (x - nx) * shift_inv;
            res.push(new_v * inv_2);
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
        let mut domain_size = self.coset.num_elements();
        let mut domain = self.coset.clone();
        let mut shift = domain.shift();
        let verifier = self.verifier.clone().unwrap();
        for _i in 0..self.log_poly_degree {
            let merkle_tree_prover =
                Self::merkle_tree_commit(self.interpolate_values.last().unwrap());
            let commit = merkle_tree_prover.commit();
            verifier.borrow_mut().receive_root(domain_size / 2, &commit);
            self.merkle_tree.push(merkle_tree_prover);

            let challenge = verifier.borrow_mut().get_challenge();
            let next_evalutation = Self::evaluation_next_domain(
                &self.interpolate_values.last().unwrap(),
                &domain,
                challenge,
            );
            self.interpolate_values.push(next_evalutation);
            shift *= shift;

            domain_size >>= 1;
            domain = Coset::new(domain_size, shift);
        }

        verifier
            .borrow_mut()
            .set_final_value(self.interpolate_values.last().unwrap()[0]);
    }

    pub fn query(&self, points: &Vec<usize>) -> Vec<(Vec<u8>, HashMap<usize, T>)> {
        let mut res = vec![];
        let mut leaf_indices = points.clone();
        for i in 0..self.log_poly_degree {
            let len = self.interpolate_values[i].len();

            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            let mut values = HashMap::new();
            for j in &leaf_indices {
                values.insert(*j, self.interpolate_values[i][*j]);
                values.insert(j + len / 2, self.interpolate_values[i][*j + len / 2]);
            }
            let proof_bytes = self.merkle_tree[i].open(&leaf_indices);
            res.push((proof_bytes, values));
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::mersenne61_ext::Mersenne61Ext;
    use rand::Rng;

    #[test]
    fn batch_inverse_and_mul() {
        let mut vec = vec![];
        for _i in 0..64 {
            vec.push(Mersenne61Ext::random_element());
        }
        let k = Mersenne61Ext::random_element();
        let mut res1 = vec![];
        for i in &vec {
            res1.push(i.inverse() * k);
        }
        let res2 = FriProver::batch_inverse_and_mul(vec, k);
        assert_eq!(res1, res2);
    }

    #[test]
    fn low_degree_test() {
        let shift = Mersenne61Ext::random_element();
        let domain = Coset::new(1 << 10, shift);
        let poly_degree_bound = 1 << 8;
        let poly = Polynomial::random_polynomial(poly_degree_bound);
        let verifier = Rc::new(RefCell::new(FriVerifier::new(&domain, 8)));
        let prover = Rc::new(RefCell::new(FriProver::from_polynomial(8, poly, &domain)));
        verifier.borrow_mut().set_prover(&prover);
        prover.borrow_mut().set_verifier(&verifier);
        prover.borrow_mut().prove();
        let mut points = vec![];
        for _i in 0..10 {
            points.push(rand::thread_rng().gen_range(0..domain.num_elements()));
        }
        let evaluation = prover.borrow().query(&points);
        assert!(verifier.borrow().verify(points, evaluation));
    }
}
