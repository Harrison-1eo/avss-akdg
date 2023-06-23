use std::{cell::RefCell, collections::HashMap, rc::Rc};

use crate::algebra::{coset::Coset, field::Field, polynomial::*};

pub struct FriProver<T: Field> {
    log_poly_degree: usize,
    coset: Coset<T>,
    interpolate_values: Vec<Vec<T>>,
    verifier: Option<Rc<RefCell<FriVerifier<T>>>>,
}

pub struct FriVerifier<T: Field> {
    coset: Coset<T>,
    log_poly_degree: usize,
    challenges: Vec<T>,
    prover: Option<Rc<RefCell<FriProver<T>>>>,
    final_value: Option<T>,
}

impl<T: Field> FriVerifier<T> {
    pub fn new(
        coset: &Coset<T>,
        log_poly_degree: usize,
    ) -> FriVerifier<T> {
        FriVerifier {
            coset: coset.clone(),
            log_poly_degree,
            challenges: vec![],
            prover: None,
            final_value: None,
        }
    }

    fn get_challenge(&mut self) -> T {
        let challenge = T::random_element();
        self.challenges.push(challenge);
        challenge
    }

    pub fn set_prover(&mut self, prover: &Rc<RefCell<FriProver<T>>>) {
        self.prover = Some(prover.clone());
    }

    pub fn verify(&self, mut points: Vec<usize>, evaluation: Vec<HashMap<usize, T>>) -> bool {
        let mut shift_inv = self.coset.shift().inverse();
        let mut generator_inv = self.coset.generator().inverse();
        let mut len = self.coset.num_elements();
        for i in 0..self.log_poly_degree {
            points = points.iter_mut().map(|v| *v % (len >> 1)).collect();
            points.sort();
            let mut query_list = vec![];
            query_list.push(points[0]);
            for j in 1..points.len() {
                if points[j] != points[j - 1] {
                    query_list.push(points[j]);
                }
            }
            points = query_list;

            for j in 0..points.len() {
                let x = evaluation[i].get(&points[j]).unwrap().clone();
                let nx = evaluation[i].get(&(points[j] + len / 2)).unwrap().clone();
                let v = x + nx + self.challenges[i] * (x - nx) * shift_inv * generator_inv.pow(points[j] as u64);
                let v = v * T::from_int(2).inverse();
                if i < self.log_poly_degree - 1 {
                    if v != *evaluation[i + 1].get(&points[j]).expect("query missing") {
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
        challenge: T
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

    pub fn prove(&mut self) {
        let mut domain_size = self.coset.num_elements();
        let mut domain = self.coset.clone();
        let mut shift = domain.shift();
        let verifier = self.verifier.clone().unwrap();
        for _i in 0..self.log_poly_degree {
            // todo: merkle tree commit
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

        verifier.borrow_mut().set_final_value(self.interpolate_values.last().unwrap()[0]);
    }

    pub fn query(&self, points: &Vec<usize>) -> Vec<HashMap<usize, T>> {
        let mut res = vec![];
        let mut points = points.clone();
        for i in 0..self.log_poly_degree {
            res.push(HashMap::new());
            let len = self.interpolate_values[i].len();

            points = points.iter_mut().map(|v| *v % (len >> 1)).collect();
            points.sort();
            let mut query_list = vec![];
            query_list.push(points[0]);
            for j in 1..points.len() {
                if points[j] != points[j - 1] {
                    query_list.push(points[j]);
                }
            }
            points = query_list;
            for j in &points {
                for k in (*j..len).step_by(len >> 1) {
                    res[i].insert(k, self.interpolate_values[i][k]);
                }
            }
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
        let verifier = Rc::new(RefCell::new(FriVerifier::new(
            &domain,
            8,
        )));
        let prover = Rc::new(RefCell::new(FriProver::from_polynomial(
            8,
            poly,
            &domain,
        )));
        verifier.borrow_mut().set_prover(&prover);
        prover.borrow_mut().set_verifier(&verifier);
        prover.borrow_mut().prove();
        let mut points = vec![];
        for _i in 0..1 {
            points.push(rand::thread_rng().gen_range(0..domain.num_elements()));
        }
        let evaluation = prover.borrow().query(&points);
        assert!(verifier.borrow().verify(points, evaluation));
    }
}
