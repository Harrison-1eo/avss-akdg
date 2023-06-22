use std::{cell::RefCell, collections::HashMap, rc::Rc};

use crate::algebra::{coset::Coset, field::Field, polynomial::*};

pub struct FriProver<T: Field> {
    parameter_array: Vec<usize>,
    coset: Coset<T>,
    interpolate_values: Vec<Vec<T>>,
    verifier: Option<Rc<RefCell<FriVerifier<T>>>>,
}

pub struct FriVerifier<T: Field> {
    parameter_array: Vec<usize>,
    coset: Coset<T>,
    poly_degree_bound: usize,
    challenges: Vec<T>,
    prover: Option<Rc<RefCell<FriProver<T>>>>,
    final_poly: Option<Polynomial<T>>,
}

impl<T: Field> FriVerifier<T> {
    pub fn new(
        parameter_array: &Vec<usize>,
        coset: &Coset<T>,
        poly_degree_bound: usize,
    ) -> FriVerifier<T> {
        FriVerifier {
            parameter_array: parameter_array.clone(),
            coset: coset.clone(),
            poly_degree_bound,
            challenges: vec![],
            prover: None,
            final_poly: None,
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
        let mut shift = self.coset.shift();
        let mut generator = self.coset.generator();
        let mut len = self.coset.num_elements();
        for i in 0..self.parameter_array.len() {
            let eta = self.parameter_array[i];
            points = points.iter_mut().map(|v| *v % (len >> eta)).collect();
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
                let mut tmp = vec![];
                for k in (points[j]..len).step_by(len >> eta) {
                    tmp.push(evaluation[i].get(&k).unwrap().clone());
                }
                let coset = Coset::new(1 << eta, shift * generator.pow(points[j] as u64));
                let poly = Polynomial::new(coset.ifft(&tmp));
                let v = poly.evaluation_at(self.challenges[i]);
                if i < self.parameter_array.len() - 1 {
                    if v != *evaluation[i + 1].get(&points[j]).expect("query missing") {
                        return false;
                    }
                } else {
                    let x = shift.pow(1 << eta as u64) * generator.pow((points[j] << eta) as u64);
                    let poly_v = self.final_poly.as_ref().unwrap().evaluation_at(x);
                    if v != poly_v {
                        return false;
                    }
                }
            }
            for _j in 0..eta {
                shift *= shift;
                generator *= generator;
            }
            len >>= eta;
        }
        true
    }

    fn set_final_poly(&mut self, final_poly: Polynomial<T>) {
        let mut d = self.poly_degree_bound;
        for i in &self.parameter_array {
            d >>= *i;
        }
        if final_poly.degree() > d {
            panic!("invalid final polynomial");
        }
        self.final_poly = Some(final_poly);
    }
}

impl<T: Field> FriProver<T> {
    pub fn new(
        interpolate_value: Vec<T>,
        parameter_array: &Vec<usize>,
        coset: &Coset<T>,
    ) -> FriProver<T> {
        FriProver {
            parameter_array: parameter_array.clone(),
            coset: coset.clone(),
            interpolate_values: vec![interpolate_value],
            verifier: None,
        }
    }

    pub fn from_polynomial(
        polynomial: Polynomial<T>,
        parameter_array: &Vec<usize>,
        coset: &Coset<T>,
    ) -> FriProver<T> {
        let interpolate_values = vec![polynomial.evaluation_over_coset(coset)];
        FriProver {
            parameter_array: parameter_array.clone(),
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
        eta: usize,
        challenge: T,
    ) -> Vec<T> {
        let coset_size = 1 << eta;
        let num_coset = current_domain.num_elements() / coset_size;
        let h_inc = current_domain.generator();
        let h_inc_to_coset_inv_plus_one = h_inc.pow(coset_size as u64).inverse() * h_inc;
        let shiftless_coset = Coset::new(coset_size, T::from_int(1));
        let g = shiftless_coset.generator();
        let g_inv = g.inverse();
        let x_to_ordet_coset = challenge.pow(coset_size as u64);
        let mut shifted_x_elements = Vec::with_capacity(coset_size);
        shifted_x_elements.push(challenge);
        for i in 1..coset_size {
            shifted_x_elements.push(shifted_x_elements[i - 1] * g_inv);
        }
        let mut cur_h = current_domain.shift();
        let first_h_to_coset_inv_plus_one = cur_h.pow(coset_size as u64).inverse() * cur_h;
        let mut cur_coset_constant_plus_h = x_to_ordet_coset * first_h_to_coset_inv_plus_one;
        let mut elements_to_invert = Vec::with_capacity(interpolate_value.len());
        let mut constant_for_each_coset = Vec::with_capacity(num_coset);

        let constant_for_all_coset = T::from_int(coset_size as u64).inverse();
        for _i in 0..num_coset {
            let coset_contant = cur_coset_constant_plus_h - cur_h;
            constant_for_each_coset.push(coset_contant);
            assert_ne!(coset_contant, T::from_int(0));
            for k in 0..coset_size {
                elements_to_invert.push(shifted_x_elements[k] - cur_h);
            }
            cur_h *= h_inc;
            cur_coset_constant_plus_h *= h_inc_to_coset_inv_plus_one;
        }
        let lagrange_coefficients =
            Self::batch_inverse_and_mul(elements_to_invert, constant_for_all_coset);
        let mut next_interpolate_value = Vec::with_capacity(num_coset);
        for j in 0..num_coset {
            let mut interpolation = T::from_int(0);
            for k in 0..coset_size {
                interpolation += interpolate_value[k * num_coset + j]
                    * lagrange_coefficients[j * coset_size + k];
            }
            interpolation *= constant_for_each_coset[j];
            next_interpolate_value.push(interpolation);
        }
        next_interpolate_value
    }

    pub fn prove(&mut self) {
        let mut domain_size = self.coset.num_elements();
        let mut domain = self.coset.clone();
        let mut shift = domain.shift();
        let verifier = self.verifier.clone().unwrap();
        for i in 0..self.parameter_array.len() {
            let eta = self.parameter_array[i];
            // todo: merkle tree commit
            let challenge = verifier.borrow_mut().get_challenge();
            let next_evalutation = Self::evaluation_next_domain(
                &self.interpolate_values.last().unwrap(),
                &domain,
                eta,
                challenge,
            );
            self.interpolate_values.push(next_evalutation);
            for _j in 0..eta {
                shift *= shift;
            }
            domain_size >>= eta;
            domain = Coset::new(domain_size, shift);
        }
        let final_poly = Polynomial::new(domain.ifft(self.interpolate_values.last().unwrap()));
        verifier.borrow_mut().set_final_poly(final_poly);
    }

    pub fn query(&self, points: &Vec<usize>) -> Vec<HashMap<usize, T>> {
        let mut res = vec![];
        let mut points = points.clone();
        for i in 0..self.parameter_array.len() {
            res.push(HashMap::new());
            let len = self.interpolate_values[i].len();
            let eta = self.parameter_array[i];

            points = points.iter_mut().map(|v| *v % (len >> eta)).collect();
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
                for k in (*j..len).step_by(len >> eta) {
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
    fn next_evalutation() {
        let polynomial = Polynomial::random_polynomial(64);
        let domain = Coset::new(512, Mersenne61Ext::random_element());
        let interpolate_value = domain.fft(polynomial.coefficients());
        let challenge = Mersenne61Ext::random_element();
        let eta = 4;
        let res1 = FriProver::evaluation_next_domain(&interpolate_value, &domain, eta, challenge);

        let num_coset = interpolate_value.len() / (1 << eta);
        let mut res2 = Vec::with_capacity(num_coset);
        let mut shift = domain.shift();
        for i in 0..num_coset {
            let mut value = vec![];
            for j in 0..(1 << eta) {
                value.push(interpolate_value[j * num_coset + i]);
            }
            let coset = Coset::new(1 << eta, shift);
            let poly = Polynomial::new(coset.ifft(&value));
            res2.push(poly.evaluation_at(challenge));
            shift *= domain.generator();
        }
        assert_eq!(res1, res2);
    }

    #[test]
    fn low_degree_test() {
        let shift = Mersenne61Ext::random_element();
        let domain = Coset::new(1 << 10, shift);
        let poly_degree_bound = 1 << 8;
        let poly = Polynomial::random_polynomial(poly_degree_bound);
        let parameter_array = vec![2, 3, 1, 2];
        let verifier = Rc::new(RefCell::new(FriVerifier::new(
            &parameter_array,
            &domain,
            poly_degree_bound,
        )));
        let prover = Rc::new(RefCell::new(FriProver::from_polynomial(
            poly,
            &parameter_array,
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
