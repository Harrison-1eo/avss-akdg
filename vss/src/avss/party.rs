use crate::one2many::verifier::One2ManyVerifier;
use std::{cell::RefCell, rc::Rc};
use util::algebra::{coset::Coset, field::Field, polynomial::MultilinearPolynomial};
use util::query_result::QueryResult;
use util::random_oracle::RandomOracle;

#[derive(Clone)]
pub struct AvssParty<T: Field> {
    pub verifier: Rc<RefCell<One2ManyVerifier<T>>>,
    open_point: Vec<T>,
    final_poly: Option<MultilinearPolynomial<T>>,
}

impl<T: Field + 'static> AvssParty<T> {
    pub fn share(&self) -> T {
        let poly = self.final_poly.as_ref().unwrap();
        let variable_num = poly.variable_num();
        let n = self.open_point.len();
        poly.evaluate(&self.open_point[n - variable_num..].to_vec())
    }

    pub fn set_share(&mut self, final_poly: &MultilinearPolynomial<T>) {
        self.final_poly = Some(final_poly.clone());
    }

    pub fn open_point(&self) -> &Vec<T> {
        &self.open_point
    }

    pub fn new(
        total_round: usize,
        interpolate_coset: &Vec<Coset<T>>,
        open_point: Vec<T>,
        oracle: &RandomOracle<T>,
    ) -> AvssParty<T> {
        AvssParty {
            verifier: Rc::new(RefCell::new(One2ManyVerifier::new_with_default_map(
                total_round,
                open_point.len(),
                interpolate_coset,
                oracle,
            ))),
            open_point,
            final_poly: None,
        }
    }

    pub fn verify(
        &self,
        folding_proofs: &Vec<QueryResult<T>>,
        function_proofs: &Vec<QueryResult<T>>,
    ) -> bool {
        self.verifier.borrow().verify_with_extra_folding(
            folding_proofs,
            function_proofs,
            &self.open_point,
            self.final_poly.as_ref().unwrap(),
        )
    }
}
