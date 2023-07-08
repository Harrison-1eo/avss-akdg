use crate::random_oracle::RandomOracle;
use crate::util::QueryResult;
use crate::{
    algebra::{coset::Coset, field::Field},
    one2many::verifier::One2ManyVerifier,
};
use std::{cell::RefCell, rc::Rc};

pub struct AvssParty<T: Field, const N: usize> {
    pub verifier: Rc<RefCell<One2ManyVerifier<T, N>>>,
    open_point: Vec<T>,
    share: Option<T>,
}

impl<T: Field + 'static, const N: usize> AvssParty<T, N> {
    pub fn share(&self) -> T {
        self.share.unwrap()
    }

    pub fn set_share(&mut self, share: T) {
        self.share = Some(share);
    }

    pub fn open_point(&self) -> &Vec<T> {
        &self.open_point
    }

    pub fn new(
        interpolate_coset: &Coset<T>,
        open_point: Vec<T>,
        oracle: &Rc<RefCell<RandomOracle<T>>>,
    ) -> AvssParty<T, N> {
        AvssParty {
            verifier: Rc::new(RefCell::new(One2ManyVerifier::new_with_default_map(
                interpolate_coset,
                oracle,
            ))),
            open_point,
            share: None,
        }
    }

    pub fn verify(
        &self,
        folding_proofs: Vec<QueryResult<T>>,
        function_proofs: Vec<QueryResult<T>>,
    ) -> bool {
        self.verifier.borrow().verify_with_extra_folding(
            folding_proofs,
            function_proofs,
            &self.open_point,
            self.share.unwrap(),
        )
    }
}
