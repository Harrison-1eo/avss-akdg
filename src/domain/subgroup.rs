use std::{rc::Rc, cell::RefCell};

use crate::fields::Field;
use super::radix2_domain::Radix2Domain;

pub struct Coset<T: Field> {
    elements: Rc<RefCell<Vec<T>>>,
    fft_cache: Rc<RefCell<Vec<T>>>,
    fft_eval_domain: Radix2Domain<T>,
    shift: T
}

impl<T: Field> Coset<T> {
    pub fn new(order: u64, shift: T) -> Self {
        assert!(!shift.is_zero());
        Coset { 
            elements: Rc::new(RefCell::new(vec![])), 
            fft_cache: Rc::new(RefCell::new(vec![])), 
            fft_eval_domain: Radix2Domain::new(order),
            shift
        }
    }

    pub fn all_elements(&self) -> Vec<T> {
        let mut elements = self.elements.borrow_mut();
        if elements.len() == 0 {
            let mut el = self.shift;
            for _i in 0..self.fft_eval_domain.order() {
                elements.push(el);
                el *= self.fft_eval_domain.omega();
            }
        }
        elements.clone()
    }

    pub fn num_elements(&self) -> usize {
        self.fft_eval_domain.order() as usize
    }

    pub fn fft(&self, evals: &Vec<T>) -> Vec<T> {
        let mut a = evals.clone();
        self.fft_eval_domain.coset_fft(&mut a, self.shift);
        a
    }

    pub fn ifft(&self, evals: &Vec<T>) -> Vec<T> {
        if evals.len() == 1 {
            return vec![evals[0]];
        };
        let mut a = evals.clone();
        self.fft_eval_domain.coset_ifft(&mut a, self.shift);
        a
    }

    pub fn shift(&self) -> T {
        self.shift
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::fp64::Fp64;

    #[test]
    fn all_elements() {
        let r = Fp64::random_element();
        let coset = Coset::new(32, r);
        let elements = coset.all_elements();
        assert_eq!(elements[0], r);
        let omega = coset.fft_eval_domain.omega();
        for i in 0..elements.len() - 1 {
            assert_eq!(elements[i] * omega, elements[i + 1]);
        }
        assert_eq!(*elements.last().unwrap() * omega, elements[0]);
    }
}
