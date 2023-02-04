use std::{rc::Rc, cell::RefCell};

use crate::fields::Field;
use super::radix2_domain::Radix2Domain;

struct Subgroup<T: Field> {
    elements: Rc<RefCell<Vec<T>>>,
    fft_cache: Rc<RefCell<Vec<T>>>,
    fft_eval_domain: Radix2Domain<T>
}

impl<T: Field> Subgroup<T> {
    fn new(order: u64) -> Self {
        Subgroup { 
            elements: Rc::new(RefCell::new(vec![])), 
            fft_cache: Rc::new(RefCell::new(vec![])), 
            fft_eval_domain: Radix2Domain::new(order)
        }
    }
}

struct Coset<T: Field> {
    subgroup: Subgroup<T>,
    shift: T
}

impl<T: Field> Coset<T> {
    pub fn new(order: u64, shift: T) -> Self {
        Coset { 
            subgroup: Subgroup::new(order), 
            shift
        }
    }

    pub fn fft(&self, evals: &Vec<T>) -> Vec<T> {
        todo!("to optimize");
        let mut a = evals.clone();
        self.subgroup.fft_eval_domain.coset_fft(&mut a, self.shift);
        a
    }

    pub fn ifft(&self, evals: &Vec<T>) -> Vec<T> {
        if evals.len() == 1 {
            return vec![evals[0]];
        };
        let mut a = evals.clone();
        self.subgroup.fft_eval_domain.coset_ifft(&mut a, self.shift);
        a
    }
}
