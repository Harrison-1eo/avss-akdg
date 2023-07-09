use crate::algebra::field::Field;
use rand::Rng;

pub struct RandomOracle<T: Field> {
    folding_challenges: Vec<T>,
    usize_elements: Option<Vec<usize>>,
}

impl<T: Field> RandomOracle<T> {
    pub fn new() -> Self {
        RandomOracle {
            folding_challenges: vec![],
            usize_elements: None,
        }
    }

    pub fn clear(&mut self) {
        self.folding_challenges.clear();
        self.usize_elements = None
    }

    pub fn query_list(&self) -> Vec<usize> {
        self.usize_elements.clone().unwrap()
    }

    pub fn generate_queries(&mut self, len: usize) {
        self.usize_elements = Some(
            (0..len)
                .into_iter()
                .map(|_| rand::thread_rng().gen())
                .collect(),
        )
    }

    pub fn get_challenge(&self, index: usize) -> T {
        self.folding_challenges[index]
    }

    pub fn generate_challenge(&mut self) -> T {
        let challenge = T::random_element();
        self.folding_challenges.push(challenge);
        challenge
    }
}
