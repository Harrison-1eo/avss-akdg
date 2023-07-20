use rand::Rng;
use crate::algebra::field::Field;

#[derive(Debug, Clone)]
pub struct RandomOracle<T: Field> {
    pub folding_challenges: Vec<T>,
    pub query_list: Vec<usize>,
}

impl<T: Field> RandomOracle<T> {
    pub fn new(total_round: usize, query_num: usize) -> Self {
        RandomOracle {
            folding_challenges: (0..total_round).into_iter().map(|_| T::random_element()).collect(),
            query_list: (0..query_num)
                .into_iter()
                .map(|_| rand::thread_rng().gen())
                .collect(),
        }
    }
}
