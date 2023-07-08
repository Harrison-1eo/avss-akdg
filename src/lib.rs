pub mod algebra {
    pub mod coset;
    pub mod field;
    pub mod polynomial;
}
pub mod merkle_tree;
pub mod one2many;
pub mod util;
pub mod avss {
    pub mod dealer;
    pub mod party;
}
pub mod random_oracle;

#[cfg(test)]
mod tests {
    use std::{cell::RefCell, rc::Rc};

    use crate::algebra::field::Field;
    use crate::algebra::polynomial::MultilinearPolynomial;
    use crate::avss::dealer::Dealer;
    use crate::avss::party::AvssParty;
    use crate::random_oracle::RandomOracle;

    use crate::algebra::coset::Coset;
    use crate::algebra::field::mersenne61_ext::Mersenne61Ext;
    #[test]
    fn test_avss() {
        let oracle = Rc::new(RefCell::new(RandomOracle::new()));
        let interpolate_coset = Coset::new(1 << 11, Mersenne61Ext::random_element());
        let polynomial = MultilinearPolynomial::random_polynomial(8);

        let x_shift = Mersenne61Ext::random_element();
        let coset_x = Coset::new(1 << 5, x_shift);
        let mut folding_parameter = vec![];
        let v = vec![4, 6, 2, 1];
        for i in &v {
            folding_parameter.push(coset_x.pow(*i).all_elements());
        }
        let y_shift = Mersenne61Ext::random_element();
        let coset_y = Coset::new(1 << 5, y_shift);
        let last_len = folding_parameter.last().unwrap().len();
        for i in v {
            folding_parameter.push(
                coset_y
                    .pow(i)
                    .all_elements()
                    .iter()
                    .map(|y| {
                        (0..last_len)
                            .into_iter()
                            .map(|_| *y)
                            .collect::<Vec<Mersenne61Ext>>()
                    })
                    .flatten()
                    .collect(),
            );
        }
        let mut parties: Vec<AvssParty<Mersenne61Ext, 8>> = vec![];
        assert_eq!(folding_parameter.last().unwrap().len(), 1024);
        for i in 0..1024 {
            let mut open_point = vec![];
            for j in 0..8 {
                open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
            }
            let x = coset_x.all_elements()[i % 32];
            assert_eq!(x, open_point[3]);
            assert_eq!(x.pow(2), open_point[2]);
            assert_eq!(x.pow(6), open_point[1]);
            assert_eq!(x.pow(4), open_point[0]);
            let y = coset_y.all_elements()[i / 32];
            assert_eq!(y, open_point[7]);
            assert_eq!(y.pow(2), open_point[6]);
            assert_eq!(y.pow(6), open_point[5]);
            assert_eq!(y.pow(4), open_point[4]);
            parties.push(AvssParty::new(&interpolate_coset, open_point, &oracle));
        }
        let mut dealer = Dealer::<Mersenne61Ext, 8>::new(
            &polynomial,
            &interpolate_coset,
            &oracle,
            &folding_parameter,
        );
        dealer.send_evaluations(&mut parties);
        dealer.commit_functions(&parties);
        dealer.prove();
        dealer.commit_foldings(&parties);
        oracle.borrow_mut().generate_queries(20);
        let (folding, function) = dealer.query();
        let mut folding715 = vec![];
        let mut function715 = vec![];
        for i in 0..8 {
            if i < 7 {
                folding715.push(folding[i][715 % folding[i].len()].clone());
            }
            function715.push(function[i][715 % function[i].len()].clone());
        }
        assert!(parties[715].verify(folding715, function715));
        assert_eq!(
            parties[715].share(),
            polynomial.evaluate(parties[715].open_point())
        );
    }
}
