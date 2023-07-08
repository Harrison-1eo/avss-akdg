// #[macro_use]
// extern crate criterion;

// use criterion::*;

// use std::{cell::RefCell, rc::Rc};

// use crate::algebra::field::Field;
// use algebra::polynomial::MultilinearPolynomial;
// use avss::dealer::Dealer;
// use avss::party::AvssParty;
// use one2many::RandomOracle;

// use crate::algebra::coset::Coset;
// use crate::algebra::field::mersenne61_ext::Mersenne61Ext;

// fn bench(c: &mut Criterion) {
//     let oracle = Rc::new(RefCell::new(RandomOracle::new()));
//     let interpolate_coset = Coset::new(1 << 11, Mersenne61Ext::random_element());
//     let polynomial = MultilinearPolynomial::random_polynomial(8);

//     let x_shift = Mersenne61Ext::random_element();
//     let coset_x = Coset::new(1 << 5, x_shift);
//     let mut folding_parameter = vec![];
//     let v = vec![4, 6, 2, 1];
//     for i in &v {
//         folding_parameter.push(coset_x.pow(*i).all_elements());
//     }
//     let y_shift = Mersenne61Ext::random_element();
//     let coset_y = Coset::new(1 << 5, y_shift);
//     let last_len = folding_parameter.last().unwrap().len();
//     for i in v {
//         folding_parameter.push(
//             coset_y
//                 .pow(i)
//                 .all_elements()
//                 .iter()
//                 .map(|y| {
//                     (0..last_len)
//                         .into_iter()
//                         .map(|_| *y)
//                         .collect::<Vec<Mersenne61Ext>>()
//                 })
//                 .flatten()
//                 .collect(),
//         );
//     }
//     let mut parties: Vec<AvssParty<Mersenne61Ext, 8>> = vec![];
//     assert_eq!(folding_parameter.last().unwrap().len(), 1024);
//     for i in 0..1024 {
//         let mut open_point = vec![];
//         for j in 0..8 {
//             open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
//         }
//         let x = coset_x.all_elements()[i % 32];
//         assert_eq!(x, open_point[3]);
//         assert_eq!(x.pow(2), open_point[2]);
//         assert_eq!(x.pow(6), open_point[1]);
//         assert_eq!(x.pow(4), open_point[0]);
//         let y = coset_y.all_elements()[i / 32];
//         assert_eq!(y, open_point[7]);
//         assert_eq!(y.pow(2), open_point[6]);
//         assert_eq!(y.pow(6), open_point[5]);
//         assert_eq!(y.pow(4), open_point[4]);
//         parties.push(AvssParty::new(&interpolate_coset, open_point, &oracle));
//     }

//     c.bench_function("with_setup", move |b| {
//         b.iter_batched_ref(
//             || data.clone(),
//             |mut data| sort(&mut data),
//             BatchSize::SmallInput,
//         )
//     });
// }

// criterion_group!(benches, bench);
// criterion_main!(benches);
