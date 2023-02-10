use crate::fields::Field;

pub struct Radix2Domain<T: Field> {
    order: u64,
    omega: T
}

impl<T: Field> Radix2Domain<T> {
    pub fn new(order: u64) -> Self {
        Radix2Domain { 
            order, 
            omega: T::get_generator(order)
        }
    }

    pub fn order(&self) -> u64 {
        self.order
    }

    pub fn omega(&self) -> T {
        self.omega
    }

    pub fn fft(&self, a: &mut Vec<T>) {
        assert_eq!(a.len() as u64, self.order);
        _fft(a, self.omega);
    }

    pub fn ifft(&self, a: &mut Vec<T>) {
        assert_eq!(a.len() as u64, self.order);
        _fft(a, self.omega.inverse());
        let t = T::from_int(self.order).inverse();
        for i in a {
            *i *= t;
        }
    }

    pub fn coset_fft(&self, a: &mut Vec<T>, shift: T) {
        multiply_by_coset(a, shift);
        self.fft(a);
    }

    pub fn coset_ifft(&self, a: &mut Vec<T>, shift: T) {
        self.ifft(a);
        multiply_by_coset(a, shift.inverse());
    }
}

fn multiply_by_coset<T: Field>(a: &mut Vec<T>, shift: T) {
    let mut t = shift;
    for i in 1..a.len() {
        a[i] *= t;
        t *= shift;
    }
}

fn bitreverse(mut x: usize, len: usize) -> usize {
    let mut r = 0;
    for _i in 0..len {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

fn _fft<T: Field>(a: &mut Vec<T>, omega: T) {
    let n = a.len();
    let log_n = (n as f64).log2() as usize;
    assert_eq!(1 << log_n, n);
    for i in 0..n {
        let rank = bitreverse(i, log_n);
        if i < rank {
            (a[i], a[rank]) = (a[rank], a[i]);
        }
    }
    let mut m = 1usize;
    for _i in 0..log_n {
        let w_m = omega.pow((n / (m * 2)) as u64);
        for j in (0..n).step_by(2*m) {
            let mut w = T::from_int(1);
            for k in 0..m {
                let t = w * a[j + k + m];
                a[j + k + m] = a[j + k] - t;
                a[j + k] += t;
                w *= w_m;
            }
        }
        m *= 2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::fields::fp64::Fp64;
    #[test]
    fn fft_and_ifft() {
        let mut a: Vec<Fp64> = vec![];
        let mut b: Vec<Fp64> = vec![];
        for _i in 0..16 {
            a.push(Fp64::random_element());
            b.push(Fp64::random_element());
        }
        for _i in 16..32 {
            a.push(Fp64::from_int(0));
            b.push(Fp64::from_int(0));
        }
        let mut a_times_b: Vec<Fp64> = vec![];
        for _i in 0..32 {
            a_times_b.push(Fp64::from_int(0));
        }
        for i in 0..16usize {
            for j in 0..16usize {
                a_times_b[i + j] += a[i] * b[j];
            }
        }
        let domain: Radix2Domain<Fp64> = Radix2Domain::new(32);
        domain.fft(&mut a);
        domain.fft(&mut b);
        for i in 0..a.len() {
            a[i] *= b[i];
        }
        domain.ifft(&mut a);
        assert_eq!(&a, &a_times_b);
        let shift = Fp64::random_element();
        domain.coset_fft(&mut a, shift);
        domain.coset_ifft(&mut a, shift);
        assert_eq!(&a, &a_times_b);
    }
}
