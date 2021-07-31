use crate::seq::*;
use num_complex::Complex;
use num_traits::Zero;
use rand::prelude::*;

pub struct FuncSeq {
    n: usize,
    f: Vec<Vec<f64>>,
}

pub fn calc_prob(e: f64, ep: f64, t: f64) -> f64 {
    let exp = -(ep - e) / t;
    exp.exp()
}

impl FuncSeq {
    pub fn new(n: usize, f: Vec<Vec<f64>>) -> Self {
        Self { n, f }
    }

    pub fn approx_k(&self, step: usize, origin: Sequence, s0: Vec<f64>, r: Vec<f64>) -> Vec<f64> {
        assert_eq!(self.f.len() > 0, true);
        // Implement simulated anealing
        let mut s = s0;
        let mut rng = thread_rng();
        let calc = |x: Vec<f64>| -> f64 {
            let mut data: Vec<Complex<f64>> = vec![Complex::zero(); self.n];
            for i in 0..self.n {
                data[i] = Complex::from(self.f[0][i].powf(x[0]));
                for j in 1..self.f.len() {
                    data[i] *= Complex::from(self.f[j][i].powf(x[j]));
                }
            }
            let q = Sequence::from(data);
            origin.rem(q).size()
        };
        let mut e = calc(s.clone());
        for k in 0..step {
            let t = 1.0 - ((k as f64) + 1.0) / (step as f64);
            let mut s_new: Vec<f64> = vec![];
            for i in 0..self.f.len() {
                s_new.push(s[i] + rng.gen_range(-r[i]..=r[i]));
            }
            let ep = calc(s_new.clone());
            if calc_prob(e, ep, t) >= rng.gen_range(0.0..=1.0) {
                s = s_new;
                e = ep;
            }
        }
        s
    }
}
