use crate::funcseq::*;
use crate::seq::*;
use num_complex::Complex;
use num_traits::Zero;
use rand::prelude::*;

const ARR_SIZE: usize = 8;
const STEP: usize = 10000;
const TESTCASES: usize = 100;
const S0: f64 = 0.0;
const R0: f64 = 1.0;
const DELTA: f64 = 0.3;

#[test]
fn log_approx_test() {
    let mut rng = thread_rng();
    for _case in 0..TESTCASES {
        let mut v = vec![Complex::zero(); ARR_SIZE];
        let c = rng.gen_range(0.0..=10000.0);
        let e = rng.gen_range(0.0..=20.0);
        for i in 0..ARR_SIZE {
            v[i] = Complex::from(c * (1.0 + (i as f64)).log10().powf(e));
        }
        let seq = Sequence::from(v);
        let mut log = vec![0.0; seq.n];
        for i in 0..seq.n {
            log[i] = (1.0 + (i as f64)).log10();
        }
        let func = FuncSeq::new(seq.n, vec![log]);
        let exp = func.approx_k(STEP, seq, vec![S0], vec![R0]);
        let abs = (e - exp[0]).abs();
        assert_eq!(abs <= DELTA, true);
    }
}

#[test]
fn twos_approx_test() {
    let mut rng = thread_rng();
    for _case in 0..TESTCASES {
        let mut v = vec![Complex::zero(); ARR_SIZE];
        let c = rng.gen_range(0.0..=1000.0);
        let e = rng.gen_range(0.0..=20.0);
        for i in 0..ARR_SIZE {
            v[i] = Complex::from(c * (2.0f64).powi(i as i32).powf(e));
        }
        let seq = Sequence::from(v);
        let mut twos = vec![0.0; seq.n];
        for i in 0..seq.n {
            twos[i] = (2.0f64).powi(i as i32);
        }
        let func = FuncSeq::new(seq.n, vec![twos]);
        let exp = func.approx_k(STEP, seq, vec![S0], vec![R0]);
        let abs = (e - exp[0]).abs();
        assert_eq!(abs <= DELTA, true);
    }
}

#[test]
fn two_func_test() {
    let mut rng = thread_rng();
    for _case in 0..TESTCASES {
        let mut v = vec![Complex::zero(); ARR_SIZE];
        let c = rng.gen_range(0.0..=1000.0);
        let e = rng.gen_range(0.0..=1.0);
        let e2 = rng.gen_range(0.0..=5.0);
        for i in 0..ARR_SIZE {
            v[i] = Complex::from(c * (2.0f64).powi(i as i32).powf(e) * (1.0 + (i as f64)).log10().powf(e2));
        }
        let seq = Sequence::from(v);
        let mut twos = vec![0.0; seq.n];
        for i in 0..seq.n {
            twos[i] = (2.0f64).powi(i as i32);
        }
        let mut log =  vec![0.0; seq.n];
        for i in 0..seq.n {
            log[i] = (1.0 + (i as f64)).log10();
        }
        let func = FuncSeq::new(seq.n, vec![twos, log]);
        let x = func.approx_k(STEP, seq, vec![S0, S0], vec![1.0, 1.0]);
        assert_eq!((e-x[0]).abs() <= DELTA, true);
        assert_eq!((e2-x[1]).abs() <= DELTA, true);
    }
}
