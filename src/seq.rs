use num_complex::Complex;
use num_traits::Zero;
use rustfft::FftPlanner;

#[derive(Clone, Debug)]
pub struct Sequence {
    pub n: usize,
    pub data: Vec<Complex<f64>>,
}

impl Sequence {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![Complex::zero(); n],
        }
    }

    pub fn from(data: Vec<Complex<f64>>) -> Self {
        Self {
            n: data.len(),
            data,
        }
    }

    pub fn fft(&self) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(self.n);
        let mut s = self.data.clone();
        s.resize(2 * self.n, Complex::zero());
        fft.process(&mut s);
        Self {
            n: 2 * self.n,
            data: s,
        }
    }

    pub fn ifft(&self) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_inverse(self.n);
        let mut s = self.data.clone();
        fft.process(&mut s);
        Self { n: self.n, data: s }
    }

    pub fn divide(&self, other: Self) -> Self {
        let mut seq = self.data.clone();
        if self.n < other.n {
            Self {
                n: self.n,
                data: self.data.clone(),
            }
        } else {
            let mut qut = vec![];
            for i in 0..=(self.n - other.n) {
                let q = seq[self.n - i - 1] / other.data[other.n - 1];
                for j in 1..=other.n {
                    seq[self.n - i - j] -= q * other.data[other.n - j];
                }
                qut.push(q);
            }
            qut.reverse();
            Self {
                n: self.n - other.n + 1,
                data: qut,
            }
        }
    }

    pub fn rem(&self, other: Self) -> Self {
        let mut seq = self.data.clone();
        if self.n < other.n {
            Self {
                n: self.n,
                data: self.data.clone(),
            }
        } else {
            for i in 0..=(self.n - other.n) {
                let q = seq[self.n - i - 1] / other.data[other.n - 1];
                for j in 1..=other.n {
                    seq[self.n - i - j] -= q * other.data[other.n - j];
                }
            }
            Self {
                n: other.n - 1,
                data: seq.get(0..(other.n - 1)).unwrap().to_vec(),
            }
        }
    }

    pub fn conv(&self, other: Self) -> Self {
        assert_eq!(self.n, other.n);
        let mut seq = vec![Complex::zero();self.n];
        for i in 0..self.n {
            for k in 0..self.n {
                seq[i] += self.data[k]*other.data[(self.n+i-k)%self.n];
            }
        }
        Self {
            n: self.n,
            data: seq,
        }
    }

    pub fn size(&self) -> f64 {
        let mut res = 0.0;
        for i in 0..self.n {
            res += self.data[i].re.powi(2) + self.data[i].im.powi(2);
        }
        res
    }
}
