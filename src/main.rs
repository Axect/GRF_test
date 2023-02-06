use rustfft::{FftPlanner, num_complex::Complex};
use peroxide::fuga::*;
use std::env::args;

fn main() {
    let n = args().nth(1).unwrap().parse::<usize>().unwrap();
    let sigma = args().nth(2).unwrap().parse::<f64>().unwrap();

    let x = linspace(0, 1, n);
    let y = grf_1(n, sigma);

    //let x_new = linspace(0, 1, 101);
    //let cs = cubic_hermite_spline(&x, &y, Quadratic);
    //let y_new = cs.eval_vec(&x_new);

    let mut df = DataFrame::new(vec![]);
    df.push("x", Series::new(x));
    df.push("y", Series::new(y));
    df.print();

    df.write_nc("test.nc").unwrap();
}

/// Gaussian Random Fields using circulant embedding method 1
///
/// * Title: An Effective Method for Simulating Gaussian Random Fields
/// * Author: Grace Chan
fn grf_1(n: usize, sigma: f64) -> Vec<f64> {
    let g = (2f64 * (n - 1) as f64).log2().ceil() as i32;
    let mut m = 2f64.powi(g) as usize;
    let qa = loop {
        let c = circulant_embedding(m, n, |x| k(x, sigma));
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(m);
        let mut c_fft = c.iter().map(|x| Complex::new(*x, 0f64)).collect::<Vec<_>>();
        fft.process(&mut c_fft);
        let c_fft = c_fft.iter().map(|x| x.re).collect::<Vec<_>>();
        let c_min = c_fft.iter().min_by(|&x, &y| x.partial_cmp(y).unwrap()).unwrap();
        c_min.print();

        if c_min >= &0f64 {
            break c_fft.fmap(|t| t.sqrt());
        } else if c_min.abs() < 1e-6 {
            break c_fft.fmap(|t| trunc(t).sqrt());
        } else {
            m *= 2;
        }
    };

    let normal = Normal(0f64, 1f64);
    let z = normal.sample(m);

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(m);
    let mut z_fft = z.iter().map(|x| Complex::new(*x, 0f64)).collect::<Vec<_>>();
    fft.process(&mut z_fft);

    let mut a = z_fft.into_iter()
        .zip(qa.into_iter())
        .map(|(x, y)| x * y)
        .collect::<Vec<_>>();

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(m);
    fft.process(&mut a);

    let y = a
        .iter()
        .map(|x| x.re / (m as f64))
        .collect::<Vec<_>>();

    y[..n].to_vec()
}

fn circulant_embedding<F: Fn(f64) -> f64>(m: usize, n: usize, kernel: F) -> Vec<f64> {
    let mut c = vec![0f64; m];
    let mid = m / 2;

    for i in 0 .. mid + 1 {
        c[i] = kernel(i as f64 / n as f64);
    }
    for i in mid + 1 .. m {
        c[i] = c[m - i];
    }

    c
}

fn trunc(x: f64) -> f64 {
    if x < 0f64 {
        0f64
    } else {
        x
    }
}

// Stationary Gaussian Kernel
fn k(dx: f64, sigma: f64) -> f64 {
    (-dx.powi(2) / (2.0 * sigma.powi(2))).exp()
}
