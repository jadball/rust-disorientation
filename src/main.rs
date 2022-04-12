extern crate nalgebra as na;

use std::f64::consts::PI;
use std::time::Instant;

use itertools::Itertools;
use na::{matrix, Matrix3, Rotation3};
use plotly::{Histogram, ImageFormat, Plot};
use plotly::histogram::HistNorm::ProbabilityDensity;
use rand::{distributions::Uniform, Rng};
use rayon::prelude::*;

fn main() {
    let now = Instant::now();
    let the_range = Uniform::from(0_f64..PI);
    let n_trials = 5000000;

    // Make a vector of 6-tuples of random numbers
    let random_vec: Vec<(f64, f64, f64, f64, f64, f64)> = (1..n_trials).map(|_|rand::thread_rng().sample_iter(&the_range).take(6).next_tuple().unwrap()).collect();
    // println!("{:?}", random_vec);
    let angles_vec = random_vec.into_par_iter().map(|x|angle_between_random(x)).collect::<Vec<f64>>();
    // println!("{:?}", angles_vec);
    let elapsed = now.elapsed();
    println!("Completed {} trials in {:.2?}", n_trials, elapsed);

    let trace = Histogram::new(angles_vec).name("h").hist_norm(ProbabilityDensity);
    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.show();
    plot.save("output.png", ImageFormat::PNG, 1280, 900, 1.0)

}

fn angle_between_random(x: (f64, f64, f64, f64, f64, f64)) -> f64{
    let (a1, a2, a3, a4, a5, a6) = x;
    let m1: Matrix3<f64> = Rotation3::from_euler_angles(a1, a2, a3).into_inner();
    let m2: Matrix3<f64> = Rotation3::from_euler_angles(a4, a5, a6).into_inner();
    let symmetries: Vec<Matrix3<f64>> = vec![
        matrix![1., 0., 0.; 0., 1., 0. ; 0., 0., 1.],
        matrix![0., 0., -1.; 0., -1., 0.; -1., 0., 0.],
        matrix![0., 0., -1.; 0., 1., 0. ; 1., 0., 0.],
        matrix![-1., 0., 0.; 0., 1., 0. ; 0., 0., -1.],
        matrix![0., 0., 1.; 0., 1., 0. ; -1., 0., 0.],
        matrix![1., 0., 0.; 0., 0., -1.; 0., 1., 0.],
        matrix![1., 0., 0.; 0., -1., 0.; 0., 0., -1.],
        matrix![1., 0., 0.; 0., 0., 1. ; 0., -1., 0.],
        matrix![0., -1., 0.; 1., 0., 0. ; 0., 0., 1.],
        matrix![-1., 0., 0.; 0., -1., 0.; 0., 0., 1.],
        matrix![0., 1., 0.; -1., 0., 0.; 0., 0., 1.],
        matrix![0., 0., 1.; 1., 0., 0. ; 0., 1., 0.],
        matrix![0., 1., 0.; 0., 0., 1. ; 1., 0., 0.],
        matrix![0., 0., -1.; -1., 0., 0.; 0., 1., 0.],
        matrix![0., -1., 0.; 0., 0., 1. ; -1., 0., 0.],
        matrix![0., 1., 0.; 0., 0., -1.; -1., 0., 0.],
        matrix![0., 0., -1.; 1., 0., 0. ; 0., -1., 0.],
        matrix![0., 0., 1.; -1., 0., 0.; 0., -1., 0.],
        matrix![0., -1., 0.; 0., 0., -1.; 1., 0., 0.],
        matrix![0., 1., 0.; 1., 0., 0. ; 0., 0., -1.],
        matrix![-1., 0., 0.; 0., 0., 1. ; 0., 1., 0.],
        matrix![0., 0., 1.; 0., -1., 0.; 1., 0., 0.],
        matrix![0., -1., 0.; -1., 0., 0.; 0., 0., -1.],
        matrix![-1., 0., 0.; 0., 0., -1.; 0., -1., 0.]
    ];
    angle_between(&m1, &m2, &symmetries)
}

fn angle_between(m1: &Matrix3<f64>, m2:  &Matrix3<f64>, symmetries: &Vec<Matrix3<f64>>) -> f64{
    // Copyright (c) 2013-2019 Henry Proudhon
    // From pymicro.crystal.microstructure.Orientation.disorientation
    // Translated from Python into Rust
    // Pymicro can be found at https://github.com/heprom/pymicro
    let mut the_angle: f64 = PI;

    let forward_tuple: (&Matrix3<f64>, &Matrix3<f64>) = (m1, m2);
    let reverse_tuple: (&Matrix3<f64>, &Matrix3<f64>) = (m2, m1);
    let tuples_list = [forward_tuple, reverse_tuple];
    for a_tuple in tuples_list {
        let (g_a, g_b) = a_tuple;
        for symm in symmetries {
            let oj: Matrix3<f64> =  symm * g_a;
            for symm_2 in symmetries {
                let oi: Matrix3<f64> = symm_2 * g_b;
                let delta:  Matrix3<f64> = oi * oj.transpose();
                let mis_angle: f64 = misorientation_angle_from_delta(&delta);
                if mis_angle < the_angle {
                    the_angle = mis_angle;
                }
            }
        }
    }
    return the_angle.to_degrees()
}

fn misorientation_angle_from_delta(delta: &Matrix3<f64>) -> f64 {
    // Copyright (c) 2013-2019 Henry Proudhon
    // From pymicro.crystal.microstructure.Orientation.disorientation
    // Translated from Python into Rust
    // Pymicro can be found at https://github.com/heprom/pymicro
    let cw: f64 = 0.5 * (delta[(0,0)] + delta[(1,1)] + delta[(2,2)] - 1.);
    let cw = if cw > 1.0 { 1.0 } else { cw };
    let omega = cw.acos();
    return omega
}