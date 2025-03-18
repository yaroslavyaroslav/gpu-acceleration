use ark_bn254::{Fr as ScalarField, G1Projective as G};
use ark_ec::Group;
use ark_ff::{Field, PrimeField, UniformRand};
use ark_std::Zero;

pub fn horners_method(points: &[G], chunk_size: usize) -> G {
    let m = ScalarField::from((1 << chunk_size) as u64);
    let mut result = points.last().unwrap().clone();
    for i in (0..points.len() - 1).rev() {
        result = result * &m + points[i];
    }
    result
}

/// A naive method to calculate:
/// Q = ∑ 2^{(j-1)s} G_j (j=1...⌈λ/s⌉)
/// where λ is 254 (scalar bit size for BN254) and s is the chunk_size
pub fn naive_method(points: &[G], chunk_size: usize) -> G {
    let lambda = ScalarField::MODULUS_BIT_SIZE as usize;
    let num_chunks = (lambda + chunk_size - 1) / chunk_size; // Ceiling of λ/s
    
    // Ensure we have enough points
    assert!(points.len() >= num_chunks, "Not enough points provided for the given chunk size");
    
    let mut result = G::zero();

    for j in 1..=num_chunks {
        // Calculate 2^{(j-1)s}
        let exponent = (j - 1) * chunk_size;
        let scalar = if exponent < 64 {
            // For small exponents, directly use u64
            ScalarField::from(1u64 << exponent)
        } else {
            // For larger exponents, use 2^exponent directly
            let power_of_two = ScalarField::from(2u64).pow([exponent as u64]);
            power_of_two
        };
        
        // Add 2^{(j-1)s} * Gj to the result
        result = result + points[j-1] * &scalar;
    }
    
    result
}

#[test]
#[serial_test::serial]
fn test_horners_method() {
    // Calculate how many points we need for the test
    let chunk_size = 16;
    let lambda = 254; // Scalar bit size for BN254
    let num_chunks = (lambda + chunk_size - 1) / chunk_size; // Ceiling of λ/s
    
    let mut rng = rand::thread_rng();
    let points = vec![G::generator() * &ScalarField::rand(&mut rng); num_chunks];
    
    let expected = naive_method(&points, chunk_size);
    let result = horners_method(&points, chunk_size);
    assert!(expected == result);
}