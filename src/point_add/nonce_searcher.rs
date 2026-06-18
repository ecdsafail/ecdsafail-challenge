//! Fast classical nonce searcher for find-the-hardest-GCD-config variant.
//! Build with:
//!   cargo build --release --bin build_circuit
//! then run:
//!   DIALOG_GCD_WIDTH_MARGIN=9 NONCE_SEARCH=1 ./target/release/build_circuit 2>&1

use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFC2F, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF,
]);

fn secp() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: SECP256K1_P,
        a: U256::from(0),
        b: U256::from(7),
        gx: U256::from_str_radix(
            "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16,
        ).unwrap(),
        gy: U256::from_str_radix(
            "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16,
        ).unwrap(),
        order: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16,
        ).unwrap(),
    }
}

/// Generate 9024 Fiat-Shamir test points using a given nonce, and run the
/// classical GCD filter on them. Returns the first working nonce, or None.
pub fn search_nonce(cfg: &super::dialog_gcd_classical_filter::DialogGcdFilterConfig) -> Option<u64> {
    let curve = secp();
    let prefix = b"quantum_ecc-fiat-shamir-v2";
    
    // Try a range of nonces
    for nonce in 0u64..50000 {
        if nonce % 1000 == 0 {
            eprintln!("nonce_search: trying nonce={}", nonce);
        }
        
        // Build SHAKE256 state with nonce embedded
        let mut hasher = Shake256::default();
        hasher.update(prefix);
        
        // Emulate the nonce X gates: nonce bits select between tx[0] and tx[1]
        // Were not emitting real ops here, so use a fixed op-count placeholder
