//! Classical numeric replays of single-inversion point-add strategies.
//!
//! Discipline: every strategy here is a *reversible-schedule simulation*.
//! Each `replay_strategy_X()` mirrors, step by step, what a quantum
//! scaffold would do — tracked as a U256 per "register". The end state
//! MUST match the reference `WeierstrassEllipticCurve::add`, across 200
//! random curve-point pairs. A strategy that doesn't pass 200/200 is
//! dead, full stop.
//!
//! No strategy graduates to reversible code in `mod.rs` unless it
//! passes and its op count beats the current 2-Kaliski scaffold.
//!
//! See `single_inv_plan.md` for the prose spec / register tables.

#![cfg(test)]

use alloy_primitives::{U256, U512};

use super::SECP256K1_P;

fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

fn neg_mod(a: U256, p: U256) -> U256 {
    sub_mod(U256::ZERO, a, p)
}

// ─────────────────────────────────────────────────────────────────────
// Reference: the established, classically-correct formula (sanity).
// ─────────────────────────────────────────────────────────────────────
pub fn single_inv_add(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let dx = sub_mod(px, qx, p);
    let dy = sub_mod(py, qy, p);
    let inv_dx = dx.inv_mod(p).expect("dx nonzero");
    let lam = dy.mul_mod(inv_dx, p);
    let lam2 = lam.mul_mod(lam, p);
    let rx = sub_mod(sub_mod(lam2, px, p), qx, p);
    let ry = sub_mod(lam.mul_mod(sub_mod(qx, rx, p), p), qy, p);
    (rx, ry)
}

// ─────────────────────────────────────────────────────────────────────
// Kaliski raw output convention (settled in 21b87fd):
//   inside a `with_kal_inv_raw(v)` body the inv_raw register holds
//     -v^{-1} * 2^{2n-1} mod p
//   (sign negative, scale 2^{2n-1} because iters = 2n-1 and Kaliski
//    doubles r every iter; the final positive-ation is skipped).
// ─────────────────────────────────────────────────────────────────────
fn kaliski_body_inv_raw(v: U256, p: U256) -> U256 {
    const TWO_2N_MINUS_1_EXP: u64 = 2 * 256 - 1;
    let two = U256::from(2);
    let scale = two.pow_mod(U256::from(TWO_2N_MINUS_1_EXP), p);
    let inv_v = v.inv_mod(p).expect("v nonzero");
    neg_mod(inv_v, p).mul_mod(scale, p)
}

fn pow2_mod(e: i64, p: U256) -> U256 {
    let two = U256::from(2);
    if e >= 0 {
        two.pow_mod(U256::from(e as u64), p)
    } else {
        two.pow_mod(U256::from((-e) as u64), p)
            .inv_mod(p)
            .expect("2 invertible")
    }
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY A
//
// Plan (from single_inv_plan.md §4):
//   1. tx = dx, ty = dy
//   2. a := tx * ty                            (fresh register)
//   3. run Kaliski on a, body entered with inv_raw_a = -a^{-1}*2^{2n-1}
//   4. inside body, lam := ty^2 * inv_raw_a * 2^{-(2n-1)} = -λ
//      (implemented as: dy_sq then mul by inv_raw then halve 2n-1 times)
//   5. Rx fold into tx: tx -= λ²; +3Qx; neg; +Qx  → tx = Rx - Qx.
//   6. ty += lam · tx; then ty -= Qy.          ← This is where the
//      ty=dy contamination has to clear. Expected: Ry mismatch.
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_a(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let mut tx = sub_mod(px, qx, p);
    let mut ty = sub_mod(py, qy, p);

    // step 2: a = tx·ty
    let a = tx.mul_mod(ty, p);

    // step 3: enter kaliski body
    let inv_raw_a = kaliski_body_inv_raw(a, p);

    // step 4: lam = ty^2 · inv_raw_a · 2^{-(2n-1)} = -λ.
    let scale_back = pow2_mod(-(2 * 256 - 1), p);
    let lam = ty
        .mul_mod(ty, p)
        .mul_mod(inv_raw_a, p)
        .mul_mod(scale_back, p);
    // sanity: lam == -dy/dx
    let lam_expected = neg_mod(ty.mul_mod(tx.inv_mod(p).unwrap(), p), p);
    assert_eq!(lam, lam_expected, "strategy A: lam must equal -λ");

    // step 5: Rx fold.  Goal: tx = Rx when done.
    //   tx := dx - λ² ; tx += 2Qx ; tx := -tx
    //      → tx = -(dx - λ² + 2Qx) = λ² - dx - 2Qx
    //          = λ² - (Px - Qx) - 2Qx = λ² - Px - Qx = Rx. ✓
    let lam2 = lam.mul_mod(lam, p);
    tx = sub_mod(tx, lam2, p);
    tx = tx.add_mod(qx.mul_mod(U256::from(2), p), p);
    tx = neg_mod(tx, p);

    // At this point tx = Rx. For Ry we need (Rx - Qx):
    let rx_minus_qx = sub_mod(tx, qx, p);

    // step 6: ty += lam·(Rx - Qx) - Qy.  lam = -λ, so
    //   lam · (Rx - Qx) = (-λ)(Rx - Qx) = λ(Qx - Rx) = Ry + Qy.
    //   ty = dy + (Ry + Qy) = (Py - Qy) + Ry + Qy = Py + Ry.
    //   ty -= Qy → ty = Py + Ry - Qy.
    // That is the predicted contamination.
    ty = ty.add_mod(lam.mul_mod(rx_minus_qx, p), p);
    ty = sub_mod(ty, qy, p);

    (tx, ty)
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY B2
//
// Plan (from §5, the version that only runs ONE Kaliski pass AND does
//       Rx/Ry computation inside the Kaliski body so lam can be
//       uncomputed via mul2-inverse before exit):
//
//   Body over tx = dx, with inv_raw = -dx^{-1} · 2^{2n-1}.
//   1. lam := ty · inv_raw                    ← lam = -dy·dx^{-1}·2^{2n-1}
//   2. halve lam (2n-1) times                 ← lam = -λ
//   3. ty += lam · tx                         ← ty := dy - λ·dx = 0 ✅
//      (at this point ty is zero and Py is gone from the state)
//   4. tx ← Rx - Qx via fold (uses lam² only)
//   5. ty += lam · tx                         ← ty := 0 + (-λ)(Rx-Qx)
//                                                  = λ(Qx - Rx) = Ry + Qy
//   6. ty -= Qy                               ← ty := Ry ✅
//   7. Now reverse-mul and reverse-halve and reverse-mul to put lam
//      back at 0 using only current state.
//
//   Step 7 details: we're still inside the Kaliski body (inv_raw, tx=dx,
//   ty=Ry, lam=-λ all live; kal_state also live). To free lam (=-λ):
//     - un-fold tx: reverse the Rx fold so tx returns to dx.
//     - subtract `ty_before_Ry_step` back from ty. We don't have it, but
//       we can recompute it: it was (Ry + Qy) prior to the -Qy. So
//       `ty += Qy` → ty = Ry + Qy. Then `ty -= lam · tx` → ty = 0. ✅
//     - But we need ty = dy at body exit for the Kaliski_backward pass
//       to close out correctly? Actually no — Kaliski_backward only
//       reverses the Kaliski state, it doesn't touch ty. ty can be
//       anything at body exit.
//     - So we can actually just LEAVE ty = Ry, skip the reverse steps
//       for ty, and only uncompute lam.
//     - To uncompute lam: reverse halves (2n-1 doublings) brings lam
//       back to `-dy·dx^{-1}·2^{2n-1} = -dy·inv_raw`. Then the inverse
//       of `lam += ty·inv_raw` is `lam -= ty·inv_raw`. But here ty has
//       changed from dy to Ry! So we need to do uncomputation BEFORE ty
//       is changed — i.e. before step 3.
//
// Revised Strategy B2:
//   1. lam := ty · inv_raw = -dy·dx^{-1}·2^{2n-1}
//   2. halve lam (2n-1) times  → lam = -λ
//   3. Take a Bennett snapshot of λ:
//        lam2 := 0
//        lam2 := lam · 1 (just cx-copy)
//   Actually cx-copy of lam into a fresh register is free in qubits
//   (just n qubits, which is what we'd pay anyway). Let's do that.
//        lam_out := lam   (fresh register, cx-copy, classical)
//   4. reverse: double lam (2n-1) times → lam = -dy·inv_raw
//   5. lam -= ty · inv_raw       → lam = 0, freeable ✅
//   6. But now we still have lam_out = -λ live. Proceed with it.
//   7. tx ← Rx - Qx fold using lam_out² (= λ²).
//   8. ty += lam_out · tx        → ty = dy + (-λ)(Rx-Qx)
//                                    = dy + λ(Qx-Rx)
//                                    = (Py-Qy) + (Ry+Qy) = Py + Ry. ✗
//
// SAME OBSTRUCTION. Ry replaces Qy but Py stays. The classical replay
// will catch this as an Ry mismatch.
//
// The ONLY way to avoid the Py contamination is to zero ty BEFORE
// computing Ry. That's what the live 2-Kaliski scaffold does via
// pair1_mul2 (ty := 0) then mul3_between_pair (ty := λ(Rx-Qx)) then
// pair2_cleanup (ty -= Qy). So inside a single Kaliski we MUST:
//   - zero ty inside the body (step 3 above), then
//   - compute Ry into the now-zero ty WITHOUT needing λ outside the body.
// That requires the Rx fold to happen inside the body too, and lam to
// STILL be live when we compute Ry.
//
// Final revised Strategy B2 (what we actually run):
//   body:
//     1. lam := ty · inv_raw                   (lam = -dy·dx^{-1}·2^{2n-1})
//     2. halve lam (2n-1)×                     (lam = -λ)
//     3. ty += lam · tx                        (ty = 0)
//     4. Rx fold in tx using lam²              (tx = Rx - Qx)
//     5. ty += lam · tx                        (ty = (-λ)(Rx-Qx) = Ry+Qy)
//     6. ty -= Qy                              (ty = Ry)
//     7. un-fold tx (reverse of step 4)         (tx = dx, lam² available)
//     8. double lam (2n-1)×                     (lam = -dy·dx^{-1}·2^{2n-1})
//     9. lam -= ty · inv_raw                    (...but ty=Ry not dy!)
//
//   Step 9 again: we need ty = dy to undo step 1, and ty = Ry now.
//   Ry - dy relation: Ry = λ(Qx-Rx) - Qy = λ(Qx - Rx) - Qy
//                     = -λ(Rx - Qx) - Qy.
//   No usable relation between Ry and dy without knowing λ and Rx.
//
//   The escape: uncompute lam via an ancilla BEFORE ty changes.
//
//   Fresh final Strategy B2:
//     body:
//      1. lam := ty · inv_raw                  (lam = -dy·inv_raw)
//      2. halve lam (2n-1)×                    (lam = -λ)
//      3. lam_copy(n) := 0 ; lam_copy ^= lam   (cx-copy; lam_copy = -λ)
//      4. double lam (2n-1)×                   (lam = -dy·inv_raw again)
//      5. lam -= ty · inv_raw                  (lam = 0, free it) ✅
//      6. ty += lam_copy · tx                  (ty = dy - λ·dx = 0) ✅
//      7. Rx fold in tx using lam_copy²        (tx = Rx - Qx)
//      8. ty += lam_copy · tx                  (ty = Ry + Qy)
//      9. ty -= Qy                             (ty = Ry)
//     10. un-fold tx                            (tx = dx briefly?)
//     11. Reverse lam_copy: …
//
//   At end we must free lam_copy. Its current value is -λ. To uncompute
//   without needing dy: compute `lam_copy -= lam_recomputed`, where
//   lam_recomputed = -dy·dx^{-1}·? ... no, dy is gone.
//
//   But we CAN recompute -λ using the current live registers:
//   tx = dx (after step 10), Qy, ty = Ry, Qx, Px (from tx+Qx).
//     Ry + Qy = λ(Qx - Rx) = -λ(Rx - Qx)
//     So λ = -(Ry + Qy) / (Rx - Qx).
//   That requires inverting (Rx - Qx). We already inverted dx. Can we
//   get (Rx - Qx)^{-1} from dx^{-1}?
//     Rx - Qx = λ² - dx - 2Qx, where λ = dy·dx^{-1}
//            = dy²/dx² - dx - 2Qx
//            = (dy² - dx³ - 2Qx·dx²) / dx²
//     So (Rx - Qx) = (dy² - dx·(dx² + 2Qx·dx)) / dx²
//              1 / (Rx - Qx) = dx² / (dy² - dx·(dx² + 2Qx·dx))
//   That's a second inversion (of dy² - dx·(dx² + 2Qx·dx)) unless we
//   have other structure. Dead for single-Kaliski.
//
// *** Therefore the classical replay below implements the REVISED
// Strategy B2 and expects Ry to fail, pinpointing the obstruction at
// step 11. ***
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_b2(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let mut tx = sub_mod(px, qx, p);
    let mut ty = sub_mod(py, qy, p);
    let dx_original = tx;

    // enter kaliski body on tx = dx.
    let inv_raw = kaliski_body_inv_raw(dx_original, p);

    // (1) lam := ty · inv_raw
    let mut lam = ty.mul_mod(inv_raw, p);
    // (2) halve lam 2n-1 times
    lam = lam.mul_mod(pow2_mod(-(2 * 256 - 1), p), p);
    // Now lam = -λ.

    // (3) cx-copy lam into lam_copy (reversibly, 0 Toffoli, just tracking)
    let lam_copy = lam;
    // (4) double lam 2n-1 times → lam = -dy·inv_raw again
    lam = lam.mul_mod(pow2_mod(2 * 256 - 1, p), p);
    // (5) lam -= ty · inv_raw  (ty is still dy here)
    lam = sub_mod(lam, ty.mul_mod(inv_raw, p), p);
    assert_eq!(lam, U256::ZERO, "lam should be zero after uncompute");

    // (6) ty += lam_copy · tx = dy + (-λ)(dx) = 0.
    ty = ty.add_mod(lam_copy.mul_mod(tx, p), p);
    assert_eq!(ty, U256::ZERO, "ty should be zero after step 6");

    // (7) Rx fold in tx using lam_copy².  Goal tx := Rx.
    let lam_sq = lam_copy.mul_mod(lam_copy, p);
    tx = sub_mod(tx, lam_sq, p);
    tx = tx.add_mod(qx.mul_mod(U256::from(2), p), p);
    tx = neg_mod(tx, p); // tx = Rx

    // For Ry we need Rx - Qx:
    let rx_minus_qx = sub_mod(tx, qx, p);

    // (8) ty += lam_copy · (Rx - Qx) = (-λ)(Rx - Qx) = λ(Qx - Rx) = Ry + Qy.
    //     But ty is 0 here (step 6), so ty := Ry + Qy.
    ty = ty.add_mod(lam_copy.mul_mod(rx_minus_qx, p), p);
    // (9) ty -= Qy → ty = Ry.
    ty = sub_mod(ty, qy, p);

    // ANCILLA LEAK: lam_copy still = -λ. Must be freed. In classical
    // replay we don't have to zero it — but as a reversible scaffold it
    // would need to. The open question is whether `lam_copy` can be
    // uncomputed from {tx=Rx, ty=Ry, ox=Qx, oy=Qy, dx_original} alone.
    // In the classical trace we just note the leak and return.
    let _ = dx_original;
    let _ = lam_copy;
    (tx, ty)
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY C — Montgomery batch: invert a single product w that yields
// both 1/dx and 1/(Rx-Qx).
//
// From §6: `dx² · (Rx - Qx) = dy² - dx² · (Px + Qx)`, so
//          `dx · (Rx - Qx) = (dy² - dx² · (Px + Qx)) / dx`.
// We can compute the numerator `u = dy² - dx² · (Px + Qx)` classically
// from live registers {dx, dy, Px_const? no Px quantum, Qx}. Px is
// quantum. But Px = dx + Qx. So:
//   u = dy² - dx²·((dx + Qx) + Qx) = dy² - dx³ - 2·Qx·dx²
// That's a polynomial in dx, dy, Qx. Fully computable from the live
// quantum registers.
//
// Letting `w = u / dx = (dy² - dx³ - 2·Qx·dx²) / dx = dy²/dx - dx² - 2·Qx·dx`
// we see w = dx·(Rx - Qx). We need w^{-1}. Then:
//   1/dx        = w^{-1} · (Rx - Qx) = w^{-1} · u/dx ...  hmm still dx
//   1/(Rx - Qx) = w^{-1} · dx
//
// We don't know Rx - Qx yet (need λ), so we can't compute w directly
// as `dx · (Rx - Qx)`. We compute it as `u / dx = (dy² - dx³ - 2Qx·dx²)/dx
//                                      = dy²·dx^{-1} - dx² - 2Qx·dx`.
// Dividing by dx needs 1/dx. Circular.
//
// Alternative: just take `w' = dx · u = dx·dy² - dx⁴ - 2Qx·dx³`. Then
// w'^{-1} = 1/(dx·u). Product of dx and u, inverted once. Then:
//   1/dx = u · w'^{-1}
//   1/u  = dx · w'^{-1}
// We get 1/dx from one inversion of a degree-4 polynomial in dx plus
// degree-2 in dy. That's two quantum muls (u, w' = dx·u) then inversion
// then two more muls to recover 1/dx. Five muls total before λ — more
// expensive than just inverting dx directly.
//
// But we ALSO get 1/u for free. And 1/u tells us λ² via:
//   u = dy² - dx²(Px+Qx)  =>  dy² = u + dx²(Px+Qx)
//   λ² = dy²/dx² = (u + dx²(Px+Qx))/dx² = u/dx² + Px + Qx
//   λ² - Px - Qx = u/dx².
//   So Rx = u/dx² and Rx - Qx = u/dx² - Qx. Hmm that's u·dx^{-2}.
//
// Simpler cleaner algebra:
//   Rx = λ² - Px - Qx
//   dx² · Rx = dy² - dx²·Px - dx²·Qx = dy² - dx²·(Px + Qx)
//   So Rx = (dy² - dx²·(Px+Qx)) / dx².
//   Let v = dy² - dx²·(Px+Qx). Then Rx = v/dx².
//   And for Ry:
//     Ry = λ(Qx - Rx) - Qy.
//     dx·Ry = dy·(Qx - Rx) - dx·Qy
//     dx³·Ry = dx²·dy·(Qx - Rx) - dx³·Qy
//            = dy·(dx²·Qx - dx²·Rx) - dx³·Qy
//            = dy·(dx²·Qx - v) - dx³·Qy
//   So Ry = (dy·(dx²·Qx - v) - dx³·Qy) / dx³.
//   Let w = dx³. Then:
//     Rx = v · w^{-1} · dx                 (since v/dx² = v·dx/dx³)
//     Ry = (dy·(dx²·Qx - v) - w·Qy) · w^{-1}
//
//   One inversion (of w = dx³), both outputs computable.
//
// This is the cleanest Strategy C candidate: invert dx³, recover
// (Rx, Ry) via muls on known-live quantum values. Counts later.
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_c(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let dx = sub_mod(px, qx, p);
    let dy = sub_mod(py, qy, p);

    // v = dy² - dx²·(Px + Qx)
    let dx2 = dx.mul_mod(dx, p);
    let dx3 = dx2.mul_mod(dx, p);
    let dy2 = dy.mul_mod(dy, p);
    let px_plus_qx = px.add_mod(qx, p);
    let v = sub_mod(dy2, dx2.mul_mod(px_plus_qx, p), p);

    // one inversion
    let w = dx3;
    let w_inv = w.inv_mod(p).expect("dx != 0");

    // Rx = v · w^{-1} · dx = v · (dx · w^{-1})
    let dx_winv = dx.mul_mod(w_inv, p);
    let rx = v.mul_mod(dx_winv, p);

    // Ry = (dy·(dx²·Qx - v) - w·Qy) · w^{-1}
    let dx2_qx = dx2.mul_mod(qx, p);
    let core = sub_mod(dx2_qx, v, p);
    let numer = sub_mod(dy.mul_mod(core, p), w.mul_mod(qy, p), p);
    let ry = numer.mul_mod(w_inv, p);

    (rx, ry)
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY E: slope-coordinate point-add permutation.
//
// This is a ground-up non-BY/non-denominator-history approach. It uses the
// slope m=dy/dx as the temporary coordinate and updates the x register by the
// involution
//
//     dx -> Rx = m² - dx - 2Qx.
//
// Then the y register is converted from slope to affine output by
//
//     m -> Ry = -m*(Rx-Qx) - Qy.
//
// Algebraically this is the cleanest one-division point-add map found so far.
// It is SOTA-shaped only if the required in-place variable multiply/divide can
// be made roughly schoolbook-cost and product-clean without inverse history.
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_e_slope_coordinate(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let dx = sub_mod(px, qx, p);
    let dy = sub_mod(py, qy, p);
    let m = dy.mul_mod(dx.inv_mod(p).expect("dx nonzero"), p);
    let rx = sub_mod(sub_mod(m.mul_mod(m, p), dx, p), qx.mul_mod(U256::from(2), p), p);
    let rx_minus_qx = sub_mod(rx, qx, p);
    let ry = sub_mod(neg_mod(m.mul_mod(rx_minus_qx, p), p), qy, p);
    (rx, ry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

    fn curve() -> WeierstrassEllipticCurve {
        WeierstrassEllipticCurve {
            modulus: SECP256K1_P,
            a: U256::from(0),
            b: U256::from(7),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            )
            .unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            )
            .unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            )
            .unwrap(),
        }
    }

    fn rand_u256(rng: &mut u64) -> U256 {
        let mut limbs = [0u64; 4];
        for l in &mut limbs {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *l = *rng;
        }
        U256::from_limbs(limbs) % SECP256K1_P
    }

    fn each_trial<F>(mut body: F)
    where
        F: FnMut(U256, U256, U256, U256, U256, U256),
    {
        let c = curve();
        let mut rng = 0xdead_beef_cafe_f00du64;
        let mut n = 0usize;
        let mut tried = 0usize;
        while n < 200 && tried < 2000 {
            tried += 1;
            let k1 = rand_u256(&mut rng);
            let k2 = rand_u256(&mut rng);
            let (px, py) = c.mul(c.gx, c.gy, k1);
            let (qx, qy) = c.mul(c.gx, c.gy, k2);
            if (px.is_zero() && py.is_zero())
                || (qx.is_zero() && qy.is_zero())
                || px == qx
            {
                continue;
            }
            let (rx_ref, ry_ref) = c.add(px, py, qx, qy);
            body(px, py, qx, qy, rx_ref, ry_ref);
            n += 1;
        }
        assert_eq!(n, 200, "needed 200 random valid trials");
    }

    #[test]
    fn reference_formula_sanity() {
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = single_inv_add(px, py, qx, qy);
            assert_eq!(rx, rx_ref);
            assert_eq!(ry, ry_ref);
        });
    }

    fn secp256k1_beta_endomorphism() -> U256 {
        U256::from_str_radix(
            "7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE",
            16,
        )
        .unwrap()
    }

    fn j0_endo_slope_numerator_xonly(x: U256, qx: U256, p: U256) -> U256 {
        x.mul_mod(x, p)
            .add_mod(qx.mul_mod(x, p), p)
            .add_mod(qx.mul_mod(qx, p), p)
    }

    #[test]
    fn secp_j0_endomorphism_slope_denominator_swap_identity_passes() {
        // secp256k1 has the j=0 automorphism (x,y)->(βx,y).  Therefore
        // (x-qx)(βx-qx)(β²x-qx)=x³-qx³=(y-qy)(y+qy), so the affine slope can
        // be written as
        //     λ = (y-qy)/(x-qx) = (x² + qx*x + qx²)/(y+qy).
        // This is a real special-curve identity: it swaps the denominator from
        // x-qx to y+qy and makes the numerator x-only quadratic.  It is kept as
        // a candidate algebraic tool, but the follow-up phase test below checks
        // that it does not by itself solve λ cleanup.
        let p = SECP256K1_P;
        let beta = secp256k1_beta_endomorphism();
        let beta2 = beta.mul_mod(beta, p);
        assert_eq!(beta.mul_mod(beta2, p), U256::from(1));
        assert_eq!(beta.add_mod(beta2, p).add_mod(U256::from(1), p), U256::ZERO);

        let mut checked = 0usize;
        let mut exceptional_y_sum = 0usize;
        each_trial(|px, py, qx, qy, _rx_ref, _ry_ref| {
            let den_y_sum = py.add_mod(qy, p);
            if den_y_sum.is_zero() {
                exceptional_y_sum += 1;
                return;
            }
            let dx = sub_mod(px, qx, p);
            let dy = sub_mod(py, qy, p);
            let lam_standard = dy.mul_mod(dx.inv_mod(p).expect("dx nonzero"), p);
            let product_form = sub_mod(beta.mul_mod(px, p), qx, p)
                .mul_mod(sub_mod(beta2.mul_mod(px, p), qx, p), p);
            let xonly_form = j0_endo_slope_numerator_xonly(px, qx, p);
            assert_eq!(product_form, xonly_form);
            let lam_endo = xonly_form.mul_mod(den_y_sum.inv_mod(p).expect("y+qy nonzero"), p);
            assert_eq!(lam_endo, lam_standard);
            checked += 1;
        });
        eprintln!(
            "secp j=0 endomorphism slope swap: checked={checked}, exceptional_y_sum={exceptional_y_sum}"
        );
        println!("METRIC endomorphism_slope_identity_samples={checked}");
        println!("METRIC endomorphism_slope_y_sum_exceptions={exceptional_y_sum}");
        assert!(checked >= 195, "random secp samples hit too many y+Qy exceptional cases");
    }

    #[test]
    fn strategy_a_rx_ok_ry_contaminated_by_dy() {
        // Strategy A is predicted DEAD. Specifically:
        //   Rx: should match reference (200/200).
        //   Ry: should equal `ref_Ry + dy` on every trial (the exact
        //       contamination predicted in the plan doc:
        //       ty_final = Py + Ry - Qy = Ry + dy).
        let mut rx_ok = 0;
        let mut ry_off_by_dy = 0;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_a(px, py, qx, qy);
            if rx == rx_ref {
                rx_ok += 1;
            }
            let dy = sub_mod(py, qy, SECP256K1_P);
            let ry_expected_with_bug = ry_ref.add_mod(dy, SECP256K1_P);
            if ry == ry_expected_with_bug {
                ry_off_by_dy += 1;
            }
        });
        eprintln!("Strategy A: rx ok {rx_ok}/200 ; ry == ref+dy in {ry_off_by_dy}/200");
        assert_eq!(rx_ok, 200, "Strategy A Rx must match all 200");
        assert_eq!(
            ry_off_by_dy, 200,
            "Strategy A Ry should be off by exactly +dy on every trial"
        );
    }

    #[test]
    fn strategy_b2_passes_but_leaks_lam_copy() {
        // Strategy B2 with the lam_copy ancilla trick: Rx and Ry both
        // pass, BUT lam_copy (= -λ) is a leaked ancilla we don't yet
        // know how to reversibly zero without a second inversion.
        // Cost accounting below ignores that leak.
        let mut rx_ok = 0;
        let mut ry_ok = 0;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_b2(px, py, qx, qy);
            if rx == rx_ref {
                rx_ok += 1;
            }
            if ry == ry_ref {
                ry_ok += 1;
            }
        });
        eprintln!(
            "Strategy B2: rx matches {rx_ok}/200, ry matches {ry_ok}/200 (lam_copy leaked)"
        );
        assert_eq!(rx_ok, 200);
        assert_eq!(ry_ok, 200);
    }

    /// Falsification test for B2: at Kaliski body exit we have the
    /// registers {tx=Rx (or Rx-Qx), ty=Ry, dx_orig=dx, inv_raw=-dx⁻¹·2^{2n-1}}
    /// plus the classical constants ox=Qx, oy=Qy.
    /// Can -λ be expressed as a small-cost polynomial combination of these?
    ///
    /// Enumerate a catalogue of candidate expressions (each one corresponds
    /// to a tiny mul/add sequence, ~1–3 muls of cost ~150k each). If ANY
    /// candidate equals -λ on 200 random trials, B2 is alive. If NONE do,
    /// B2 is dead.
    #[test]
    fn strategy_b2_lam_copy_uncompute_falsification() {
        let p = SECP256K1_P;
        // Candidate functions: take (rx, ry, qx, qy, dx, dx_inv) → guess for -λ.
        type Cand = fn(U256, U256, U256, U256, U256, U256, U256) -> U256;
        let neg = |a: U256| sub_mod(U256::ZERO, a, p);
        let _ = neg;

        // Polynomial in {rx, ry, qx, qy, dx, dx_inv, dx_inv^2}. dx_inv^2 is
        // "free" relative to adding a new inversion: it's one extra mul.
        let candidates: &[(&str, Cand)] = &[
            ("(Ry+Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("-(Ry+Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(U256::ZERO, ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P), SECP256K1_P)
            }),
            ("(Qy-Ry)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(qy, ry, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("(Ry-Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(ry, qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("(Rx-Qx)·dx_inv", |rx, _ry, qx, _qy, _dx, dx_inv, _| {
                sub_mod(rx, qx, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("dx_inv_sq", |_rx, _ry, _qx, _qy, _dx, _dx_inv, dx_inv_sq| dx_inv_sq),
            ("(Ry+Qy)·dx_inv_sq", |_rx, ry, _qx, qy, _dx, _dx_inv, dx_inv_sq| {
                ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("(Rx-Qx)·dx_inv_sq", |rx, _ry, qx, _qy, _dx, _dx_inv, dx_inv_sq| {
                sub_mod(rx, qx, SECP256K1_P).mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("dx·dx_inv_sq", |_rx, _ry, _qx, _qy, dx, _dx_inv, dx_inv_sq| {
                dx.mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("(Rx+Qx)·dx_inv", |rx, _ry, qx, _qy, _dx, dx_inv, _| {
                rx.add_mod(qx, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
        ];

        let mut hits = vec![0usize; candidates.len()];
        let mut total = 0usize;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            total += 1;
            let dx = sub_mod(px, qx, p);
            let dy = sub_mod(py, qy, p);
            let lam = dy.mul_mod(dx.inv_mod(p).unwrap(), p);
            let neg_lam = sub_mod(U256::ZERO, lam, p);
            let dx_inv = dx.inv_mod(p).unwrap();
            let dx_inv_sq = dx_inv.mul_mod(dx_inv, p);
            for (i, (_name, f)) in candidates.iter().enumerate() {
                let got = f(rx_ref, ry_ref, qx, qy, dx, dx_inv, dx_inv_sq);
                if got == neg_lam {
                    hits[i] += 1;
                }
            }
        });
        let mut any_matched_all = false;
        for (i, (name, _)) in candidates.iter().enumerate() {
            eprintln!("  candidate {name}: {}/{} matches -λ", hits[i], total);
            if hits[i] == total {
                any_matched_all = true;
            }
        }
        assert!(
            !any_matched_all,
            "B2 rescued: a low-cost polynomial expression in {{Rx, Ry, Qx, Qy, dx, dx⁻¹, dx⁻²}} \
             matches -λ on all trials; re-examine the obstruction"
        );
        eprintln!("B2 falsification: no low-cost polynomial in {{Rx, Ry, Qx, Qy, dx, dx⁻¹, dx⁻²}} equals -λ. B2 DEAD.");
    }

    #[test]
    fn strategy_c_passes_200() {
        // Strategy C: invert w = dx³, recover both Rx and Ry from it.
        // Everything is classical-reversible — only question is cost.
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_c(px, py, qx, qy);
            assert_eq!(rx, rx_ref);
            assert_eq!(ry, ry_ref);
        });
    }

    #[test]
    fn strategy_e_slope_coordinate_formula_passes_200() {
        // New ground-up attempt: convert to the line slope m=dy/dx, update
        // x by the involution dx -> Rx, then convert m to Ry. This validates
        // the algebra before any circuit work.
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_e_slope_coordinate(px, py, qx, qy);
            assert_eq!(rx, rx_ref);
            assert_eq!(ry, ry_ref);
        });
    }

    #[test]
    fn chord_product_identity_does_not_batch_the_two_affine_inversions() {
        // Tempting Montgomery-trick idea: the two affine denominators are
        //   dx = Px-Qx
        //   bx = Rx-Qx
        // If their product were available before the first inversion cleanup,
        // maybe one inversion of dx*bx could replace the second Kaliski.
        // The chord cubic gives an exact identity
        //   dx * (Rx-Qx) = 3 Qx^2 - 2 lambda Qy,
        // so the product is only linear in the already-computed slope.  But to
        // clean lambda we need 1/(Rx-Qx), and
        //   1/(Rx-Qx) = dx / (dx*(Rx-Qx)).
        // This has merely moved the second inversion from bx to the product d,
        // and then adds a variable multiply by dx.  It does not batch the two
        // inversions into one.
        let p = SECP256K1_P;
        let mut samples = 0usize;
        each_trial(|px, py, qx, qy, rx_ref, _ry_ref| {
            samples += 1;
            let dx = sub_mod(px, qx, p);
            let dy = sub_mod(py, qy, p);
            let lam = dy.mul_mod(dx.inv_mod(p).unwrap(), p);
            let bx = sub_mod(rx_ref, qx, p);
            let chord_product = dx.mul_mod(bx, p);
            let three_qx2 = qx.mul_mod(qx, p).mul_mod(U256::from(3u64), p);
            let two_lam_qy = lam.mul_mod(qy, p).mul_mod(U256::from(2u64), p);
            let derivative_identity = sub_mod(three_qx2, two_lam_qy, p);
            assert_eq!(
                chord_product, derivative_identity,
                "chord derivative identity for dx*(Rx-Qx) failed"
            );

            let inv_bx_direct = bx.inv_mod(p).unwrap();
            let inv_bx_via_product = dx.mul_mod(chord_product.inv_mod(p).unwrap(), p);
            assert_eq!(inv_bx_direct, inv_bx_via_product);
        });

        let current_second_inverse = 1_600_000usize;
        let added_variable_multiply_floor = 149_889usize;
        let product_trick_cost_floor = current_second_inverse + added_variable_multiply_floor;
        eprintln!(
            "Chord product identity checked on {samples} samples; replacing inv(Rx-Qx) by inv(dx*(Rx-Qx)) adds at least one variable multiply, floor={product_trick_cost_floor}"
        );
        println!("METRIC chord_product_identity_samples={samples}");
        println!("METRIC chord_product_second_inverse_moved_not_removed=1");
        println!("METRIC chord_product_alt_extra_mul_floor={added_variable_multiply_floor}");
        assert!(product_trick_cost_floor > current_second_inverse);
    }

    #[test]
    fn strategy_e_slope_coordinate_budget_requires_new_inplace_variable_multiply() {
        // The slope-coordinate map has one division plus one in-place variable
        // multiplication m -> -m*(Rx-Qx)-Qy. Known reversible ways to make
        // that multiplication product-clean are equivalent to the pair2
        // product-clean primitive already measured. This budget is the early
        // invalidation gate: current primitives miss SOTA, while a genuinely
        // new schoolbook-like in-place variable multiply would be worth wiring.
        let non_div_scaffold_after_one_div = 942_750.0;
        let compact_div_target = 900_000.0;
        let known_product_clean = 1_145_760.0;
        let schoolbook_like_product_target = 180_000.0;
        let current_known_total = non_div_scaffold_after_one_div + compact_div_target + known_product_clean;
        let target_if_new_mul = non_div_scaffold_after_one_div + compact_div_target + schoolbook_like_product_target;
        eprintln!(
            "Strategy E slope-coordinate budget: current_known≈{current_known_total:.0}, if_new_inplace_mul≈{target_if_new_mul:.0}, need_new_mul_saving≈{:.0}",
            known_product_clean - schoolbook_like_product_target
        );
        assert!(current_known_total > 2_700_000.0, "known product-clean primitive would already be SOTA-shaped; wire Strategy E");
        assert!(target_if_new_mul < 2_100_000.0, "even a schoolbook-cost in-place variable multiply would not make Strategy E worthwhile");
    }

    fn destructive_montgomery_step(t: u64, a: u64, bit: u64, p: u64) -> u64 {
        let mut u = t + bit * a;
        if (u & 1) != 0 {
            u += p;
        }
        u >> 1
    }

    fn destructive_montgomery_block(mut t: u64, a: u64, bits_lsb_first: u64, k: usize, p: u64) -> u64 {
        for i in 0..k {
            t = destructive_montgomery_step(t, a, (bits_lsb_first >> i) & 1, p);
        }
        t
    }

    #[test]
    fn destructive_montgomery_product_is_algebraically_promising_but_not_locally_reversible() {
        // Attempt after Strategy E: make the missing in-place multiply by
        // destructively scanning the multiplier bits through a Montgomery
        // add-and-halve accumulator.  Forward algebra is promising: for an
        // n-bit prime p, n steps output a*b*2^-n (mod p), up to a final p
        // subtraction.  If consumed multiplier bits were recoverable from the
        // accumulator, this would be a schoolbook-like product-clean primitive.
        let p = 251u64;
        let n = 8usize;
        let a = 173u64;
        let b = 123u64;
        let t = destructive_montgomery_block(0, a, b, n, p);
        let r_inv = U256::from(1u64 << n).inv_mod(U256::from(p)).unwrap();
        let expected = U256::from(a)
            .mul_mod(U256::from(b), U256::from(p))
            .mul_mod(r_inv, U256::from(p));
        assert_eq!(U256::from(t % p), expected);

        // Fast invalidation: after an 8-bit destructive window, the post-window
        // accumulator does NOT determine the consumed input bits and prior
        // accumulator.  For this concrete reachable poststate there are 512
        // valid (old_t, consumed_bits) predecessors.  Therefore a reversible
        // circuit must keep history/checkpoints or compute a nonlocal inverse;
        // the hoped-for local bit clearing is dead.
        let post = destructive_montgomery_block(0, a, 0b1011_0110, n, p);
        let mut preimages = 0usize;
        for old_t in 0..(2 * p) {
            for bits in 0..(1u64 << n) {
                if destructive_montgomery_block(old_t, a, bits, n, p) == post {
                    preimages += 1;
                }
            }
        }
        eprintln!("destructive Montgomery window post={post}, preimages={preimages}");
        assert_eq!(post, 223);
        assert_eq!(preimages, 512);
    }

    fn quotient_phase_truth_table_anf_stats(n: usize, p: u16, mask: u16) -> (usize, usize) {
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let mut inv = vec![0u16; p as usize];
        for x in 1..p {
            for y in 1..p {
                if ((x as u32) * (y as u32)) % (p as u32) == 1 {
                    inv[x as usize] = y;
                    break;
                }
            }
        }
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let x = (idx as u16) & limb_mask;
            let z = ((idx >> n) as u16) & limb_mask;
            let q = if x != 0 && x < p && z < p {
                ((z as u32 * inv[x as usize] as u32) % p as u32) as u16
            } else {
                0
            };
            anf[idx] = ((q & mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn montgomery_q_history_xy(n: usize, p: u16, x: u16, y: u16) -> (u16, u16) {
        let mut t = 0u32;
        let mut q_hist = 0u16;
        for i in 0..n {
            if ((y >> i) & 1) != 0 {
                t += x as u32;
            }
            let q = (t & 1) as u16;
            q_hist |= q << i;
            if q != 0 {
                t += p as u32;
            }
            t >>= 1;
        }
        ((t % p as u32) as u16, q_hist)
    }

    fn montgomery_q_history_phase_anf_stats(n: usize, p: u16, mask: u16) -> (usize, usize) {
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let x = (idx as u16) & limb_mask;
            let y = ((idx >> n) as u16) & limb_mask;
            let q_hist = if x != 0 && x < p && y < p {
                montgomery_q_history_xy(n, p, x, y).1
            } else {
                0
            };
            anf[idx] = ((q_hist & mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    #[test]
    fn montgomery_q_history_mbuc_phase_is_high_degree_but_structured() {
        // Primitive-level rescue attempt for Strategy E: compute an in-place
        // product with the bit-serial Montgomery loop, MBUC-measure the
        // quotient/history q_i bits, and phase-correct those q_i instead of
        // dividing by the product source.  Unlike the final quotient z/x phase,
        // q-history is not dense on the n=8 toy field, but it still reaches
        // maximal algebraic degree.  This keeps it as a narrow research lead,
        // not an implementation target: we need a closed-form sparse phase
        // circuit before touching point-add wiring.
        let (degree, density) = montgomery_q_history_phase_anf_stats(8, 251, 0b1010_0101);
        eprintln!("Montgomery q-history phase ANF: degree={degree}, density={density}/65536");
        assert!(degree >= 15);
        assert!(density > 2_000);
    }

    fn montgomery_q_history_phase_anf_stats_xz(n: usize, p: u16, mask: u16) -> (usize, usize) {
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let r_mod = ((1u32 << n) % p as u32) as u16;
        let mut inv = vec![0u16; p as usize];
        for x in 1..p {
            for y in 1..p {
                if ((x as u32) * (y as u32)) % (p as u32) == 1 {
                    inv[x as usize] = y;
                    break;
                }
            }
        }
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let x = (idx as u16) & limb_mask;
            let z = ((idx >> n) as u16) & limb_mask;
            let q_hist = if x != 0 && x < p && z < p {
                // Montgomery output z = x*y*R^-1, so y = z*R*x^-1.
                let y = (((z as u32) * (r_mod as u32) * (inv[x as usize] as u32)) % p as u32) as u16;
                let (check_z, q) = montgomery_q_history_xy(n, p, x, y);
                assert_eq!(check_z, z);
                q
            } else {
                0
            };
            anf[idx] = ((q_hist & mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn inv_mod_u16_for_phase_test(a: u16, p: u16) -> u16 {
        for x in 1..p {
            if ((a as u32) * (x as u32)) % (p as u32) == 1 {
                return x;
            }
        }
        0
    }

    fn sub_mod_u16_for_phase_test(a: u16, b: u16, p: u16) -> u16 {
        ((a as u32 + p as u32 - b as u32) % p as u32) as u16
    }

    fn add_mod_u16_for_phase_test(a: u16, b: u16, p: u16) -> u16 {
        ((a as u32 + b as u32) % p as u32) as u16
    }

    fn mul_mod_u16_for_phase_test(a: u16, b: u16, p: u16) -> u16 {
        (((a as u32) * (b as u32)) % (p as u32)) as u16
    }

    fn curve_rhs_u16_for_phase_test(x: u16, p: u16) -> u16 {
        add_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(x, x, p), x, p), 7, p)
    }

    fn is_curve_point_u16_for_phase_test(x: u16, y: u16, p: u16) -> bool {
        mul_mod_u16_for_phase_test(y, y, p) == curve_rhs_u16_for_phase_test(x, p)
    }

    fn endomorphism_output_lambda_phase_anf_stats(n: usize, p: u16, mask: u16) -> (usize, usize) {
        // For j=0 curves with β in the base field, the same slope cleanup can
        // be expressed from output R as
        //     λ = (Rx² + Qx*Rx + Qx²)/(Qy - Ry)
        // because the line through Q and -R has slope λ.  If the endomorphism
        // denominator swap were a cheap cleanup, this phase should be sparse.
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let (qx, qy) = first_curve_point_u16_for_phase_test(p);
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let rx = (idx as u16) & limb_mask;
            let ry = ((idx >> n) as u16) & limb_mask;
            let lam = if rx < p && ry < p && is_curve_point_u16_for_phase_test(rx, ry, p) {
                let den = sub_mod_u16_for_phase_test(qy, ry, p);
                if den == 0 {
                    0
                } else {
                    let numer = add_mod_u16_for_phase_test(
                        add_mod_u16_for_phase_test(
                            mul_mod_u16_for_phase_test(rx, rx, p),
                            mul_mod_u16_for_phase_test(qx, rx, p),
                            p,
                        ),
                        mul_mod_u16_for_phase_test(qx, qx, p),
                        p,
                    );
                    mul_mod_u16_for_phase_test(numer, inv_mod_u16_for_phase_test(den, p), p)
                }
            } else {
                0
            };
            anf[idx] = ((lam & mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn point_sub_const_u16_for_phase_test(rx: u16, ry: u16, qx: u16, qy: u16, p: u16) -> Option<(u16, u16)> {
        if !is_curve_point_u16_for_phase_test(rx, ry, p) {
            return None;
        }
        // P = R - Q = R + (Qx,-Qy).
        let nqy = if qy == 0 { 0 } else { p - qy };
        if rx == qx && ry == qy {
            return None;
        }
        let dx = sub_mod_u16_for_phase_test(rx, qx, p);
        if dx == 0 {
            return None;
        }
        let dy = sub_mod_u16_for_phase_test(ry, nqy, p);
        let lam = mul_mod_u16_for_phase_test(dy, inv_mod_u16_for_phase_test(dx, p), p);
        let lam2 = mul_mod_u16_for_phase_test(lam, lam, p);
        let px = sub_mod_u16_for_phase_test(sub_mod_u16_for_phase_test(lam2, rx, p), qx, p);
        let py = sub_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(lam, sub_mod_u16_for_phase_test(qx, px, p), p), nqy, p);
        Some((px, py))
    }

    fn point_add_const_u16_for_phase_test(px: u16, py: u16, qx: u16, qy: u16, p: u16) -> Option<(u16, u16)> {
        if !is_curve_point_u16_for_phase_test(px, py, p) {
            return None;
        }
        if px == qx && add_mod_u16_for_phase_test(py, qy, p) == 0 {
            return None;
        }
        let lam = if px == qx && py == qy {
            if py == 0 {
                return None;
            }
            let num = mul_mod_u16_for_phase_test(3, mul_mod_u16_for_phase_test(px, px, p), p);
            let den = mul_mod_u16_for_phase_test(2, py, p);
            mul_mod_u16_for_phase_test(num, inv_mod_u16_for_phase_test(den, p), p)
        } else {
            let dx = sub_mod_u16_for_phase_test(qx, px, p);
            if dx == 0 {
                return None;
            }
            let dy = sub_mod_u16_for_phase_test(qy, py, p);
            mul_mod_u16_for_phase_test(dy, inv_mod_u16_for_phase_test(dx, p), p)
        };
        let rx = sub_mod_u16_for_phase_test(sub_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(lam, lam, p), px, p), qx, p);
        let ry = sub_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(lam, sub_mod_u16_for_phase_test(px, rx, p), p), py, p);
        Some((rx, ry))
    }

    fn top_level_measured_input_phase_anf_stats(n: usize, p: u16, qx: u16, qy: u16, mask: u16) -> (usize, usize) {
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let rx = (idx as u16) & limb_mask;
            let ry = ((idx >> n) as u16) & limb_mask;
            let phase_word = if rx < p && ry < p {
                point_sub_const_u16_for_phase_test(rx, ry, qx, qy, p)
                    .map(|(px, py)| px ^ (py << n))
                    .unwrap_or(0)
            } else {
                0
            };
            anf[idx] = ((phase_word & mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    #[test]
    fn top_level_mbuc_of_old_point_requires_dense_point_subtraction_phase() {
        // Another single-inversion escape would compute R out-of-place, measure
        // the old input point P, and phase-correct from the surviving output R
        // using P = R - Q.  On a toy secp256k1-shaped curve this phase oracle is
        // already high-degree/dense.  So generic top-level MBUC does not avoid
        // the affine reversibility wall; it just asks for a phase version of
        // point subtraction.
        let p = 251u16;
        let (qx, qy) = first_curve_point_u16_for_phase_test(p);
        let (degree, density) = top_level_measured_input_phase_anf_stats(8, p, qx, qy, 0b1010_0101_0101_1010);
        eprintln!(
            "Top-level old-point MBUC phase ANF: q=({qx},{qy}), degree={degree}, density={density}/65536"
        );
        assert!(degree >= 14);
        assert!(density > 10_000);
    }

    #[test]
    fn endomorphism_slope_swap_cleanup_phase_is_still_dense() {
        // The j=0 denominator swap is algebraically real and removes dy from
        // the *forward* slope numerator, but cleaning a measured/live λ from the
        // output still asks for a quotient-like phase:
        //     λ = (Rx² + Qx*Rx + Qx²)/(Qy - Ry).
        // Full-domain toy ANFs are already high-degree/dense, so this special
        // secp256k1 automorphism is not a free one-inversion cleanup.
        let cases = [
            (8usize, 241u16, 0b1010_0101u16, 15usize, 20_000usize),
            (10usize, 1009u16, 0b10_1001_0101u16, 18usize, 250_000usize),
        ];
        for &(n, p, mask, min_degree, min_density) in &cases {
            let (degree, density) = endomorphism_output_lambda_phase_anf_stats(n, p, mask);
            let table = 1usize << (2 * n);
            eprintln!(
                "Endomorphism slope-cleanup phase ANF: n={n}, p={p}, degree={degree}/{}, density={density}/{table}",
                2 * n
            );
            if n == 10 {
                println!("METRIC endomorphism_slope_phase_degree_n10={degree}");
                println!("METRIC endomorphism_slope_phase_density_n10={density}");
            }
            assert!(degree >= min_degree);
            assert!(density >= min_density);
        }
    }

    #[test]
    fn endomorphism_slope_support_degree_still_grows() {
        // The full-domain ANF above is pessimistic.  Repeat the interpolation
        // only on valid curve outputs for j=0 toy primes (p≡1 mod 3).  The
        // minimum degree still grows with n, matching the earlier ordinary
        // lambda cleanup story rather than revealing a constant-degree
        // automorphism phase.
        let cases = [
            (4usize, 13u16, 0b1010u16, 3usize),
            (6usize, 61u16, 0b10_1010u16, 4usize),
            (8usize, 241u16, 0b1010_0101u16, 5usize),
            (10usize, 1009u16, 0b10_1001_0101u16, 6usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 6usize),
        ];
        let mut last_min = 0usize;
        for &(n, p, mask, expected_upper) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = endomorphism_curve_support_lambda_phase_min_degree(
                n,
                p,
                qx,
                qy,
                mask,
                expected_upper,
            )
            .expect("endomorphism support phase should interpolate by expected degree");
            eprintln!(
                "Endomorphism support-restricted slope phase: n={n}, p={p}, q=({qx},{qy}), min_degree={min_degree}"
            );
            if n == 12 {
                println!("METRIC endomorphism_slope_support_min_degree_n12={min_degree}");
            }
            assert!(min_degree >= last_min, "support degree unexpectedly decreased");
            last_min = min_degree;
        }
    }

    fn monomial_masks_for_curve_phase_test(vars: usize, max_degree: usize) -> Vec<u32> {
        fn rec(out: &mut Vec<u32>, vars: usize, start: usize, left: usize, acc: u32) {
            if left == 0 {
                out.push(acc);
                return;
            }
            for bit in start..=vars - left {
                rec(out, vars, bit + 1, left - 1, acc | (1u32 << bit));
            }
        }
        let mut masks = vec![0u32];
        for degree in 1..=max_degree {
            rec(&mut masks, vars, 0, degree, 0);
        }
        masks
    }

    fn gf2_rank_bitrows_for_curve_phase_test(rows: &mut [Vec<u64>], bits: usize) -> usize {
        let mut rank = 0usize;
        for col in 0..bits {
            let word = col / 64;
            let mask = 1u64 << (col % 64);
            let pivot = (rank..rows.len()).find(|&r| (rows[r][word] & mask) != 0);
            if let Some(p) = pivot {
                rows.swap(rank, p);
                let pivot_row = rows[rank].clone();
                for r in 0..rows.len() {
                    if r != rank && (rows[r][word] & mask) != 0 {
                        for w in word..rows[r].len() {
                            rows[r][w] ^= pivot_row[w];
                        }
                    }
                }
                rank += 1;
                if rank == rows.len() {
                    break;
                }
            }
        }
        rank
    }

    fn curve_support_old_point_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let vars = 2 * n;
        let masks = monomial_masks_for_curve_phase_test(vars, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut rows = Vec::new();
        for rx in 0..p {
            for ry in 0..p {
                if !is_curve_point_u16_for_phase_test(rx, ry, p) {
                    continue;
                }
                let Some((px, py)) = point_sub_const_u16_for_phase_test(rx, ry, qx, qy, p) else {
                    continue;
                };
                let idx = (rx as u32) | ((ry as u32) << n);
                let phase_word = px ^ (py << n);
                let mut row = vec![0u64; chunks];
                for (col, &m) in masks.iter().enumerate() {
                    if (idx & m) == m {
                        row[col / 64] |= 1u64 << (col % 64);
                    }
                }
                if ((phase_word & phase_mask).count_ones() & 1) != 0 {
                    row[cols / 64] |= 1u64 << (cols % 64);
                }
                rows.push(row);
            }
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn curve_support_old_point_phase_min_degree(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            curve_support_old_point_phase_has_degree_at_most(n, p, qx, qy, phase_mask, d)
        })
    }

    fn first_curve_point_u16_for_phase_test(p: u16) -> (u16, u16) {
        for x in 1..p {
            for y in 1..p {
                if is_curve_point_u16_for_phase_test(x, y, p) {
                    return (x, y);
                }
            }
        }
        panic!("toy curve had no point")
    }

    fn second_curve_point_distinct_x_u16_for_phase_test(p: u16, qx0: u16) -> (u16, u16) {
        for x in 1..p {
            if x == qx0 {
                continue;
            }
            for y in 1..p {
                if is_curve_point_u16_for_phase_test(x, y, p) {
                    return (x, y);
                }
            }
        }
        panic!("toy curve had no second point with distinct x")
    }

    fn slope_tag_retarget_phase_anf_stats(n: usize, p: u16, qa: (u16, u16), qb: (u16, u16), phase_mask: u16) -> (usize, usize) {
        // Full-domain ANF for changing the carried slope tag from being
        // relative to qa to being relative to qb:
        //   lambda_b = (qa_y - qb_y + lambda_a*(x - qa_x))/(x - qb_x).
        // This is the coordinate-conversion cost that a slope-carried
        // accumulator would pay when the next window selects a different table
        // point.  Invalid bit patterns are mapped to 0, as in the other dense
        // phase probes in this file.
        let vars = 2 * n;
        let size = 1usize << vars;
        let limb_mask = (1u16 << n) - 1;
        let mut anf = vec![0u8; size];
        for idx in 0..size {
            let x = (idx as u16) & limb_mask;
            let lambda_a = ((idx >> n) as u16) & limb_mask;
            let phase_word = if x < p && lambda_a < p && x != qb.0 {
                let y = add_mod_u16_for_phase_test(
                    qa.1,
                    mul_mod_u16_for_phase_test(lambda_a, sub_mod_u16_for_phase_test(x, qa.0, p), p),
                    p,
                );
                let lambda_b = mul_mod_u16_for_phase_test(
                    sub_mod_u16_for_phase_test(y, qb.1, p),
                    inv_mod_u16_for_phase_test(sub_mod_u16_for_phase_test(x, qb.0, p), p),
                    p,
                );
                lambda_b
            } else {
                0
            };
            anf[idx] = ((phase_word & phase_mask).count_ones() & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn slope_tag_retarget_support_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qa: (u16, u16),
        qb: (u16, u16),
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let vars = 2 * n;
        let masks = monomial_masks_for_curve_phase_test(vars, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut rows = Vec::new();
        for x in 0..p {
            if x == qa.0 || x == qb.0 {
                continue;
            }
            for y in 0..p {
                if !is_curve_point_u16_for_phase_test(x, y, p) {
                    continue;
                }
                let lambda_a = mul_mod_u16_for_phase_test(
                    sub_mod_u16_for_phase_test(y, qa.1, p),
                    inv_mod_u16_for_phase_test(sub_mod_u16_for_phase_test(x, qa.0, p), p),
                    p,
                );
                let lambda_b = mul_mod_u16_for_phase_test(
                    sub_mod_u16_for_phase_test(y, qb.1, p),
                    inv_mod_u16_for_phase_test(sub_mod_u16_for_phase_test(x, qb.0, p), p),
                    p,
                );
                let idx = (x as u32) | ((lambda_a as u32) << n);
                let mut row = vec![0u64; chunks];
                for (col, &m) in masks.iter().enumerate() {
                    if (idx & m) == m {
                        row[col / 64] |= 1u64 << (col % 64);
                    }
                }
                if ((lambda_b & phase_mask).count_ones() & 1) != 0 {
                    row[cols / 64] |= 1u64 << (cols % 64);
                }
                rows.push(row);
            }
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn slope_tag_retarget_support_phase_min_degree(
        n: usize,
        p: u16,
        qa: (u16, u16),
        qb: (u16, u16),
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            slope_tag_retarget_support_phase_has_degree_at_most(n, p, qa, qb, phase_mask, d)
        })
    }

    fn lambda_from_output_const_u16_for_phase_test(rx: u16, ry: u16, qx: u16, qy: u16, p: u16) -> Option<u16> {
        if !is_curve_point_u16_for_phase_test(rx, ry, p) {
            return None;
        }
        let dx = sub_mod_u16_for_phase_test(rx, qx, p);
        if dx == 0 {
            return None;
        }
        // From R = P + Q, Ry = -lambda*(Rx-Qx) - Qy.
        let num = if add_mod_u16_for_phase_test(ry, qy, p) == 0 {
            0
        } else {
            p - add_mod_u16_for_phase_test(ry, qy, p)
        };
        Some(mul_mod_u16_for_phase_test(num, inv_mod_u16_for_phase_test(dx, p), p))
    }

    fn endomorphism_lambda_from_output_const_u16_for_phase_test(
        rx: u16,
        ry: u16,
        qx: u16,
        qy: u16,
        p: u16,
    ) -> Option<u16> {
        if !is_curve_point_u16_for_phase_test(rx, ry, p) {
            return None;
        }
        let den = sub_mod_u16_for_phase_test(qy, ry, p);
        if den == 0 {
            return None;
        }
        let numer = add_mod_u16_for_phase_test(
            add_mod_u16_for_phase_test(
                mul_mod_u16_for_phase_test(rx, rx, p),
                mul_mod_u16_for_phase_test(qx, rx, p),
                p,
            ),
            mul_mod_u16_for_phase_test(qx, qx, p),
            p,
        );
        Some(mul_mod_u16_for_phase_test(numer, inv_mod_u16_for_phase_test(den, p), p))
    }

    fn endomorphism_curve_support_lambda_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let vars = 2 * n;
        let masks = monomial_masks_for_curve_phase_test(vars, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut rows = Vec::new();
        for rx in 0..p {
            for ry in 0..p {
                let Some(lambda) = endomorphism_lambda_from_output_const_u16_for_phase_test(rx, ry, qx, qy, p) else {
                    continue;
                };
                let idx = (rx as u32) | ((ry as u32) << n);
                let mut row = vec![0u64; chunks];
                for (col, &m) in masks.iter().enumerate() {
                    if (idx & m) == m {
                        row[col / 64] |= 1u64 << (col % 64);
                    }
                }
                if ((lambda & phase_mask).count_ones() & 1) != 0 {
                    row[cols / 64] |= 1u64 << (cols % 64);
                }
                rows.push(row);
            }
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn endomorphism_curve_support_lambda_phase_min_degree(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            endomorphism_curve_support_lambda_phase_has_degree_at_most(n, p, qx, qy, phase_mask, d)
        })
    }

    fn curve_support_lambda_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let vars = 2 * n;
        let masks = monomial_masks_for_curve_phase_test(vars, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut rows = Vec::new();
        for rx in 0..p {
            for ry in 0..p {
                let Some(lambda) = lambda_from_output_const_u16_for_phase_test(rx, ry, qx, qy, p) else {
                    continue;
                };
                let idx = (rx as u32) | ((ry as u32) << n);
                let mut row = vec![0u64; chunks];
                for (col, &m) in masks.iter().enumerate() {
                    if (idx & m) == m {
                        row[col / 64] |= 1u64 << (col % 64);
                    }
                }
                if ((lambda & phase_mask).count_ones() & 1) != 0 {
                    row[cols / 64] |= 1u64 << (cols % 64);
                }
                rows.push(row);
            }
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn curve_support_lambda_phase_min_degree(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            curve_support_lambda_phase_has_degree_at_most(n, p, qx, qy, phase_mask, d)
        })
    }

    fn curve_x_support_inverse_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qx: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let masks = monomial_masks_for_curve_phase_test(n, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut quadratic_residue = vec![false; p as usize];
        for y in 0..p {
            quadratic_residue[((y as u32 * y as u32) % p as u32) as usize] = true;
        }
        let mut rows = Vec::new();
        for x in 0..p {
            if x == qx {
                continue;
            }
            let x2 = mul_mod_u16_for_phase_test(x, x, p);
            let rhs = add_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(x2, x, p), 7, p);
            if !quadratic_residue[rhs as usize] {
                continue;
            }
            let denom = sub_mod_u16_for_phase_test(x, qx, p);
            let inv = inv_mod_u16_for_phase_test(denom, p);
            let mut row = vec![0u64; chunks];
            let idx = x as u32;
            for (col, &m) in masks.iter().enumerate() {
                if (idx & m) == m {
                    row[col / 64] |= 1u64 << (col % 64);
                }
            }
            if ((inv & phase_mask).count_ones() & 1) != 0 {
                row[cols / 64] |= 1u64 << (cols % 64);
            }
            rows.push(row);
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn curve_x_support_inverse_phase_min_degree(
        n: usize,
        p: u16,
        qx: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            curve_x_support_inverse_phase_has_degree_at_most(n, p, qx, phase_mask, d)
        })
    }

    fn curve_y_support_inverse_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qy: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let masks = monomial_masks_for_curve_phase_test(n, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut y_support = vec![false; p as usize];
        for x in 0..p {
            let x2 = mul_mod_u16_for_phase_test(x, x, p);
            let rhs = add_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(x2, x, p), 7, p);
            for y in 0..p {
                if ((y as u32 * y as u32) % p as u32) as u16 == rhs {
                    y_support[y as usize] = true;
                }
            }
        }
        let mut rows = Vec::new();
        for y in 0..p {
            if !y_support[y as usize] {
                continue;
            }
            let denom = add_mod_u16_for_phase_test(y, qy, p);
            if denom == 0 {
                continue;
            }
            let inv = inv_mod_u16_for_phase_test(denom, p);
            let mut row = vec![0u64; chunks];
            let idx = y as u32;
            for (col, &m) in masks.iter().enumerate() {
                if (idx & m) == m {
                    row[col / 64] |= 1u64 << (col % 64);
                }
            }
            if ((inv & phase_mask).count_ones() & 1) != 0 {
                row[cols / 64] |= 1u64 << (cols % 64);
            }
            rows.push(row);
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn curve_y_support_inverse_phase_min_degree(
        n: usize,
        p: u16,
        qy: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            curve_y_support_inverse_phase_has_degree_at_most(n, p, qy, phase_mask, d)
        })
    }

    fn variable_mul_old_multiplier_cleanup_anf_stats(n: usize, p: u16, phase_mask: u16) -> (usize, usize) {
        // Phase for X-measuring the old multiplier b after a hypothetical
        // in-place variable multiply has kept (a, t=a*b).  Cleaning b needs
        // b=t/a mod p, i.e. modular division in the output frame.
        let vars = 2 * n;
        let size = 1usize << vars;
        let mut anf = vec![0u8; size];
        for a in 1..p {
            let inv_a = inv_mod_u16_for_phase_test(a, p);
            for t in 0..p {
                let old_b = mul_mod_u16_for_phase_test(t, inv_a, p);
                let idx = (a as usize) | ((t as usize) << n);
                anf[idx] = ((old_b & phase_mask).count_ones() & 1) as u8;
            }
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&v| v != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn curve_x_membership_anf_stats(n: usize, p: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut quadratic_residue = vec![false; p as usize];
        for y in 0..p {
            quadratic_residue[((y as u32 * y as u32) % p as u32) as usize] = true;
        }
        let mut anf = vec![0u8; size];
        for x in 0..p {
            let x2 = mul_mod_u16_for_phase_test(x, x, p);
            let rhs = add_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(x2, x, p), 7, p);
            anf[x as usize] = quadratic_residue[rhs as usize] as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&v| v != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn pencil_slope_root_choice_anf_stats(n: usize, p: u16, qx: u16, qy: u16, phase_mask: u16) -> (usize, usize, usize) {
        let size = 1usize << n;
        let mut anf = vec![0u8; size];
        let mut supported_slopes = 0usize;
        for lambda_idx in 0..size {
            let mut phase = 0u8;
            if (lambda_idx as u16) < p {
                let lambda = lambda_idx as u16;
                let mut roots = Vec::new();
                for x in 0..p {
                    let dx = sub_mod_u16_for_phase_test(x, qx, p);
                    let y = add_mod_u16_for_phase_test(qy, mul_mod_u16_for_phase_test(lambda, dx, p), p);
                    if x != qx && is_curve_point_u16_for_phase_test(x, y, p) {
                        roots.push(x);
                    }
                }
                if roots.len() == 2 {
                    supported_slopes += 1;
                    let canonical_root = roots[0].min(roots[1]);
                    phase = ((canonical_root & phase_mask).count_ones() & 1) as u8;
                }
            }
            anf[lambda_idx] = phase;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density, supported_slopes)
    }

    #[test]
    fn endomorphism_y_denominator_support_inverse_still_grows() {
        // The j=0 endomorphism identity swaps the usual slope denominator
        // (x-Qx) for (y+Qy), with an x-only quadratic numerator.  That would be
        // much more useful if reciprocal on the curve y-coordinate support were
        // low-degree.  On secp-shaped toy primes p=1 mod 3, the y-support is
        // smaller than the x-support, but the interpolation degree still grows;
        // the endomorphism only changes which hard inverse we face.
        let cases = [
            (10usize, 1021u16, 0b10_1001_0101u16, 4usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 5usize),
            (14usize, 16381u16, 0b10_1010_0101_0101u16, 6usize),
        ];
        let mut last = 0usize;
        for &(n, p, mask, expected) in &cases {
            let (_qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = curve_y_support_inverse_phase_min_degree(n, p, qy, mask, expected)
                .expect("curve-y inverse phase should interpolate at expected threshold");
            eprintln!(
                "Curve-y-support inverse phase: n={n}, p={p}, qy={qy}, min_degree={min_degree}"
            );
            if n == 14 {
                println!("METRIC endomorphism_y_support_inv_min_degree_n14={min_degree}");
            }
            assert!(min_degree >= last);
            assert_eq!(min_degree, expected);
            last = min_degree;
        }
        assert!(last >= 6);
    }

    #[test]
    fn in_place_variable_multiply_cleanup_is_division_dense() {
        // Strategy E and several local-IMUL fantasies need an in-place variable
        // multiply near the cost of schoolbook multiplication.  If an
        // out-of-place product t=a*b is swapped into b and the old b is
        // X-measured, the phase correction from surviving (a,t) is exactly
        // b=t/a.  Toy ANFs show this cleanup is an ordinary dense modular
        // division, not a cheap kickmix side effect.
        let cases = [
            (6usize, 61u16, 0b10_1010u16),
            (8usize, 251u16, 0b1010_0101u16),
            (10usize, 1021u16, 0b10_1001_0101u16),
        ];
        for &(n, p, mask) in &cases {
            let (degree, density) = variable_mul_old_multiplier_cleanup_anf_stats(n, p, mask);
            let table = 1usize << (2 * n);
            eprintln!(
                "In-place variable-mul old-b cleanup ANF: n={n}, p={p}, degree={degree}, density={density}/{table}"
            );
            if n == 10 {
                println!("METRIC inplace_mul_cleanup_degree_n10={degree}");
                println!("METRIC inplace_mul_cleanup_density_n10={density}");
            }
            assert!(degree + 2 >= 2 * n);
            assert!(density > table / 4);
        }
    }

    #[test]
    fn curve_x_support_does_not_make_inverse_low_degree() {
        // A narrower version of the curve-support escape hatch: the affine
        // denominator dx = Px-Qx is not uniformly arbitrary; Px is an
        // x-coordinate occurring on y^2=x^3+7.  If inversion on this half-sized
        // x-support had a very low-degree phase extension, one could imagine
        // measuring quotient/inverse garbage and kickmixing it from x alone.
        // Exhaustive toy interpolation says no: the minimum degree follows the
        // coding-theory threshold for a ~2^(n-1) support set, already n/2 at
        // n=14.  This closes the "curve x-set makes reciprocal easy" variant
        // independently of the full (x,y) point-add cleanup tests.
        let cases = [
            (8usize, 251u16, 0b1010_0101u16, 4usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 5usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 6usize),
            (14usize, 16381u16, 0b10_1010_0101_0101u16, 7usize),
        ];
        let mut last = 0usize;
        for &(n, p, mask, expected) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = curve_x_support_inverse_phase_min_degree(n, p, qx, mask, expected)
                .expect("curve-x inverse phase should interpolate at expected threshold");
            eprintln!(
                "Curve-x-support inverse phase: n={n}, p={p}, q=({qx},{qy}), min_degree={min_degree}"
            );
            if n == 14 {
                println!("METRIC curve_x_support_inv_min_degree_n14={min_degree}");
            }
            assert!(min_degree >= last);
            assert_eq!(min_degree, expected);
            last = min_degree;
        }
        assert!(last >= 7);
    }

    #[test]
    fn curve_membership_oracle_itself_is_dense() {
        // Any branch decoder that tries to use "which predecessor stays on the
        // curve?" must coherently test whether x^3+7 is a square.  That
        // Legendre-symbol-style predicate is not a small local side condition:
        // on toy secp-shaped fields its full-domain ANF is maximal/near-maximal
        // degree and about half dense.  This is why the low curve-collision rate
        // cannot be counted as a free Kaliski/old-point cleanup oracle.
        let cases = [(8usize, 251u16), (10usize, 1021u16), (12usize, 4093u16), (14usize, 16381u16)];
        for &(n, p) in &cases {
            let (degree, density) = curve_x_membership_anf_stats(n, p);
            let table = 1usize << n;
            eprintln!(
                "Curve x-membership ANF: n={n}, p={p}, degree={degree}, density={density}/{table}"
            );
            if n == 14 {
                println!("METRIC curve_x_membership_degree_n14={degree}");
                println!("METRIC curve_x_membership_density_n14={density}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 4);
        }
    }

    #[test]
    fn pencil_slope_coordinate_needs_dense_root_choice_phase() {
        // A tempting coordinate change is the pencil of lines through fixed Q:
        // store the slope lambda of the line through Q and P.  But one slope
        // generally corresponds to two non-Q curve intersections; choosing the
        // right root is a square-root/discriminant problem.  A canonical root
        // bit is already dense/high-degree as a function of lambda on toy
        // secp-shaped curves, so this does not give a cheap affine point-add
        // coordinate system or MBUC cleanup.
        let cases = [
            (4usize, 13u16, 0b1010u16),
            (6usize, 61u16, 0b10_1010u16),
            (8usize, 251u16, 0b1010_0101u16),
            (10usize, 1021u16, 0b10_1001_0101u16),
            (12usize, 4093u16, 0b1010_0101_0101u16),
        ];
        for &(n, p, mask) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let (degree, density, support) = pencil_slope_root_choice_anf_stats(n, p, qx, qy, mask);
            let table = 1usize << n;
            eprintln!(
                "Pencil-slope root-choice phase: n={n}, p={p}, q=({qx},{qy}), support_slopes={support}, degree={degree}, density={density}/{table}"
            );
            assert!(degree + 1 >= n);
            assert!(density > table / 4);
        }
    }

    #[test]
    fn slope_carried_coordinate_retargeting_is_dense_division() {
        // A more structural single-inversion escape is to stop insisting that
        // the accumulator is affine.  Store a point as (x, lambda_Q), where
        // lambda_Q = (y-Q_y)/(x-Q_x) is the slope of the line from a fixed table
        // point Q to the accumulator.  For repeated addition by the SAME Q this
        // carries the hard slope instead of erasing it.  But a windowed ECDLP
        // point-add must retarget the slope channel when the next selected table
        // point is Q' instead of Q:
        //   lambda_Q' = (Q_y - Q'_y + lambda_Q*(x-Q_x))/(x-Q'_x).
        // This is a fresh variable division.  Full-domain toy ANFs are nearly
        // maximal degree/density, and even curve-support interpolation follows
        // the same growing threshold as direct lambda MBUC.  So a single slope
        // channel is not a universal 2n-bit coordinate system for Google's
        // three-lookup/window point-add architecture.
        let (qa8, qb8) = {
            let qa = first_curve_point_u16_for_phase_test(251);
            let qb = second_curve_point_distinct_x_u16_for_phase_test(251, qa.0);
            (qa, qb)
        };
        let (deg8, dens8) = slope_tag_retarget_phase_anf_stats(8, 251, qa8, qb8, 0b1010_0101);
        eprintln!(
            "slope-tag retarget full ANF n=8 qa={qa8:?} qb={qb8:?}: degree={deg8}, density={dens8}/65536"
        );
        assert_eq!(deg8, 15);
        assert_eq!(dens8, 32_320);

        let (qa10, qb10) = {
            let qa = first_curve_point_u16_for_phase_test(1021);
            let qb = second_curve_point_distinct_x_u16_for_phase_test(1021, qa.0);
            (qa, qb)
        };
        let (deg10, dens10) = slope_tag_retarget_phase_anf_stats(10, 1021, qa10, qb10, 0b10_1001_0101);
        println!("METRIC slope_tag_retarget_full_degree_n10={deg10}");
        println!("METRIC slope_tag_retarget_full_density_n10={dens10}");
        eprintln!(
            "slope-tag retarget full ANF n=10 qa={qa10:?} qb={qb10:?}: degree={deg10}, density={dens10}/1048576"
        );
        assert_eq!(deg10, 19);
        assert_eq!(dens10, 522_204);

        let support_cases = [
            (6usize, 61u16, 0b10_1010u16, 2usize),
            (8usize, 251u16, 0b1010_0101u16, 3usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 3usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 4usize),
        ];
        let mut last = 0usize;
        for &(n, p, mask, expected) in &support_cases {
            let qa = first_curve_point_u16_for_phase_test(p);
            let qb = second_curve_point_distinct_x_u16_for_phase_test(p, qa.0);
            let got = slope_tag_retarget_support_phase_min_degree(n, p, qa, qb, mask, expected)
                .expect("slope retarget support phase should interpolate by expected degree");
            eprintln!(
                "slope-tag retarget support phase: n={n}, p={p}, qa={qa:?}, qb={qb:?}, min_degree={got}"
            );
            if n == 12 {
                println!("METRIC slope_tag_retarget_support_min_degree_n12={got}");
            }
            assert!(got >= last);
            last = got;
        }
        assert!(last >= 4);
    }

    #[test]
    fn measuring_lambda_after_affine_add_still_needs_growing_degree_phase() {
        // Another way to delete the second affine inversion would be to keep the
        // computed output R, X-measure the slope lambda, and phase-correct from
        // the surviving output using lambda = -(Ry+Qy)/(Rx-Qx).  Restricted to
        // valid curve outputs this still follows the same growing-degree
        // interpolation threshold as old-point MBUC.  It is a division phase in
        // disguise, not a tiny kickmix cleanup.
        let cases = [
            (4usize, 13u16, 0b1010u16, 2usize),
            (6usize, 61u16, 0b10_1010u16, 3usize),
            (8usize, 251u16, 0b1010_0101u16, 3usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 4usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 4usize),
        ];
        let mut last_min = 0usize;
        for &(n, p, mask, expected_upper) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = curve_support_lambda_phase_min_degree(n, p, qx, qy, mask, expected_upper)
                .expect("lambda phase should interpolate by expected upper degree");
            eprintln!(
                "Support-restricted lambda phase: n={n}, p={p}, q=({qx},{qy}), min_degree={min_degree}"
            );
            assert!(min_degree >= last_min);
            last_min = min_degree;
        }
        assert!(last_min >= 4);
    }

    #[test]
    fn curve_support_mbuc_phase_still_scales_not_constant_degree() {
        // The full-domain ANF above is intentionally pessimistic: after a
        // correct point-add, the surviving output is on the elliptic curve.
        // This support-restricted interpolation asks whether that caveat saves
        // generic MBUC.  It does not look like a constant-degree/sparse phase:
        // the minimum degree follows the coding-theory dimension threshold as n
        // grows, already requiring degree 4 at n=12.  Extrapolating
        // sum_i<=d C(2n,i) >= ~2^n puts the generic real-curve extension near
        // d≈0.22n (≈56 for n=256), before any sparsity cost.
        let cases = [
            (4usize, 13u16, 0b1010u16, 2usize),
            (6usize, 61u16, 0b10_1010u16, 3usize),
            (8usize, 251u16, 0b1010_0101u16, 3usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 4usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 4usize),
        ];
        let mut last_min = 0usize;
        for &(n, p, mask, expected_upper) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = curve_support_old_point_phase_min_degree(n, p, qx, qy, mask, expected_upper)
                .expect("phase should interpolate by expected upper degree");
            eprintln!(
                "Support-restricted old-point phase: n={n}, p={p}, q=({qx},{qy}), min_degree={min_degree}"
            );
            assert!(min_degree >= last_min);
            last_min = min_degree;
        }
        assert!(last_min >= 4);
    }

    fn sequential_old_y_phase_min_degree_with_old_x_live(
        n: usize,
        p: u16,
        qx: u16,
        qy: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        let vars = 3 * n;
        assert!(vars <= 32, "test mask helper is u32-backed");
        for degree in 0..=max_degree {
            let masks = monomial_masks_for_curve_phase_test(vars, degree);
            let cols = masks.len();
            let chunks = (cols + 1 + 63) / 64;
            let mut rows = Vec::new();
            for px in 0..p {
                for py in 0..p {
                    if !is_curve_point_u16_for_phase_test(px, py, p) {
                        continue;
                    }
                    let Some((rx, ry)) = point_add_const_u16_for_phase_test(px, py, qx, qy, p) else {
                        continue;
                    };
                    let idx = (px as u32) | ((rx as u32) << n) | ((ry as u32) << (2 * n));
                    let mut row = vec![0u64; chunks];
                    for (col, &m) in masks.iter().enumerate() {
                        if (idx & m) == m {
                            row[col / 64] |= 1u64 << (col % 64);
                        }
                    }
                    if ((py & phase_mask).count_ones() & 1) != 0 {
                        row[cols / 64] |= 1u64 << (cols % 64);
                    }
                    rows.push(row);
                }
            }
            let mut rows_a = rows.clone();
            for row in &mut rows_a {
                row[cols / 64] &= !(1u64 << (cols % 64));
            }
            let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
            let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
            if rank_a == rank_aug {
                return Some(degree);
            }
        }
        None
    }

    #[test]
    fn sequential_old_coordinate_mbuc_still_has_growing_phase_degree() {
        // Try to rescue top-level MBUC by measuring old coordinates one at a
        // time.  If old x remains live when old y is X-measured, the phase only
        // needs to be a function of (old_x, R_x, R_y), not R alone.  This extra
        // n-bit side information lowers the degree threshold, but does not make
        // the phase constant-degree: support interpolation still grows from
        // degree 1 to 3 over tiny fields.  Dimension extrapolation for 3n live
        // variables and ~2^n supported points is still around degree 49 at
        // secp256k1, i.e. not a cheap kickmix correction.
        let cases = [
            (4usize, 13u16, 0b1010u16, 1usize),
            (6usize, 61u16, 0b10_1010u16, 2usize),
            (8usize, 251u16, 0b1010_0101u16, 2usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 3usize),
        ];
        let mut last = 0usize;
        for &(n, p, mask, expected) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let got = sequential_old_y_phase_min_degree_with_old_x_live(n, p, qx, qy, mask, expected)
                .expect("phase should interpolate by expected degree");
            eprintln!(
                "Sequential old-y MBUC with old-x live: n={n}, p={p}, q=({qx},{qy}), min_degree={got}"
            );
            assert!(got >= last);
            last = got;
        }
        assert!(last >= 3);
    }

    #[test]
    fn montgomery_q_history_phase_in_output_frame_is_dense_dead() {
        // The promising sparse q-history phase above is in the (x, old-y)
        // frame.  For an in-place multiplier after old y has been replaced by
        // product z, MBUC phase correction must be a function of (x,z).  In that
        // output frame it collapses back to quotient-like dense inversion.
        let (degree, density) = montgomery_q_history_phase_anf_stats_xz(8, 251, 0b1010_0101);
        eprintln!("Montgomery q-history output-frame ANF: degree={degree}, density={density}/65536");
        assert!(degree >= 14);
        assert!(density > 20_000);
    }

    #[test]
    fn montgomery_q_history_phase_growth_is_not_obviously_exponential_dense() {
        // Scaling probe for the q-history MBUC lead.  The full truth table is
        // still exponential, so this is only for n<=10.  The useful signal is
        // that densities stay far below half the table, unlike z/x quotient
        // cleanup.  The bad signal is that degree reaches the maximum each
        // time.  Treat this as "structured but not solved".
        let cases = [
            (4usize, 13u16, 0b1010u16),
            (6usize, 61u16, 0b10_1010u16),
            (8usize, 251u16, 0b1010_0101u16),
            (10usize, 1021u16, 0b10_1001_0101u16),
        ];
        let mut last_density = 0usize;
        for &(n, p, mask) in &cases {
            let (degree, density) = montgomery_q_history_phase_anf_stats(n, p, mask);
            let table = 1usize << (2 * n);
            eprintln!(
                "Montgomery q-history phase growth: n={n}, degree={degree}, density={density}/{table}"
            );
            assert!(degree >= 2 * n - 1);
            assert!(density < table / 4, "q-history became quotient-like dense at n={n}");
            assert!(density >= last_density);
            last_density = density;
        }
    }

    fn destructive_montgomery_reverse_frontier_sizes(final_t: u64, a: u64, n: usize, p: u64) -> Vec<usize> {
        let mut states = std::collections::BTreeSet::from([final_t]);
        let mut sizes = Vec::new();
        for _ in (0..n).rev() {
            let mut prev = std::collections::BTreeSet::new();
            for &next_t in &states {
                for bit in 0..=1u64 {
                    for q in 0..=1u64 {
                        let old = 2i128 * next_t as i128 - bit as i128 * a as i128 - q as i128 * p as i128;
                        if !(0..(2 * p as i128)).contains(&old) {
                            continue;
                        }
                        let old = old as u64;
                        if ((old + bit * a) & 1) == q {
                            prev.insert(old);
                        }
                    }
                }
            }
            states = prev;
            sizes.push(states.len());
        }
        sizes
    }

    #[test]
    fn efficient_curve_model_transforms_need_missing_torsion() {
        // Another architectural escape would move secp256k1 into a model with
        // cheaper complete addition laws (Montgomery/Edwards/Hessian), do the
        // point-add there, then convert back.  Over the base field, birational
        // maps preserve rational torsion.  Montgomery/Edwards models require a
        // rational 2-torsion point; Hessian/twisted-Hessian models require a
        // rational 3-torsion point.  secp256k1 has prime odd order not
        // divisible by 3, so these base-field model changes are unavailable for
        // the exact affine benchmark.
        let order = U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap();
        let two = U256::from(2u64);
        let three = U256::from(3u64);
        eprintln!(
            "secp256k1 order torsion check: order mod 2 = {}, order mod 3 = {}",
            order % two,
            order % three
        );
        assert_eq!(order % two, U256::from(1u64));
        assert_eq!(order % three, U256::from(1u64));
    }

    fn sqrt_phase_anf_stats_for_lambda_cleanup_test(n: usize, p: u16, mask: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut sqrt = vec![0u16; p as usize];
        let mut seen = vec![false; p as usize];
        for y in 0..p {
            let a = mul_mod_u16_for_phase_test(y, y, p) as usize;
            let neg_y = if y == 0 { 0 } else { p - y };
            let canonical = y.min(neg_y);
            if !seen[a] || canonical < sqrt[a] {
                sqrt[a] = canonical;
                seen[a] = true;
            }
        }
        let mut anf = vec![0u8; size];
        for x in 0..size {
            let y = if x < p as usize { sqrt[x] } else { 0 };
            anf[x] = ((y & mask).count_ones() & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    #[test]
    fn lambda_square_cleanup_would_require_dense_sqrt_phase() {
        // Another one-Kaliski escape: preserve enough old denominator data so
        // that after Rx is known we know λ² = Rx + dx + 2Qx, then recover λ
        // by a square root instead of a second division.  On p≡3 mod 4 this is
        // an exponentiation; as a Boolean phase/function it is already dense on
        // toy fields.  This is not the missing low-cost cleanup.
        let cases = [
            (8usize, 251u16, 0b1010_0101u16),
            (10usize, 1021u16, 0b10_1001_0101u16),
            (12usize, 4093u16, 0b1010_0101_0101u16),
        ];
        for &(n, p, mask) in &cases {
            let (degree, density) = sqrt_phase_anf_stats_for_lambda_cleanup_test(n, p, mask);
            let table = 1usize << n;
            eprintln!(
                "canonical sqrt phase for lambda cleanup: n={n}, p={p}, degree={degree}, density={density}/{table}"
            );
            assert!(degree + 1 >= n);
            assert!(density > table / 3);
        }
    }

    #[test]
    fn destructive_montgomery_reverse_trellis_needs_field_sized_state() {
        // Global uniqueness of y from (x,z) does not by itself make the
        // destructive Montgomery product reversible with small local state.
        // Reverse-stepping the recurrence without the consumed y_i/q_i history
        // creates an almost full [0,2p) frontier on tiny fields.  Enforcing the
        // final t0=0 condition would require a nonlocal search/quotient oracle,
        // i.e. the dense cleanup already killed above.
        let cases = [
            (8usize, 251u64, 125u64, 178u64, 183u64),
            (10usize, 1021u64, 238u64, 280u64, 432u64),
            (12usize, 4093u64, 2899u64, 1154u64, 3217u64),
        ];
        for &(n, p, a, b, expected_final) in &cases {
            let final_t = destructive_montgomery_block(0, a, b, n, p);
            assert_eq!(final_t, expected_final);
            let sizes = destructive_montgomery_reverse_frontier_sizes(final_t, a, n, p);
            let max_frontier = *sizes.iter().max().unwrap();
            eprintln!(
                "destructive Montgomery reverse frontier: n={n}, p={p}, max={max_frontier}, sizes_from_output={sizes:?}"
            );
            assert!(max_frontier >= (2 * p - 2) as usize);
        }
    }

    fn u256_bit_len(mut x: U256) -> usize {
        let mut n = 0usize;
        while !x.is_zero() {
            x >>= 1;
            n += 1;
        }
        n
    }

    fn u512_from_u256_for_halfgcd_test(x: U256) -> U512 {
        let l = x.as_limbs();
        U512::from_limbs([l[0], l[1], l[2], l[3], 0, 0, 0, 0])
    }

    fn u512_bit_len_for_halfgcd_test(x: U512) -> usize {
        if x.is_zero() { 0 } else { 512 - x.leading_zeros() as usize }
    }

    #[derive(Clone, Copy)]
    struct SignedMagU512ForHalfGcdTest {
        neg: bool,
        mag: U512,
    }

    fn smag_for_halfgcd_test(neg: bool, mag: U512) -> SignedMagU512ForHalfGcdTest {
        SignedMagU512ForHalfGcdTest { neg: neg && !mag.is_zero(), mag }
    }

    fn signed_add_for_halfgcd_test(
        a: SignedMagU512ForHalfGcdTest,
        b: SignedMagU512ForHalfGcdTest,
    ) -> SignedMagU512ForHalfGcdTest {
        if a.mag.is_zero() { return b; }
        if b.mag.is_zero() { return a; }
        if a.neg == b.neg {
            smag_for_halfgcd_test(a.neg, a.mag + b.mag)
        } else if a.mag >= b.mag {
            smag_for_halfgcd_test(a.neg, a.mag - b.mag)
        } else {
            smag_for_halfgcd_test(b.neg, b.mag - a.mag)
        }
    }

    fn signed_sub_scaled_for_halfgcd_test(
        a: SignedMagU512ForHalfGcdTest,
        q: U256,
        b: SignedMagU512ForHalfGcdTest,
    ) -> SignedMagU512ForHalfGcdTest {
        let prod = smag_for_halfgcd_test(b.neg, b.mag * u512_from_u256_for_halfgcd_test(q));
        let neg_prod = smag_for_halfgcd_test(!prod.neg, prod.mag);
        signed_add_for_halfgcd_test(a, neg_prod)
    }

    fn signed_neg_for_halfgcd_test(x: SignedMagU512ForHalfGcdTest) -> SignedMagU512ForHalfGcdTest {
        smag_for_halfgcd_test(!x.neg, x.mag)
    }

    fn signed_mul_mag_for_halfgcd_test(
        x: SignedMagU512ForHalfGcdTest,
        q_neg: bool,
        q: U512,
    ) -> SignedMagU512ForHalfGcdTest {
        smag_for_halfgcd_test(x.neg ^ q_neg, x.mag * q)
    }

    fn usize_bit_len_for_payload_test(x: usize) -> usize {
        if x == 0 { 1 } else { usize::BITS as usize - x.leading_zeros() as usize }
    }

    fn plusminus_odd_gcd_payload_for_divisor(x: U256, p: U256) -> (usize, usize, usize) {
        // A from-scratch non-BY idea: strip powers of two from x, then run an
        // ordered odd GCD where each step replaces (u>=v) by the ordered pair
        // {v, (u-v)/2^k}.  Optimistically, the shift counts k are the only
        // arithmetic payload.  Reversibility still needs the direction bit
        // telling which ordered output component was v and which was the new
        // difference; without it the predecessor is often ambiguous.
        assert!(!x.is_zero());
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut shift_payload = usize_bit_len_for_payload_test(initial_twos);
        let mut direction_payload = shift_payload;
        let mut steps = 0usize;
        if u < v {
            core::mem::swap(&mut u, &mut v);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let k_bits = usize_bit_len_for_payload_test(k);
            shift_payload += k_bits;
            direction_payload += k_bits;
            if d != v {
                direction_payload += 1;
            }
            steps += 1;
            if v >= d {
                u = v;
                v = d;
            } else {
                u = d;
            }
        }
        assert_eq!(u, U512::from(1u64));
        (shift_payload, direction_payload, steps)
    }

    fn plusminus_k_sequence_for_divisor(x: U256, p: U256) -> Vec<usize> {
        assert!(!x.is_zero());
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut out = vec![initial_twos];
        if u < v {
            core::mem::swap(&mut u, &mut v);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            out.push(k);
            if v >= d {
                u = v;
                v = d;
            } else {
                u = d;
            }
        }
        out
    }

    fn plusminus_kseq_dirs_for_toy(x: u16, p: u16) -> (Vec<usize>, Vec<u8>) {
        let mut u = p as u32;
        let mut v = (x as u32) >> x.trailing_zeros();
        let mut ks = vec![x.trailing_zeros() as usize];
        let mut dirs = Vec::new();
        if u < v {
            core::mem::swap(&mut u, &mut v);
        }
        while u != v {
            let diff = u - v;
            let k = diff.trailing_zeros();
            let d = diff >> k;
            ks.push(k as usize);
            dirs.push(if v >= d { 1u8 } else { 0u8 });
            if v >= d {
                u = v;
                v = d;
            } else {
                u = d;
            }
        }
        (ks, dirs)
    }

    fn plusminus_raw_k_bits_for_toy(ks: &[usize]) -> String {
        let mut out = String::new();
        for &k in ks {
            if k == 0 {
                out.push('0');
            } else {
                out.push_str(&format!("{k:b}"));
            }
        }
        out
    }

    fn plusminus_raw_k_rank_anf_stats(n: usize, p: u16) -> (usize, usize, usize) {
        use std::collections::{BTreeMap, BTreeSet};
        let size = 1usize << n;
        let mut by_raw: BTreeMap<String, BTreeSet<(Vec<usize>, Vec<u8>)>> = BTreeMap::new();
        let mut data = Vec::new();
        for x in 1..p {
            let (ks, dirs) = plusminus_kseq_dirs_for_toy(x, p);
            let raw = plusminus_raw_k_bits_for_toy(&ks);
            by_raw.entry(raw.clone()).or_default().insert((ks.clone(), dirs.clone()));
            data.push((x as usize, raw, ks, dirs));
        }
        let max_multiplicity = by_raw.values().map(|v| v.len()).max().unwrap_or(1);
        let ranked: BTreeMap<String, Vec<(Vec<usize>, Vec<u8>)>> = by_raw
            .into_iter()
            .map(|(raw, set)| (raw, set.into_iter().collect()))
            .collect();
        let mut anf = vec![0u8; size];
        for (x, raw, ks, dirs) in data {
            let entries = ranked.get(&raw).unwrap();
            let rank = entries.iter().position(|entry| entry.0 == ks && entry.1 == dirs).unwrap();
            anf[x] = (rank & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&v| v != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density, max_multiplicity)
    }

    fn plusminus_kseq_direction_sidecar_bits_for_toy(p: u16) -> (usize, usize) {
        use std::collections::{BTreeMap, BTreeSet};
        let mut by_kseq: BTreeMap<Vec<usize>, BTreeSet<Vec<u8>>> = BTreeMap::new();
        for x in 1..p {
            let mut u = p as u32;
            let mut v = (x as u32) >> x.trailing_zeros();
            let mut ks = vec![x.trailing_zeros() as usize];
            let mut dirs = Vec::new();
            if u < v {
                core::mem::swap(&mut u, &mut v);
            }
            while u != v {
                let diff = u - v;
                let k = diff.trailing_zeros();
                let d = diff >> k;
                ks.push(k as usize);
                dirs.push(if v >= d { 1u8 } else { 0u8 });
                if v >= d {
                    u = v;
                    v = d;
                } else {
                    u = d;
                }
            }
            by_kseq.entry(ks).or_default().insert(dirs);
        }
        let mut bits_by_x = Vec::new();
        let mut max_bits = 0usize;
        for x in 1..p {
            let mut u = p as u32;
            let mut v = (x as u32) >> x.trailing_zeros();
            let mut ks = vec![x.trailing_zeros() as usize];
            if u < v {
                core::mem::swap(&mut u, &mut v);
            }
            while u != v {
                let diff = u - v;
                let k = diff.trailing_zeros();
                let d = diff >> k;
                ks.push(k as usize);
                if v >= d {
                    u = v;
                    v = d;
                } else {
                    u = d;
                }
            }
            let mult = by_kseq.get(&ks).unwrap().len();
            let bits = if mult <= 1 { 0 } else { usize_bit_len_for_payload_test(mult - 1) };
            max_bits = max_bits.max(bits);
            bits_by_x.push(bits);
        }
        bits_by_x.sort_unstable();
        let p99 = bits_by_x[bits_by_x.len() * 99 / 100];
        (max_bits, p99)
    }

    fn plusminus_k_only_reverse_ambiguity_for_toy(p: u16) -> (usize, usize) {
        use std::collections::BTreeMap;
        let mut seen: BTreeMap<(u16, u16, u8), u8> = BTreeMap::new();
        let mut ambiguous: BTreeMap<(u16, u16, u8), bool> = BTreeMap::new();
        let mut total = 0usize;
        for x in 1..p {
            let mut u = p as u32;
            let mut v = (x as u32) >> x.trailing_zeros();
            if u < v {
                core::mem::swap(&mut u, &mut v);
            }
            while u != v {
                let diff = u - v;
                let k = diff.trailing_zeros();
                let d = diff >> k;
                let direction = if v >= d { 1u8 } else { 0u8 };
                let key = (v.max(d) as u16, v.min(d) as u16, k as u8);
                if let Some(&old) = seen.get(&key) {
                    if old != direction {
                        ambiguous.insert(key, true);
                    }
                } else {
                    seen.insert(key, direction);
                }
                total += 1;
                if v >= d {
                    u = v;
                    v = d;
                } else {
                    u = d;
                }
            }
        }
        let mut ambiguous_occurrences = 0usize;
        for x in 1..p {
            let mut u = p as u32;
            let mut v = (x as u32) >> x.trailing_zeros();
            if u < v {
                core::mem::swap(&mut u, &mut v);
            }
            while u != v {
                let diff = u - v;
                let k = diff.trailing_zeros();
                let d = diff >> k;
                let key = (v.max(d) as u16, v.min(d) as u16, k as u8);
                if ambiguous.contains_key(&key) {
                    ambiguous_occurrences += 1;
                }
                if v >= d {
                    u = v;
                    v = d;
                } else {
                    u = d;
                }
            }
        }
        (ambiguous_occurrences, total)
    }

    fn centered_euclid_abs_quotients_for_divisor(x: U256, p: U256) -> Vec<U512> {
        let mut u = smag_for_halfgcd_test(false, u512_from_u256_for_halfgcd_test(p));
        let mut v = smag_for_halfgcd_test(false, u512_from_u256_for_halfgcd_test(x));
        let mut out = Vec::new();
        while !v.mag.is_zero() {
            let numerator = (u.mag << 1usize) + v.mag;
            let denominator = v.mag << 1usize;
            let q = numerator / denominator;
            let q_neg = u.neg ^ v.neg;
            let qv = signed_mul_mag_for_halfgcd_test(v, q_neg, q);
            let r = signed_add_for_halfgcd_test(u, signed_neg_for_halfgcd_test(qv));
            out.push(q);
            u = v;
            v = r;
        }
        assert_eq!(u.mag, U512::from(1u64));
        out
    }

    fn half_gcd_matrix_parity_anf_stats(n: usize, p: u16, phase_mask: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut anf = vec![0u8; size];
        for x in 1..p {
            let mut u = p as i128;
            let mut v = x as i128;
            let mut a = 1i128;
            let mut b = 0i128;
            let mut c = 0i128;
            let mut d = 1i128;
            while v != 0 && ((u as u128).ilog2().max((v as u128).ilog2()) as usize + 1) > n / 2 {
                let q = u / v;
                let rem = u - q * v;
                (a, b, c, d) = (c, d, a - q * c, b - q * d);
                u = v;
                v = rem;
            }
            let word = ((a.unsigned_abs() as u16) ^ (b.unsigned_abs() as u16)
                ^ (c.unsigned_abs() as u16) ^ (d.unsigned_abs() as u16)) & phase_mask;
            anf[x as usize] = (word.count_ones() & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&v| v != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn euclid_quotients_for_divisor(x: U256, p: U256) -> Vec<U256> {
        assert!(!x.is_zero());
        let mut u = p;
        let mut v = x;
        let mut out = Vec::new();
        while !v.is_zero() {
            let q = u / v;
            let rem = u - q * v;
            out.push(q);
            u = v;
            v = rem;
        }
        assert_eq!(u, U256::from(1u64));
        out
    }

    fn replay_euclid_quotient_division(x: U256, y: U256, p: U256) -> (U256, U256, Vec<U256>) {
        let qs = euclid_quotients_for_divisor(x, p);
        // Apply the same unimodular quotient matrices to a data row.  If
        // (u,v)=(p,x) maps to (1,0), then (r,s)=(0,y) maps to
        // (y/x,0) modulo p.
        let mut r = U256::ZERO;
        let mut s = y;
        for &q in &qs {
            let old_r = r;
            r = s;
            s = sub_mod(old_r, (q % p).mul_mod(s, p), p);
        }
        (r, s, qs)
    }

    #[test]
    fn half_gcd_matrix_checkpoint_still_has_tail_history_problem() {
        // Recursive/half-GCD is the principled way to avoid a monolithic
        // quotient stream: stop Euclid near sqrt(p), store the 2x2 transform
        // matrix, then recurse on the residual pair.  The first checkpoint is
        // tantalizing because the matrix entries are only ~128 bits each.  But
        // it is not a 600-scratch local DIV primitive by itself: the residual
        // state is another ~256 bits, and even the raw bitlength of the tail
        // quotient payload puts matrix+tail above 600 before exact parsing,
        // recursive matrix multiplication, or cleanup is charged.
        let p = SECP256K1_P;
        let samples = 2048usize;
        let mut rng = 0x9e37_79b9_7f4a_7c15u64;
        let mut matrix_bits = Vec::with_capacity(samples);
        let mut residual_bits = Vec::with_capacity(samples);
        let mut matrix_plus_tail_payload = Vec::with_capacity(samples);
        let mut matrix_plus_residual = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = p;
            let mut v = x;
            let mut a = smag_for_halfgcd_test(false, U512::from(1u64));
            let mut b = smag_for_halfgcd_test(false, U512::ZERO);
            let mut c = smag_for_halfgcd_test(false, U512::ZERO);
            let mut d = smag_for_halfgcd_test(false, U512::from(1u64));
            while !v.is_zero() && u256_bit_len(u).max(u256_bit_len(v)) > 128 {
                let q = u / v;
                let rem = u - q * v;
                let na = c;
                let nb = d;
                let nc = signed_sub_scaled_for_halfgcd_test(a, q, c);
                let nd = signed_sub_scaled_for_halfgcd_test(b, q, d);
                u = v;
                v = rem;
                a = na;
                b = nb;
                c = nc;
                d = nd;
            }
            let mb = [a, b, c, d]
                .iter()
                .map(|z| u512_bit_len_for_halfgcd_test(z.mag))
                .sum::<usize>();
            let rb = u256_bit_len(u) + u256_bit_len(v);
            let mut tail_payload = 0usize;
            let mut tu = u;
            let mut tv = v;
            while !tv.is_zero() {
                let q = tu / tv;
                tail_payload += u256_bit_len(q);
                let rem = tu - q * tv;
                tu = tv;
                tv = rem;
            }
            matrix_bits.push(mb);
            residual_bits.push(rb);
            matrix_plus_residual.push(mb + rb);
            matrix_plus_tail_payload.push(mb + tail_payload);
        }
        matrix_bits.sort_unstable();
        residual_bits.sort_unstable();
        matrix_plus_residual.sort_unstable();
        matrix_plus_tail_payload.sort_unstable();
        let p99 = samples * 99 / 100;
        let matrix_p99 = matrix_bits[p99];
        let residual_p99 = residual_bits[p99];
        let matrix_residual_p99 = matrix_plus_residual[p99];
        let matrix_tail_p99 = matrix_plus_tail_payload[p99];
        eprintln!(
            "half-GCD checkpoint: matrix_p99={matrix_p99}, residual_p99={residual_p99}, matrix+residual_p99={matrix_residual_p99}, matrix+tail_raw_p99={matrix_tail_p99}"
        );
        println!("METRIC halfgcd_matrix_bits_p99={matrix_p99}");
        println!("METRIC halfgcd_residual_bits_p99={residual_p99}");
        println!("METRIC halfgcd_matrix_residual_bits_p99={matrix_residual_p99}");
        println!("METRIC halfgcd_matrix_tail_raw_bits_p99={matrix_tail_p99}");
        assert!(matrix_p99 < 540, "first half-GCD matrix should be compact enough to be tempting");
        assert!(matrix_residual_p99 > 760, "matrix plus live residual state exceeds 600 scratch");
        assert!(matrix_tail_p99 > 680, "matrix plus even raw tail payload exceeds 600 scratch");
    }

    #[test]
    fn plusminus_odd_gcd_shift_stream_fits_only_if_direction_is_free() {
        // Build a non-BY DIV candidate from first principles: keep only odd
        // residuals and replace (u>=v) by {v, (u-v)/2^k}.  The k stream is much
        // smaller than binary-GCD branch history and, by itself, would fit the
        // 600-scratch model.  But the ordered poststate does not say which row
        // was v and which row was the shifted difference.  One direction bit per
        // nontrivial step is the exact reversible payload, and that alone moves
        // the p99 scratch estimate back above 600 before delimiters or phase
        // cleanup are charged.
        let p = SECP256K1_P;
        let samples = 4096usize;
        let mut rng = 0x9175_6d1f_f00d_cafeu64;
        let mut shift_payloads = Vec::with_capacity(samples);
        let mut exact_payloads = Vec::with_capacity(samples);
        let mut steps = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let (shift_payload, exact_payload, step_count) = plusminus_odd_gcd_payload_for_divisor(x, p);
            shift_payloads.push(shift_payload);
            exact_payloads.push(exact_payload);
            steps.push(step_count);
        }
        shift_payloads.sort_unstable();
        exact_payloads.sort_unstable();
        steps.sort_unstable();
        let p99 = samples * 99 / 100;
        let shift_p99 = shift_payloads[p99];
        let exact_p99 = exact_payloads[p99];
        let steps_p99 = steps[p99];
        let (ambiguous, total) = plusminus_k_only_reverse_ambiguity_for_toy(4093);
        let ambiguity_frac = ambiguous as f64 / total as f64;
        eprintln!(
            "plus-minus odd GCD stream: shift_p99={shift_p99}, exact_p99={exact_p99}, steps_p99={steps_p99}, k-only ambiguity n12={ambiguous}/{total} ({ambiguity_frac:.3})"
        );
        println!("METRIC plusminus_shift_payload_p99={shift_p99}");
        println!("METRIC plusminus_shift_scratch_p99={}", 256 + shift_p99);
        println!("METRIC plusminus_exact_payload_p99={exact_p99}");
        println!("METRIC plusminus_exact_scratch_p99={}", 256 + exact_p99);
        println!("METRIC plusminus_steps_p99={steps_p99}");
        println!("METRIC plusminus_k_only_ambiguity_frac_n12={ambiguity_frac:.6}");
        assert!(256 + shift_p99 < 600, "k-only plus-minus stream should be tempting");
        assert!(ambiguity_frac > 0.25, "k-only reverse should not determine the direction bit");
        assert!(256 + exact_p99 > 730, "exact direction history should exceed the 600-scratch target");
    }

    #[test]
    fn plusminus_k_sequence_compresses_direction_but_not_parser() {
        // Correct the per-step direction pessimism: if the entire k-sequence is
        // known, the ordered plus-minus GCD usually fixes the direction sequence;
        // on exhaustive n=16 toys, a rank sidecar of at most two bits suffices.
        // That still does not make a usable 600-scratch primitive, because the
        // raw binary k payload is not self-delimiting.  Boundary bits or even an
        // empirical entropy code for the geometric k alphabet exceed the budget.
        use std::collections::BTreeMap;
        let (max_dir_bits, p99_dir_bits) = plusminus_kseq_direction_sidecar_bits_for_toy(65521);
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x51de_c0de_a11c_e55u64;
        let mut seqs = Vec::with_capacity(samples);
        let mut freq: BTreeMap<usize, usize> = BTreeMap::new();
        let mut total = 0usize;
        let mut raw_payloads = Vec::with_capacity(samples);
        let mut counts = Vec::with_capacity(samples);
        let mut unary_payloads = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let raw = ks.iter().map(|&k| usize_bit_len_for_payload_test(k)).sum::<usize>();
            let unary = ks.iter().sum::<usize>();
            for &k in &ks {
                *freq.entry(k).or_insert(0) += 1;
                total += 1;
            }
            counts.push(ks.len());
            raw_payloads.push(raw);
            unary_payloads.push(unary);
            seqs.push(ks);
        }
        raw_payloads.sort_unstable();
        counts.sort_unstable();
        unary_payloads.sort_unstable();
        let p99 = samples * 99 / 100;
        let raw_p99 = raw_payloads[p99];
        let count_p99 = counts[p99];
        let unary_p99 = unary_payloads[p99];
        let boundary_scratch_p99 = 256 + raw_p99 + count_p99;
        let unary_scratch_p99 = 256 + unary_p99;
        let log_total = (total as f64).log2();
        let mut entropy_lengths = Vec::with_capacity(samples);
        for ks in &seqs {
            let mut bits = 0.0f64;
            for &k in ks {
                let f = *freq.get(&k).unwrap() as f64;
                bits += log_total - f.log2();
            }
            entropy_lengths.push(bits);
        }
        entropy_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let entropy_scratch_p99 = 256.0 + entropy_lengths[p99];
        eprintln!(
            "plus-minus k-sequence: dir_max_bits_n16={max_dir_bits}, dir_p99_bits_n16={p99_dir_bits}, raw_p99={raw_p99}, count_p99={count_p99}, boundary_scratch_p99={boundary_scratch_p99}, unary_scratch_p99={unary_scratch_p99}, entropy_scratch_p99={entropy_scratch_p99:.1}"
        );
        println!("METRIC plusminus_kseq_direction_max_bits_n16={max_dir_bits}");
        println!("METRIC plusminus_kseq_direction_p99_bits_n16={p99_dir_bits}");
        println!("METRIC plusminus_kseq_boundary_scratch_p99={boundary_scratch_p99}");
        println!("METRIC plusminus_kseq_unary_scratch_p99={unary_scratch_p99}");
        println!("METRIC plusminus_kseq_entropy_scratch_p99={entropy_scratch_p99:.3}");
        assert!(max_dir_bits <= 2, "whole k-sequence should nearly determine directions on toys");
        assert!(boundary_scratch_p99 > 740, "explicit k boundaries should miss scratch");
        assert!(unary_scratch_p99 > 630, "unary/self-delimiting shifts should miss scratch");
        assert!(entropy_scratch_p99 > 630.0, "empirical entropy-coded k parser should still miss scratch");
    }

    #[test]
    fn plusminus_raw_k_rank_decoder_is_dense() {
        // The clever rescue for plus-minus is to concatenate raw binary k values
        // and store only a tiny rank among the valid parses.  Information-wise
        // this is tempting: toy max multiplicity is tiny and secp scratch would
        // fit if the rank decoder were local.  But the rank bit is a global
        // parsing function of x; its ANF is already maximal/half-dense on toy
        // fields.  This mirrors the rolling-hash Kaliski lesson: compressed
        // history without a cheap branch-pop decoder is not a reversible DIV.
        let cases = [(8usize, 251u16), (10, 1021), (12, 4093), (14, 16381)];
        for &(n, p) in &cases {
            let (degree, density, max_mult) = plusminus_raw_k_rank_anf_stats(n, p);
            let table = 1usize << n;
            eprintln!(
                "plus-minus raw-k rank ANF: n={n}, degree={degree}, density={density}/{table}, max_multiplicity={max_mult}"
            );
            if n == 14 {
                println!("METRIC plusminus_rawk_rank_degree_n14={degree}");
                println!("METRIC plusminus_rawk_rank_density_n14={density}");
                println!("METRIC plusminus_rawk_rank_max_multiplicity_n14={max_mult}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 3);
        }
    }

    #[test]
    fn centered_euclid_raw_stream_fits_but_parser_entropy_does_not() {
        // Replace ordinary Euclid by a signed/centered remainder step using the
        // nearest quotient.  Signs are derivable from the signed residual pair,
        // and the absolute quotients are much smaller: raw magnitude payload is
        // below the 600-scratch line when paired with one 256-bit data register.
        // Unfortunately that is not a self-contained reversible stream.  As
        // soon as we charge even one boundary bit per quotient, or an empirical
        // prefix-code length for the quotient alphabet, the parser state is back
        // above the budget.  This is a sharper version of the quotient-stream
        // lesson: centered quotients are interesting only if a live-state parser
        // can consume raw magnitude bits without separate denominator work regs.
        use std::collections::BTreeMap;
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0xced0_e0c1_1d5e_5eedu64;
        let mut seqs = Vec::with_capacity(samples);
        let mut freq: BTreeMap<U512, usize> = BTreeMap::new();
        let mut total = 0usize;
        let mut raw_payloads = Vec::with_capacity(samples);
        let mut counts = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let qs = centered_euclid_abs_quotients_for_divisor(x, p);
            let raw = qs.iter().map(|&q| u512_bit_len_for_halfgcd_test(q)).sum::<usize>();
            for &q in &qs {
                *freq.entry(q).or_insert(0) += 1;
                total += 1;
            }
            counts.push(qs.len());
            raw_payloads.push(raw);
            seqs.push(qs);
        }
        raw_payloads.sort_unstable();
        counts.sort_unstable();
        let p99 = samples * 99 / 100;
        let raw_p99 = raw_payloads[p99];
        let count_p99 = counts[p99];
        let boundary_scratch_p99 = 256 + raw_p99 + count_p99;
        let log_total = (total as f64).log2();
        let mut entropy_lengths = Vec::with_capacity(samples);
        for qs in &seqs {
            let mut bits = 0.0f64;
            for &q in qs {
                let f = *freq.get(&q).unwrap() as f64;
                bits += log_total - f.log2();
            }
            entropy_lengths.push(bits);
        }
        entropy_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let entropy_p99 = entropy_lengths[p99];
        let entropy_scratch_p99 = 256.0 + entropy_p99;
        eprintln!(
            "centered Euclid quotients: raw_p99={raw_p99}, count_p99={count_p99}, boundary_scratch_p99={boundary_scratch_p99}, entropy_scratch_p99={entropy_scratch_p99:.1}"
        );
        println!("METRIC centered_euclid_raw_payload_p99={raw_p99}");
        println!("METRIC centered_euclid_raw_scratch_p99={}", 256 + raw_p99);
        println!("METRIC centered_euclid_boundary_scratch_p99={boundary_scratch_p99}");
        println!("METRIC centered_euclid_entropy_scratch_p99={entropy_scratch_p99:.3}");
        assert!(256 + raw_p99 < 600, "raw centered quotient magnitudes should be tantalizing");
        assert!(boundary_scratch_p99 > 700, "boundary bits should kill self-contained raw packing");
        assert!(entropy_scratch_p99 > 690.0, "empirical prefix-code parser should still exceed scratch");
    }

    #[test]
    fn centered_euclid_live_x_parser_recompute_still_gate_dead() {
        // The only possible rescue for the raw centered stream is a parser that
        // uses the live denominator x to infer quotient boundaries.  The naive
        // reversible version recomputes the centered Euclid prefix from x around
        // every quotient use.  Even with the much smaller centered quotients,
        // the optimistic bitlength*width proxy is still millions of controlled
        // trials for one DIV, so the parser cannot be a simple live-x replay.
        let p = SECP256K1_P;
        let samples = 1024usize;
        let mut rng = 0xc0de_17e0_c1d5_eedu64;
        let mut weights = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = smag_for_halfgcd_test(false, u512_from_u256_for_halfgcd_test(p));
            let mut v = smag_for_halfgcd_test(false, u512_from_u256_for_halfgcd_test(x));
            let mut prefix = 0usize;
            let mut total = 0usize;
            while !v.mag.is_zero() {
                let q = ((u.mag << 1usize) + v.mag) / (v.mag << 1usize);
                let q_bits = u512_bit_len_for_halfgcd_test(q);
                prefix += q_bits * u512_bit_len_for_halfgcd_test(u.mag).max(1);
                total += 2 * prefix;
                let q_neg = u.neg ^ v.neg;
                let qv = signed_mul_mag_for_halfgcd_test(v, q_neg, q);
                let r = signed_add_for_halfgcd_test(u, signed_neg_for_halfgcd_test(qv));
                u = v;
                v = r;
            }
            weights.push(total);
        }
        weights.sort_unstable();
        let mean = weights.iter().sum::<usize>() as f64 / samples as f64;
        let p99 = weights[samples * 99 / 100];
        eprintln!("centered Euclid live-x parser recompute weight: mean={mean:.1}, p99={p99}");
        println!("METRIC centered_euclid_live_recompute_weight_mean={mean:.3}");
        println!("METRIC centered_euclid_live_recompute_weight_p99={p99}");
        assert!(mean > 5_000_000.0);
    }

    #[test]
    fn half_gcd_checkpoint_matrix_mbu_phase_is_dense_too() {
        // If the first half-GCD matrix is too large to carry with the tail, a
        // natural kickmix thought is to X-measure the matrix checkpoint and
        // phase-correct it from the original denominator x.  A representative
        // parity of the checkpoint matrix is already high-degree/half-dense on
        // toy fields, so matrix checkpoints are not a free MBUC object either.
        let cases = [
            (8usize, 251u16, 0b1010_0101u16),
            (10usize, 1021u16, 0b10_1001_0101u16),
            (12usize, 4093u16, 0b1010_0101_0101u16),
            (14usize, 16381u16, 0b10_1010_0101_0101u16),
        ];
        for &(n, p, mask) in &cases {
            let (degree, density) = half_gcd_matrix_parity_anf_stats(n, p, mask);
            let table = 1usize << n;
            eprintln!(
                "half-GCD checkpoint matrix parity ANF: n={n}, p={p}, degree={degree}, density={density}/{table}"
            );
            if n == 14 {
                println!("METRIC halfgcd_matrix_mbu_degree_n14={degree}");
                println!("METRIC halfgcd_matrix_mbu_density_n14={density}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 4);
        }
    }

    #[test]
    fn quotient_stream_division_is_algebraically_good_but_packing_blocks_scratch600() {
        // Ground-up DIV attempt, independent of Kaliski/BY microsteps:
        // compute the ordinary Euclidean quotient stream of (p,x), uncompute
        // the denominator pass, then replay the quotient matrices on
        // (r,s)=(0,y).  Algebraically this is beautiful: the data row ends as
        // (y/x,0), so after a swap it is exactly an in-place division.  The
        // scratch accounting is also tantalizing: denominator generation needs
        // one 256-bit u-register plus the quotient stream, and coefficient
        // replay needs one 256-bit r-register plus the same stream.
        //
        // The catch is that the quotient stream is variable-length.  The raw
        // payload is already right at the 600-scratch edge, and any reversible
        // self-delimiting/pointer scheme pushes it over.  Without a cheap
        // stack/packing primitive this becomes another hidden history channel,
        // not a local DIV circuit.
        let p = SECP256K1_P;
        let mut rng = 0x1234_5678_9abc_def0u64;
        for _ in 0..200 {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let y = rand_u256(&mut rng);
            let (q, z, _qs) = replay_euclid_quotient_division(x, y, p);
            assert!(z.is_zero());
            assert_eq!(q, y.mul_mod(x.inv_mod(p).unwrap(), p));
        }

        let samples = 4096usize;
        let mut payload_bits = Vec::with_capacity(samples);
        let mut one_boundary_bit_bits = Vec::with_capacity(samples);
        let mut counts = Vec::with_capacity(samples);
        let mut longdiv_weight = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = p;
            let mut v = x;
            let mut payload = 0usize;
            let mut count = 0usize;
            let mut weighted_trials = 0usize;
            while !v.is_zero() {
                let q = u / v;
                let q_bits = u256_bit_len(q);
                payload += q_bits;
                count += 1;
                weighted_trials += q_bits * u256_bit_len(u);
                let rem = u - q * v;
                u = v;
                v = rem;
            }
            payload_bits.push(payload);
            counts.push(count);
            // Unrealistically optimistic packing model: raw quotient payload
            // plus only one boundary bit per quotient.  A real prefix/rank
            // code costs more, so if even this misses scratch the dynamic
            // stream is not locally solved.
            one_boundary_bit_bits.push(payload + count);
            longdiv_weight.push(weighted_trials);
        }
        payload_bits.sort_unstable();
        one_boundary_bit_bits.sort_unstable();
        counts.sort_unstable();
        longdiv_weight.sort_unstable();
        let p99 = samples * 99 / 100;
        let payload_mean = payload_bits.iter().sum::<usize>() as f64 / samples as f64;
        let longdiv_mean = longdiv_weight.iter().sum::<usize>() as f64 / samples as f64;
        let payload_p99 = payload_bits[p99];
        let payload_max = *payload_bits.last().unwrap();
        let count_p99 = counts[p99];
        let one_boundary_p99 = one_boundary_bit_bits[p99];
        let longdiv_p99 = longdiv_weight[p99];
        let raw_scratch_max = 256 + payload_max;
        let one_boundary_scratch_p99 = 256 + one_boundary_p99;
        eprintln!(
            "Euclid quotient DIV stream: payload_mean={payload_mean:.1}, payload_p99={payload_p99}, payload_max={payload_max}, count_p99={count_p99}, one_boundary_p99={one_boundary_p99}, longdiv_mean={longdiv_mean:.1}, longdiv_p99={longdiv_p99}, raw_scratch_max={raw_scratch_max}, one_boundary_scratch_p99={one_boundary_scratch_p99}"
        );
        println!("METRIC euclid_div_replay_samples=200");
        println!("METRIC euclid_quotient_payload_mean_bits={payload_mean:.3}");
        println!("METRIC euclid_quotient_payload_p99_bits={payload_p99}");
        println!("METRIC euclid_quotient_payload_max_bits={payload_max}");
        println!("METRIC euclid_quotient_count_p99={count_p99}");
        println!("METRIC euclid_quotient_raw_scratch_max={raw_scratch_max}");
        println!("METRIC euclid_quotient_one_boundary_scratch_p99={one_boundary_scratch_p99}");
        println!("METRIC euclid_longdiv_weight_mean={longdiv_mean:.3}");
        println!("METRIC euclid_longdiv_weight_p99={longdiv_p99}");
        assert!(payload_p99 < 360, "payload should stay close to the 344-bit sidecar target");
        assert!(raw_scratch_max > 600, "raw quotient payload unexpectedly fits worst sampled scratch");
        assert!(one_boundary_scratch_p99 > 760, "even one-boundary-bit quotient packing unexpectedly fits scratch");
    }

    #[test]
    fn euclid_quotient_stream_entropy_also_exceeds_scratch600() {
        // Follow-up to the raw-payload quotient-stream DIV test.  The tempting
        // objection is that a clever prefix/arithmetic code could pack the
        // quotients near the raw bitlength.  But the quotient values themselves
        // have Gauss-Kuzmin-like entropy; on secp-sized samples an idealized
        // empirical arithmetic code is already ~513 bits.  With the mandatory
        // 256-bit data partner for coefficient replay, this is far beyond the
        // 600-scratch model before decoder cost, pointer state, or exact
        // worst-case tails are charged.
        use std::collections::BTreeMap;
        let p = SECP256K1_P;
        let mut rng = 0xface_feed_dead_beefu64;
        let samples = 8192usize;
        let mut seqs: Vec<Vec<U256>> = Vec::with_capacity(samples);
        let mut freq: BTreeMap<U256, usize> = BTreeMap::new();
        let mut total = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let qs = euclid_quotients_for_divisor(x, p);
            for &q in &qs {
                *freq.entry(q).or_insert(0) += 1;
                total += 1;
            }
            seqs.push(qs);
        }
        let log_total = (total as f64).log2();
        let mut ideal_lengths = Vec::with_capacity(samples);
        for qs in &seqs {
            let mut bits = 0.0f64;
            for &q in qs {
                let f = *freq.get(&q).unwrap() as f64;
                bits += log_total - f.log2();
            }
            ideal_lengths.push(bits);
        }
        ideal_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mean = ideal_lengths.iter().sum::<f64>() / samples as f64;
        let p99 = ideal_lengths[samples * 99 / 100];
        let max = *ideal_lengths.last().unwrap();
        let scratch_p99 = 256.0 + p99;
        eprintln!(
            "Euclid quotient empirical entropy: mean={mean:.1}, p99={p99:.1}, max={max:.1}, scratch_p99={scratch_p99:.1}, alphabet={}",
            freq.len()
        );
        println!("METRIC euclid_quotient_entropy_mean_bits={mean:.3}");
        println!("METRIC euclid_quotient_entropy_p99_bits={p99:.3}");
        println!("METRIC euclid_quotient_entropy_scratch_p99={scratch_p99:.3}");
        assert!(mean > 500.0);
        assert!(scratch_p99 > 770.0);
    }

    fn centered_euclid_qseq_for_toy(x: u16, p: u16) -> Vec<usize> {
        let mut u = p as i128;
        let mut v = x as i128;
        let mut out = Vec::new();
        while v != 0 {
            let q_mag = ((2 * u.unsigned_abs()) + v.unsigned_abs()) / (2 * v.unsigned_abs());
            let q_signed = if (u < 0) ^ (v < 0) { -(q_mag as i128) } else { q_mag as i128 };
            out.push(q_mag as usize);
            let rem = u - q_signed * v;
            u = v;
            v = rem;
        }
        out
    }

    fn raw_binary_bits_from_usizes_for_test(xs: &[usize]) -> String {
        let mut out = String::new();
        for &x in xs {
            if x == 0 {
                out.push('0');
            } else {
                out.push_str(&format!("{x:b}"));
            }
        }
        out
    }

    fn centered_euclid_raw_q_rank_anf_stats(n: usize, p: u16) -> (usize, usize, usize) {
        use std::collections::{BTreeMap, BTreeSet};
        let size = 1usize << n;
        let mut by_raw: BTreeMap<String, BTreeSet<Vec<usize>>> = BTreeMap::new();
        let mut data = Vec::new();
        for x in 1..p {
            let qs = centered_euclid_qseq_for_toy(x, p);
            let raw = raw_binary_bits_from_usizes_for_test(&qs);
            by_raw.entry(raw.clone()).or_default().insert(qs.clone());
            data.push((x as usize, raw, qs));
        }
        let max_multiplicity = by_raw.values().map(|v| v.len()).max().unwrap_or(1);
        let ranked: BTreeMap<String, Vec<Vec<usize>>> = by_raw
            .into_iter()
            .map(|(raw, set)| (raw, set.into_iter().collect()))
            .collect();
        let mut anf = vec![0u8; size];
        for (x, raw, qs) in data {
            let entries = ranked.get(&raw).unwrap();
            let rank = entries.iter().position(|entry| *entry == qs).unwrap();
            anf[x] = (rank & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density, max_multiplicity)
    }

    fn centered_euclid_quotient_payload_parity_anf_stats(n: usize, p: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut anf = vec![0u8; size];
        for x in 1..p {
            let mut u = p as i128;
            let mut v = x as i128;
            let mut parity = 0u8;
            while v != 0 {
                let q_mag = ((2 * u.unsigned_abs()) + v.unsigned_abs()) / (2 * v.unsigned_abs());
                parity ^= (q_mag.count_ones() as u8) & 1;
                let q_signed = if (u < 0) ^ (v < 0) { -(q_mag as i128) } else { q_mag as i128 };
                let rem = u - q_signed * v;
                u = v;
                v = rem;
            }
            anf[x as usize] = parity;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn euclid_quotient_payload_parity_anf_stats(n: usize, p: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut anf = vec![0u8; size];
        for x in 1..p {
            let mut u = p as u32;
            let mut v = x as u32;
            let mut parity = 0u8;
            while v != 0 {
                let q = u / v;
                parity ^= (q.count_ones() as u8) & 1;
                let rem = u - q * v;
                u = v;
                v = rem;
            }
            anf[x as usize] = parity;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn two_adic_inverse_correction_phase_anf_stats(n: usize, c: u16, mask: u16) -> (usize, usize) {
        let p = ((1u32 << n) - c as u32) as u16;
        let size = 1usize << n;
        let modulus = 1u32 << n;
        let mut anf = vec![0u8; size];
        for x in 1..p {
            let tz = x.trailing_zeros() as usize;
            let odd = x >> tz;
            let y0 = inv_mod_u16_for_power_two_odd(odd, n);
            let mut chosen = None;
            for rep in 0..4u32 {
                let yy = y0 as u32 + rep * modulus;
                let d = ((odd as u32) * yy % (p as u32)) as u16;
                if d != 0 {
                    let corr = inv_mod_u16_for_phase_test(d, p);
                    chosen = Some(corr);
                    break;
                }
            }
            let corr = chosen.expect("one of a few 2-adic representatives should avoid d=0");
            anf[x as usize] = ((corr & mask).count_ones() & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&v| v != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn inv_mod_u16_for_power_two_odd(a: u16, n: usize) -> u16 {
        debug_assert_eq!(a & 1, 1);
        let modulus = 1u32 << n;
        for y in (1u32..modulus).step_by(2) {
            if ((a as u32) * y) % modulus == 1 { return y as u16; }
        }
        unreachable!()
    }

    #[test]
    fn two_adic_inverse_still_needs_dense_field_correction() {
        // Another tempting inversion primitive for the pseudo-Mersenne prime
        // p=2^n-c: invert the odd part of x modulo 2^n by Hensel lifting, then
        // correct from the 2-adic inverse to the mod-p inverse.  For odd a,
        // a*y0 = 1 (mod 2^n), so over p we have a*y0 = d = 1+c*t and need a
        // second inverse d^-1 mod p.  Factoring powers of two out of even x is
        // only a known constant correction; the hard part is this dense d^-1.
        let cases = [
            (8usize, 5u16, 0b1010_0101u16),   // p=251
            (10usize, 3u16, 0b10_1001_0101u16), // p=1021
            (12usize, 3u16, 0b1010_0101_0101u16), // p=4093
        ];
        for &(n, c, mask) in &cases {
            let (degree, density) = two_adic_inverse_correction_phase_anf_stats(n, c, mask);
            let table = 1usize << n;
            eprintln!(
                "2-adic inverse correction phase: n={n}, c={c}, degree={degree}, density={density}/{table}"
            );
            if n == 12 {
                println!("METRIC two_adic_inv_correction_degree_n12={degree}");
                println!("METRIC two_adic_inv_correction_density_n12={density}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 4);
        }
    }

    #[test]
    fn live_x_recompute_of_euclid_quotients_is_gate_dead() {
        // The only remaining way for quotient-stream DIV to avoid storing a
        // self-contained ~500-bit history is to use the live denominator x as
        // side information while replaying the coefficient transform.  The
        // naive version recomputes the Euclidean prefix from x for every
        // quotient bit/step, then uncomputes it.  Even in the optimistic
        // long-division weight units from the previous test this is millions
        // of bit-width trials, far beyond a one-inversion budget before any
        // real controlled add/comparator constants are charged.
        let p = SECP256K1_P;
        let mut rng = 0x0ddc_0ffe_e15e_u64;
        let samples = 1024usize;
        let mut weights = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = p;
            let mut v = x;
            let mut prefix = 0usize;
            let mut total = 0usize;
            while !v.is_zero() {
                let q = u / v;
                prefix += u256_bit_len(q) * u256_bit_len(u);
                total += 2 * prefix; // compute and uncompute this prefix around one use.
                let rem = u - q * v;
                u = v;
                v = rem;
            }
            weights.push(total);
        }
        weights.sort_unstable();
        let mean = weights.iter().sum::<usize>() as f64 / samples as f64;
        let p99 = weights[samples * 99 / 100];
        eprintln!("live-x Euclid quotient recompute weight: mean={mean:.1}, p99={p99}");
        println!("METRIC euclid_live_x_recompute_weight_mean={mean:.3}");
        println!("METRIC euclid_live_x_recompute_weight_p99={p99}");
        assert!(mean > 8_000_000.0);
    }

    #[test]
    fn centered_euclid_raw_quotient_rank_decoder_is_dense() {
        // Same objection as the plus-minus k stream: perhaps raw concatenated
        // quotient bits plus a tiny rank among valid continued-fraction parses
        // avoids explicit boundaries.  The sidecar is indeed information-small
        // on toys, but the rank decoder is a global function of x and is already
        // high-degree/half-dense.  This closes the smart-boundary version of the
        // centered Euclid raw stream.
        let cases = [(8usize, 251u16), (10, 1021), (12, 4093), (14, 16381)];
        for &(n, p) in &cases {
            let (degree, density, max_mult) = centered_euclid_raw_q_rank_anf_stats(n, p);
            let table = 1usize << n;
            eprintln!(
                "centered Euclid raw-q rank ANF: n={n}, degree={degree}, density={density}/{table}, max_multiplicity={max_mult}"
            );
            if n == 14 {
                println!("METRIC centered_euclid_rawq_rank_degree_n14={degree}");
                println!("METRIC centered_euclid_rawq_rank_density_n14={density}");
                println!("METRIC centered_euclid_rawq_rank_max_multiplicity_n14={max_mult}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 3);
        }
    }

    #[test]
    fn mbuc_of_centered_euclid_quotient_stream_is_dense_too() {
        // Centered quotients made the raw payload tantalizingly small, so check
        // the other standard escape: measure the quotient stream and kickmix it
        // from x.  A representative parity of the centered quotient magnitudes
        // is already maximal/near-half-dense on toy fields.  Thus the centered
        // raw stream is not cheaply measurable; it requires an actual reversible
        // parser/consumer if it is ever revived.
        let cases = [(8usize, 251u16), (10, 1021), (12, 4093), (14, 16381)];
        for &(n, p) in &cases {
            let (degree, density) = centered_euclid_quotient_payload_parity_anf_stats(n, p);
            let table = 1usize << n;
            eprintln!(
                "Centered Euclid quotient payload parity ANF: n={n}, degree={degree}, density={density}/{table}"
            );
            if n == 14 {
                println!("METRIC centered_euclid_mbu_degree_n14={degree}");
                println!("METRIC centered_euclid_mbu_density_n14={density}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 3);
        }
    }

    #[test]
    fn mbuc_of_euclid_quotient_history_is_dense_too() {
        // If quotient history is too large, another standard escape is to
        // measure it and pay only the MBUC phase correction.  A representative
        // measurement mask (xor of all quotient payload bits) is already
        // essentially maximal-degree/half-dense as a function of x on toy
        // fields.  So Euclidean quotient history cannot be made cheap by
        // generic MBUC/kickmix either.
        let cases = [(8usize, 251u16), (10, 1021), (12, 4093)];
        for &(n, p) in &cases {
            let (degree, density) = euclid_quotient_payload_parity_anf_stats(n, p);
            let table = 1usize << n;
            eprintln!(
                "Euclid quotient-history payload parity ANF: n={n}, degree={degree}, density={density}/{table}"
            );
            if n == 12 {
                println!("METRIC euclid_quotient_mbu_degree_n12={degree}");
                println!("METRIC euclid_quotient_mbu_density_n12={density}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 3);
        }
    }

    #[test]
    fn mbuc_product_cleanup_phase_oracle_is_not_low_degree_on_toy_field() {
        // Another possible rescue for Strategy E: compute product into a clean
        // accumulator, X-measure the old multiplier, and apply only the MBUC
        // phase correction instead of reversibly dividing by the product
        // source.  The required phase is a known-mask bit of z/x mod p.
        // On even an 8-bit toy field this quotient phase function has almost
        // maximal algebraic degree and about half of all ANF monomials, so the
        // hoped-for cheap low-degree phase oracle is not present.
        let (degree, density) = quotient_phase_truth_table_anf_stats(8, 251, 0b1010_0101);
        eprintln!("quotient phase ANF: degree={degree}, density={density}/65536");
        assert!(degree >= 14);
        assert!(density > 30_000);
    }
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY D: one-invocation Kaliski on w = dx³, with Strategy C output
// formulas, using the REVERSIBILITY of the point-add map itself to clean
// up without a second inversion.
//
// Register schedule (aim for exactly 1 kaliski_inv invocation):
//
//    tx: Px → dx → dx → dx (preserved during inversion) → Rx (final)
//    ty: Py → dy → dy → Ry
//    Ancillary registers allocated inside:
//      w_reg      n         = dx³ (computed via mul+sq from tx)
//      (Kaliski state: u, v_w, r, s, m_hist)   ← standard
//      (inv_raw resides inside Kaliski r)
//      lam_aux    n         = -(some scaled form of Ry-ish) or a clean-up aux
//
// Because we preserve tx = dx through the whole inversion, the Kaliski
// backward uncompute of w_reg and its internal state is FREE (tx still
// holds dx, so the reverse-squaring/cubing that zeros w_reg is clean).
//
// The remaining question is: can we reversibly compute (Rx, Ry) into tx,
// ty using the Strategy C formulas, AND reverse any ancilla used, BEFORE
// the Kaliski backward? Or can the backward absorb the cleanup?
//
// Classical schedule (no scale factors — we'll add those later):
//
//   1. tx := Px - Qx = dx                     (classical)
//   2. ty := Py - Qy = dy                     (classical)
//   3. dx2_reg := tx * tx                     (sq, fresh reg)
//   4. w_reg  := dx2_reg * tx                 (mul, fresh reg: w_reg = dx³)
//   5. Kaliski forward on w_reg → inv_raw = w⁻¹ in Kaliski r
//   6. dy2_reg := ty * ty                     (sq, fresh reg: dy²)
//   7. v_reg := dy2_reg - dx2_reg * (tx + 2*Qx)    (sub+mul, fresh reg or reuse)
//      Here tx + 2*Qx = Px + Qx (classical offset).
//      So v = dy² - dx²·(Px + Qx).
//   8. Strategy C Rx: Rx - Qx = v · dx² · w⁻¹ = v · (dx²·w⁻¹) = v·dx⁻¹
//         (since dx²·w⁻¹ = dx²·dx⁻³ = dx⁻¹)
//      Wait, let me recheck. w = dx³, w⁻¹ = dx⁻³.
//      Strategy C says Rx = v·w⁻¹·dx = v·dx⁻²·dx = v·dx⁻¹. So Rx = v/dx.
//      But v/dx = (dy² - dx²·(Px+Qx))/dx = dy²/dx - dx·(Px+Qx).
//      Check: dy²/dx² = λ², so dy²/dx = λ² · dx. Rx = λ²·dx - dx·(Px+Qx) =
//      dx·(λ² - Px - Qx). Hmm, that's dx·Rx actually. Off by a factor.
//
// Let me re-derive correctly.
// Rx = λ² - Px - Qx. λ = dy/dx. λ² = dy²/dx². So Rx = dy²/dx² - Px - Qx.
// Rx - Qx = dy²/dx² - Px - 2Qx.
// We want: Rx - Qx = (dy² - dx²·Px - 2·dx²·Qx) / dx² = (dy² - dx²·(Px + 2Qx)) / dx²
// Let v = dy² - dx²·(Px + 2Qx). Then Rx - Qx = v / dx² = v · dx⁻².
// In Strategy C's own test code, Rx = v · dx⁻² (with Rx being Rx not Rx-Qx).
// Hmm, need to reconcile. Actually Strategy C's test shows Rx = v · dx_winv,
// where dx_winv = dx · w⁻¹ = dx · dx⁻³ = dx⁻². So Rx = v · dx⁻².
// That gives Rx (not Rx - Qx). OK.
//
// But wait: does Rx = v·dx⁻² match Rx = λ² - Px - Qx?
// v = dy² - dx²·(Px + Qx) (from strategy_c code).
// v·dx⁻² = dy²·dx⁻² - Px - Qx = λ² - Px - Qx. ✓ YES.
//
// So Strategy C uses: v_C = dy² - dx²·(Px+Qx), and Rx = v_C · dx⁻².
// We need dx²·(Px+Qx). Px+Qx = dx + 2Qx. So dx²·(dx + 2Qx) = dx³ + 2Qx·dx².
//
// v_C = dy² - dx³ - 2Qx·dx² = dy² - w - 2Qx·dx².
// (Using w = dx³ for convenience.)
//
// Now for Ry:
//   Ry = λ(Px - Rx) - Py. Using classical Qy: Ry - Qy = λ(Px - Rx) - Py - Qy
//     = λ(Px - Rx) - (ty + Qy + Qy) = λ(Px - Rx) - ty - 2Qy.
//   Messy. Use Strategy C's direct formula:
//   Ry = (dy·(dx²·Qx - v_C) - w·Qy) · w⁻¹
//      = dy·(dx²·Qx - v_C)·w⁻¹ - Qy
//   Verify: dx²·Qx - v_C = dx²·Qx - dy² + dx²·(Px+Qx) = dx²·(Px + 2Qx) - dy².
//   And dy·(dx²·(Px+2Qx) - dy²)·w⁻¹ = dy·(Px+2Qx)·dx²·dx⁻³ - dy³·dx⁻³
//                                   = dy·(Px+2Qx)/dx - dy³/dx³
//                                   = λ·(Px+2Qx) - λ³
//   So Ry + Qy = λ·(Px+2Qx) - λ³. And λ·(Px+2Qx) - λ³ = λ·(Px+2Qx - λ²)
//                                                    = λ·(Px+2Qx - (dy/dx)²).
//   Hmm, that should equal λ(Px - Rx) I think. Let's verify:
//   λ(Px - Rx) = λ·Px - λ·Rx. Rx = λ² - Px - Qx. So λ·Rx = λ³ - λ·Px - λ·Qx.
//   ∴ λ(Px - Rx) = λ·Px - λ³ + λ·Px + λ·Qx = 2λ·Px + λ·Qx - λ³ = λ·(2Px+Qx) - λ³.
//   And our derived formula gives: λ·(Px + 2Qx) - λ³. That's λ·(Px + 2Qx - λ²),
//   but λ(Px - Rx) = λ·(2Px + Qx - λ²).   DISAGREEMENT: (Px + 2Qx) vs (2Px + Qx).
//
//   Let me recheck strategy_c code more carefully.
//   `core = sub_mod(dx2_qx, v, p);` where dx2_qx = dx²·qx, v = dy² - dx²·(px+qx).
//   So core = dx²·qx - dy² + dx²·(px+qx) = dx²·(2qx + px) - dy².
//   Then numer = dy·core - w·qy.  Ry = numer · w⁻¹.
//   = (dy · (dx²·(2qx + px) - dy²) - dx³·qy) / dx³
//   = dy·(2qx + px)/dx - dy³/dx³ - qy
//   = λ·(2qx + px) - λ³ - qy
//   But standard curve add: Ry = λ(Px - Rx) - Py = λ·(2Px + Qx - λ²) - Py.
//                              = λ·(2Px + Qx) - λ³ - Py
//   Strategy C gives: λ·(Px + 2Qx) - λ³ - Qy. Standard: λ·(2Px + Qx) - λ³ - Py.
//   These differ in whether it's Qy or Py at the end, AND whether it's
//   (Px + 2Qx) or (2Px + Qx) in the linear term. Let me check:
//   λ·(Px + 2Qx) - Qy vs λ·(2Px + Qx) - Py.
//   Difference = λ·(2Qx - 2Px - Qx + Px) + (Py - Qy)
//              = λ·(Qx - Px) + (Py - Qy)
//              = -λ·dx + dy = -dy + dy = 0.  ✓
//   They're equal. Good. Strategy C's formula checks out.
//
// Back to the scheduling question. We compute Ry = strategy_c_Ry using:
//   - live: tx=dx, ty=dy, inv_raw = w⁻¹ (inside Kaliski r), Qx and Qy classical.
//   - needed temps: dy² (fresh n), dx² (fresh n), v (fresh n),
//                   dx²·(px+qx) (can reuse dx² after use), core (fresh n),
//                   dy·core (fresh n), dx³·qy (classical-scaled, cheap).
//   - outputs: Rx into fresh rx_reg, Ry into fresh ry_reg.
//
// Qubit cost during inside-of-Kaliski body: +4n for dx², dy², v, core,
// + 2n for rx_reg, ry_reg = +6n = 1536 qubits. Peak = Kaliski_state (1025)
// + tx (256) + ty (256) + w_reg (256) + dx² (256) + dy² (256) + v (256)
// + core (256) + rx_reg (256) + ry_reg (256) ≈ **3325 qubits**. Over cap!
//
// Need to fuse temps. For example, compute v INTO dy² register (since we
// can do `dy²_reg -= dx²·(px+qx)` in place). Compute core inside dx².
// Fold dy·core INTO a reused buffer. Multiply result by w⁻¹ directly into
// ry_reg.
//
// Optimistic inside-body peak: Kaliski state (~1025) + tx (256) + ty (256)
// + w_reg (256) + shared-scratch-reg (256) + rx_reg (256) + ry_reg (256)
// = 2561 qubits. Still high but below our current 2716.
//
// BUT: after Kaliski backward zeros Kaliski state, we still have rx_reg
// and ry_reg alive with the values we want to end up in tx, ty. How do
// we zero rx_reg and ry_reg?
//
// Option: after Kaliski backward, swap tx ↔ rx_reg and ty ↔ ry_reg (free
// in qubits, 2n CX for CNOT swap). Now tx = Rx, ty = Ry, rx_reg = dx,
// ry_reg = dy. Now we need to zero rx_reg and ry_reg. But dx = Px - Qx
// and Px is now gone from tx. So we can't reconstruct dx classically.
//
// Unless we CLASSICALLY add Qx, Qy to tx, ty before the swap to restore
// Px, Py there... no, that would re-overwrite Rx, Ry.
//
// Alternative: at end, do `rx_reg -= tx + Qx` after swap? Then rx_reg = dx -
// Rx - Qx = (Px - Qx) - Rx - Qx = Px - Rx - 2Qx. Nonzero.
//
// The problem persists: dx is an independent quantity from Rx after
// swap. Zeroing a quantum register holding dx requires knowing dx via
// some live quantum source + classical constants. The only live source
// of dx-worthy info at circuit end would be Px, which isn't there.
//
// CONCLUSION FROM CLASSICAL SCHEDULING:
//   A 1-Kaliski-invocation scaffold that keeps tx = Px and produces
//   (Rx, Ry) in fresh output registers fails to zero those output
//   registers at circuit end. It DOES NOT have the cleanup property
//   we need.
//
// The only way forward with Strategy C+fresh-output is: after computing
// (rx_reg, ry_reg), we use (Rx, Ry) to reconstruct (Px, Py) via the
// reversed add (which is algebraically a point-subtract), which again
// requires ONE inversion (of (Rx - Qx), a different quantity). That's
// back to 2 inversions.
//
// Takeaway: the classical schedule shows the 1-Kaliski-invocation goal
// is genuinely blocked by reversibility, not by our implementation. The
// two escapes left:
//   (a) Google's undisclosed lookup/windowing trick.
//   (b) A fundamentally different algebra (e.g. a representation where
//       cleanup IS free, like Montgomery form or projective coords, but
//       with exact affine output preserved).
//
// Strategy D is therefore NOT implemented. Keeping this note as a
// ground-truth artifact.
