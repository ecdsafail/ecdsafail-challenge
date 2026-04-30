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

    fn kaliski_terminal_iters_for_numeric_test(x: U256, max_iters: usize) -> usize {
        let (_hist, snaps) = crate::point_add::kaliski_classical_replay::kaliski_run(x, SECP256K1_P, max_iters);
        snaps.iter()
            .position(|&(_, _, _, _, f)| f == 0)
            .unwrap_or(max_iters)
    }

    #[test]
    fn approximate_iteration_trimming_is_not_sota_scale() {
        // Google's ZK statement allows approximate point-add correctness, so
        // check whether simply reducing Kaliski iterations could explain the
        // 2.7M/2.1M gap.  It cannot: 99%-ish thresholds are only a handful of
        // iterations below the exact 404/401 settings, giving tens of thousands
        // of Toffolis, not the needed >1M structural reduction.
        let c = curve();
        let mut rng = 0x99ac_c01d_17ee_5eedu64;
        let samples = 4096usize;
        let mut pair1_terms = Vec::with_capacity(samples);
        let mut pair2_terms = Vec::with_capacity(samples);
        while pair1_terms.len() < samples {
            let k1 = rand_u256(&mut rng);
            let k2 = rand_u256(&mut rng);
            let (px, py) = c.mul(c.gx, c.gy, k1);
            let (qx, qy) = c.mul(c.gx, c.gy, k2);
            if (px.is_zero() && py.is_zero()) || (qx.is_zero() && qy.is_zero()) || px == qx {
                continue;
            }
            let (rx, _ry) = c.add(px, py, qx, qy);
            let dx1 = sub_mod(px, qx, SECP256K1_P);
            let dx2 = sub_mod(qx, rx, SECP256K1_P);
            pair1_terms.push(kaliski_terminal_iters_for_numeric_test(dx1, 512));
            pair2_terms.push(kaliski_terminal_iters_for_numeric_test(dx2, 512));
        }
        pair1_terms.sort_unstable();
        pair2_terms.sort_unstable();
        let q = |v: &[usize], num: usize, den: usize| -> usize { v[(v.len() * num) / den] };
        let p99_1 = q(&pair1_terms, 99, 100);
        let p999_1 = q(&pair1_terms, 999, 1000);
        let max_1 = *pair1_terms.last().unwrap();
        let p99_2 = q(&pair2_terms, 99, 100);
        let p999_2 = q(&pair2_terms, 999, 1000);
        let max_2 = *pair2_terms.last().unwrap();
        let fail_at_392_1 = pair1_terms.iter().filter(|&&t| t > 392).count();
        let fail_at_392_2 = pair2_terms.iter().filter(|&&t| t > 392).count();
        let fail_at_384_1 = pair1_terms.iter().filter(|&&t| t > 384).count();
        let fail_at_384_2 = pair2_terms.iter().filter(|&&t| t > 384).count();
        let fail392 = (fail_at_392_1 + fail_at_392_2) as f64 / (2 * samples) as f64;
        let fail384 = (fail_at_384_1 + fail_at_384_2) as f64 / (2 * samples) as f64;
        eprintln!(
            "Kaliski terminal iters: pair1 p99={p99_1}, p999={p999_1}, max={max_1}; pair2 p99={p99_2}, p999={p999_2}, max={max_2}; fail392={fail392:.4}, fail384={fail384:.4}"
        );
        println!("METRIC kaliski_pair1_terminal_p99_iters={p99_1}");
        println!("METRIC kaliski_pair1_terminal_p999_iters={p999_1}");
        println!("METRIC kaliski_pair1_terminal_max_iters={max_1}");
        println!("METRIC kaliski_pair2_terminal_p99_iters={p99_2}");
        println!("METRIC kaliski_pair2_terminal_p999_iters={p999_2}");
        println!("METRIC kaliski_pair2_terminal_max_iters={max_2}");
        println!("METRIC kaliski_iter_trim_fail_frac_392={fail392:.6}");
        println!("METRIC kaliski_iter_trim_fail_frac_384={fail384:.6}");
        assert!(p99_1 > 380 && p99_2 > 380);
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

    fn curve_xy_support_inverse_phase_has_degree_at_most(
        n: usize,
        p: u16,
        qx: u16,
        phase_mask: u16,
        degree: usize,
    ) -> bool {
        let vars = 2 * n;
        let masks = monomial_masks_for_curve_phase_test(vars, degree);
        let cols = masks.len();
        let chunks = (cols + 1 + 63) / 64;
        let mut roots = vec![Vec::<u16>::new(); p as usize];
        for y in 0..p {
            roots[((y as u32 * y as u32) % p as u32) as usize].push(y);
        }
        let mut rows = Vec::new();
        for x in 0..p {
            if x == qx {
                continue;
            }
            let x2 = mul_mod_u16_for_phase_test(x, x, p);
            let rhs = add_mod_u16_for_phase_test(mul_mod_u16_for_phase_test(x2, x, p), 7, p);
            let denom = sub_mod_u16_for_phase_test(x, qx, p);
            let inv = inv_mod_u16_for_phase_test(denom, p);
            for &y in &roots[rhs as usize] {
                let idx = (x as u32) | ((y as u32) << n);
                let mut row = vec![0u64; chunks];
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
        }
        let mut rows_a = rows.clone();
        for row in &mut rows_a {
            row[cols / 64] &= !(1u64 << (cols % 64));
        }
        let rank_a = gf2_rank_bitrows_for_curve_phase_test(&mut rows_a, cols);
        let rank_aug = gf2_rank_bitrows_for_curve_phase_test(&mut rows, cols + 1);
        rank_a == rank_aug
    }

    fn curve_xy_support_inverse_phase_min_degree(
        n: usize,
        p: u16,
        qx: u16,
        phase_mask: u16,
        max_degree: usize,
    ) -> Option<usize> {
        (0..=max_degree).find(|&d| {
            curve_xy_support_inverse_phase_has_degree_at_most(n, p, qx, phase_mask, d)
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
    fn keeping_curve_y_live_only_moves_inverse_to_support_interpolation() {
        // Stronger objection than the x-support test: in point addition the old
        // y-coordinate is still live, so maybe the cleanup phase for 1/(x-Qx)
        // is a low-degree function of the full curve point (x,y), not of x
        // alone.  The support has only ~p rows in 2n variables, so interpolation
        // is easier, but the minimum degree still follows the support-dimension
        // threshold and grows with n; it is not a local selector/cleanup.
        let cases = [
            (8usize, 251u16, 0b1010_0101u16, 3usize),
            (10usize, 1021u16, 0b10_1001_0101u16, 4usize),
            (12usize, 4093u16, 0b1010_0101_0101u16, 4usize),
            (14usize, 16381u16, 0b10_1010_0101_0101u16, 5usize),
        ];
        let mut last = 0usize;
        for &(n, p, mask, max_degree) in &cases {
            let (qx, qy) = first_curve_point_u16_for_phase_test(p);
            let min_degree = curve_xy_support_inverse_phase_min_degree(n, p, qx, mask, max_degree)
                .expect("curve-(x,y) inverse phase should interpolate within toy bound");
            eprintln!(
                "Curve-(x,y)-support inverse phase: n={n}, p={p}, q=({qx},{qy}), min_degree={min_degree}"
            );
            if n == 14 {
                println!("METRIC curve_xy_support_inv_min_degree_n14={min_degree}");
            }
            assert!(min_degree >= last);
            last = min_degree;
        }
        assert!(last >= 4);
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

    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
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

    fn plusminus_scaled_coeff_width_for_divisor(x: U256, p: U256) -> (usize, usize, usize, usize, SignedMagU512ForHalfGcdTest) {
        // Keep coefficients over the integers with a global denominator 2^S:
        //   value = coeff * x / 2^S (mod p).
        // A plus-minus step d=(u-v)/2^k updates cd=cu-cv while any retained
        // old v gets multiplied by 2^k to share the new global scale.  If these
        // scaled coefficients stayed near 256 bits, unary shifts could be cheap
        // relabels instead of modular controlled halvings.  If they grow toward
        // the total shift count, the representation spends the same state the
        // unary stream saved.
        assert!(!x.is_zero());
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut scale = initial_twos;
        let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
        let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
        let mut max_bits = 1usize;
        let mut steps = 0usize;
        if u < v {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut cu, &mut cv);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
            let cv_scaled = smag_for_halfgcd_test(cv.neg, cv.mag << k);
            scale += k;
            max_bits = max_bits
                .max(u512_bit_len_for_halfgcd_test(cd.mag))
                .max(u512_bit_len_for_halfgcd_test(cv_scaled.mag));
            steps += 1;
            if v >= d {
                u = v;
                v = d;
                cu = cv_scaled;
                cv = cd;
            } else {
                u = d;
                cu = cd;
                cv = cv_scaled;
            }
        }
        assert_eq!(u, U512::from(1u64));
        (max_bits, scale, steps, initial_twos, cv)
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

    fn plusminus_toy_max_scale_steps(p: u16) -> (usize, usize) {
        let mut max_scale = 0usize;
        let mut max_steps = 0usize;
        for x in 1..p {
            let (ks, _dirs) = plusminus_kseq_dirs_for_toy(x, p);
            max_scale = max_scale.max(ks.iter().sum::<usize>());
            max_steps = max_steps.max(ks.len());
        }
        (max_scale, max_steps)
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
    fn plusminus_barrel_step_budget_kills_naive_physical_integration() {
        // The barrel repair fixes the W^2 bug but still charges two variable
        // shifts plus unary->binary conversion per step.  Before wiring 257-bit
        // lanes, check whether that physical step budget alone already exceeds
        // the Google low-qubit target.
        let step_forward_257 = 14_355usize; // from plusminus_barrel_step_extrap257_forward_ccx.
        let steps_public = 202usize; // sampled public envelope / active-chain stress scale.
        let two_div_step_only = 2usize * step_forward_257 * steps_public;
        let gap = two_div_step_only as isize - 2_700_000isize;
        eprintln!("plus-minus barrel physical step-only budget: step257={step_forward_257}, steps={steps_public}, two_div_step_only={two_div_step_only}, gap={gap}");
        println!("METRIC plusminus_barrel_budget_step257_ccx={step_forward_257}");
        println!("METRIC plusminus_barrel_budget_steps={steps_public}");
        println!("METRIC plusminus_barrel_budget_two_div_step_only={two_div_step_only}");
        println!("METRIC plusminus_barrel_budget_step_only_gap_to_2700k={gap}");
        assert!(gap > 0, "naive physical barrel shifts would not kill the budget; revisit");
    }

    #[test]
    fn plusminus_active_chain_generator_makes_conservative_bound_fit() {
        // The active prefix chain already equals the unary history bits, so the
        // separate unary-output ANDs in the first generator model were double
        // counting.  Recompute the conservative S<=512, steps<=256 budget with
        // the active-chain generator.
        let gen_ccx = trailing_zero_active_chain_cost_for_plusminus(256);
        let (scale_dp, _chunks) = solinas_history_carry_scale_dp_for_plusminus(512);
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
        let one_div = 2 * (256 * cint_add_ccx + 512 * cshift_ccx) + 256 * step_tax + scale_dp[512];
        let projected = 642_716usize + 2 * one_div;
        let gap = projected as isize - 2_700_000isize;
        eprintln!("plus-minus active-chain conservative budget: gen={gen_ccx}, projected={projected}, gap={gap}");
        println!("METRIC plusminus_active_chain_generator_ccx={gen_ccx}");
        println!("METRIC plusminus_active_chain_conservative_projected={projected}");
        println!("METRIC plusminus_active_chain_conservative_gap={gap}");
        assert!(gap < 0, "active-chain generator still misses conservative plus-minus bound");
    }

    #[test]
    fn plusminus_conservative_bound_budget_model() {
        // If a later proof gives coarse bounds like total scale S<=2n and
        // steps<=n, does the plus-minus/Solinas model still fit?  This guards
        // against over-trusting sampled p99/max tails.
        let max_scale = 512usize;
        let (scale_dp, _chunks) = solinas_history_carry_scale_dp_for_plusminus(max_scale);
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
        let one_div = |steps: usize, scale: usize| -> usize {
            2 * (steps * cint_add_ccx + scale * cshift_ccx) + steps * step_tax + scale_dp[scale]
        };
        let conservative_one = one_div(256, 512);
        let conservative_projected = 642_716usize + 2 * conservative_one;
        let conservative_gap = conservative_projected as isize - 2_700_000isize;
        let sampled_steps_bound_one = one_div(202, 512);
        let sampled_steps_bound_projected = 642_716usize + 2 * sampled_steps_bound_one;
        let sampled_steps_bound_gap = sampled_steps_bound_projected as isize - 2_700_000isize;
        eprintln!(
            "plus-minus conservative budget: steps256_scale512_projected={conservative_projected}, gap={conservative_gap}; steps202_scale512_projected={sampled_steps_bound_projected}, gap={sampled_steps_bound_gap}"
        );
        println!("METRIC plusminus_conservative_steps256_scale512_projected={conservative_projected}");
        println!("METRIC plusminus_conservative_steps256_scale512_gap={conservative_gap}");
        println!("METRIC plusminus_conservative_steps202_scale512_projected={sampled_steps_bound_projected}");
        println!("METRIC plusminus_conservative_steps202_scale512_gap={sampled_steps_bound_gap}");
        assert!(sampled_steps_bound_gap < 0, "S<=512 alone with sampled step bound would not fit");
    }

    #[test]
    fn plusminus_toy_worstcase_scale_steps_grow_linearly() {
        // Exact exhaustive toy maxima for the ordered plus-minus k stream.  This
        // is not a secp proof, but it is the first check on whether sampled
        // S≈400 tails might hide an exponential/worst-case disaster.
        let cases = [(8usize, 251u16), (12usize, 4093u16), (16usize, 65521u16)];
        let mut last_scale = 0usize;
        for &(n, p) in &cases {
            let (max_scale, max_steps) = plusminus_toy_max_scale_steps(p);
            eprintln!("plus-minus toy worst case: n={n}, p={p}, max_scale={max_scale}, max_steps={max_steps}");
            if n == 16 {
                println!("METRIC plusminus_toy_n16_max_scale={max_scale}");
                println!("METRIC plusminus_toy_n16_max_steps={max_steps}");
            }
            assert!(max_scale >= last_scale);
            last_scale = max_scale;
        }
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
        let raw_max = *raw_payloads.last().unwrap();
        let count_max = *counts.last().unwrap();
        let unary_max = *unary_payloads.last().unwrap();
        let boundary_scratch_p99 = 256 + raw_p99 + count_p99;
        let unary_scratch_p99 = 256 + unary_p99;
        let boundary_scratch_max = 256 + raw_max + count_max;
        let unary_scratch_max = 256 + unary_max;
        let unary_over_google = unary_payloads.iter().filter(|&&u| 256 + u > 663).count();
        let unary_over_google_frac = unary_over_google as f64 / samples as f64;
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
        let entropy_scratch_max = 256.0 + *entropy_lengths.last().unwrap();
        let entropy_over_google = entropy_lengths.iter().filter(|&&e| 256.0 + e > 663.0).count();
        let entropy_over_google_frac = entropy_over_google as f64 / samples as f64;
        eprintln!(
            "plus-minus k-sequence: dir_max_bits_n16={max_dir_bits}, dir_p99_bits_n16={p99_dir_bits}, raw_p99={raw_p99}, count_p99={count_p99}, boundary_scratch_p99={boundary_scratch_p99}, unary_scratch_p99={unary_scratch_p99}, unary_scratch_max={unary_scratch_max}, entropy_scratch_p99={entropy_scratch_p99:.1}, entropy_scratch_max={entropy_scratch_max:.1}"
        );
        println!("METRIC plusminus_kseq_direction_max_bits_n16={max_dir_bits}");
        println!("METRIC plusminus_kseq_direction_p99_bits_n16={p99_dir_bits}");
        println!("METRIC plusminus_kseq_boundary_scratch_p99={boundary_scratch_p99}");
        println!("METRIC plusminus_kseq_unary_scratch_p99={unary_scratch_p99}");
        println!("METRIC plusminus_kseq_entropy_scratch_p99={entropy_scratch_p99:.3}");
        println!("METRIC plusminus_kseq_boundary_scratch_max={boundary_scratch_max}");
        println!("METRIC plusminus_kseq_unary_scratch_max={unary_scratch_max}");
        println!("METRIC plusminus_kseq_unary_over_google_frac={unary_over_google_frac:.6}");
        println!("METRIC plusminus_kseq_entropy_scratch_max={entropy_scratch_max:.3}");
        println!("METRIC plusminus_kseq_entropy_over_google_frac={entropy_over_google_frac:.6}");
        assert!(max_dir_bits <= 2, "whole k-sequence should nearly determine directions on toys");
        assert!(boundary_scratch_p99 > 740, "explicit k boundaries should miss scratch");
        assert!(unary_scratch_p99 > 630, "unary/self-delimiting shifts should miss strict-600 scratch");
        assert!(unary_scratch_max <= 663, "sampled unary parser exceeds Google-low-qubit slack; demote plus-minus again");
        assert!(entropy_scratch_p99 > 630.0, "empirical entropy-coded k parser should still miss strict-600 scratch");
    }

    #[test]
    fn plusminus_unary_google663_replay_proxy_is_sota_shaped_but_unproven() {
        // Under the strict 600q filter the self-delimiting unary k stream was
        // dead on state size alone.  Under Google's actual 663 scratch allowance
        // the sampled unary stream fits, so charge the most optimistic local
        // replay floor: one n-bit add/sub-like operation per odd-GCD step and
        // one modular shift-like operation per unary shift bit.  This is NOT a
        // circuit claim; it deliberately excludes exact direction sidecar,
        // sign/range handling, phase-clean cleanup, and worst-case proof.  It
        // only decides whether the reopened state shape is worth a real parser
        // circuit experiment.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x6630_0ddc_0ffe_e123u64;
        let mut scratches = Vec::with_capacity(samples);
        let mut one_div_proxy = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            scratches.push(256 + unary);
            // 255 CCX is the measured order for our n=256 add/sub/halve family.
            // Double it to cover a denominator-like state update plus one
            // coefficient/product channel, still before cleanup/history tax.
            one_div_proxy.push(2 * 255usize * (unary + steps));
        }
        scratches.sort_unstable();
        one_div_proxy.sort_unstable();
        let p99 = samples * 99 / 100;
        let p999 = samples * 999 / 1000;
        let scratch_p99 = scratches[p99];
        let scratch_p999 = scratches[p999];
        let scratch_max = *scratches.last().unwrap();
        let over_google = scratches.iter().filter(|&&s| s > 663).count();
        let over_google_frac = over_google as f64 / samples as f64;
        let one_div_p99 = one_div_proxy[p99];
        let one_div_max = *one_div_proxy.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let scaffold_after_div = 642_716usize; // from low-scratch DIV budget after deleting pair1 muls.
        let projected_p99 = scaffold_after_div + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus unary Google663 replay proxy: scratch_p99={scratch_p99}, scratch_p999={scratch_p999}, scratch_max={scratch_max}, over663={over_google_frac:.6}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_unary_google_scratch_p99={scratch_p99}");
        println!("METRIC plusminus_unary_google_scratch_p999={scratch_p999}");
        println!("METRIC plusminus_unary_google_scratch_max={scratch_max}");
        println!("METRIC plusminus_unary_google_over663_frac={over_google_frac:.6}");
        println!("METRIC plusminus_unary_proxy_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_unary_proxy_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_unary_proxy_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_unary_proxy_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_unary_proxy_gap_p99_to_2700k={gap_p99}");
        assert_eq!(over_google, 0, "sampled unary plus-minus stream exceeded Google 663 scratch");
        assert!(gap_p99 < 0, "optimistic unary plus-minus replay proxy is not SOTA-shaped even before cleanup tax");
    }

    #[test]
    fn plusminus_unary_controlled_parser_tax_with_existing_primitives() {
        // Now charge the first obvious quantum-control tax for a unary-scan
        // parser.  The stream bits are denominator-derived quantum history, not
        // classical compile-time constants, so each delimiter-controlled add and
        // each unary-controlled shift needs controlled modular arithmetic.  This
        // uses the existing exact controlled add/double/halve primitives as a
        // pessimistic but concrete baseline.  If this is too high, a plus-minus
        // revival needs a Solinas/signed-representation controlled-shift
        // breakthrough, not just the good state-size result above.
        let p = SECP256K1_P;
        let ccx_count = |ops: &[crate::circuit::Op]| -> usize {
            ops.iter()
                .filter(|o| matches!(o.kind, crate::circuit::OperationType::CCX | crate::circuit::OperationType::CCZ))
                .count()
        };

        let mut b_add = super::super::B::new();
        let acc = b_add.alloc_qubits(256);
        let addend = b_add.alloc_qubits(256);
        let ctrl = b_add.alloc_qubit();
        let start = b_add.ops.len();
        super::super::cmod_add_qq(&mut b_add, &acc, &addend, ctrl, p);
        let cadd_ccx = ccx_count(&b_add.ops[start..]);

        let mut b_dbl = super::super::B::new();
        let v = b_dbl.alloc_qubits(256);
        let ctrl_d = b_dbl.alloc_qubit();
        let start = b_dbl.ops.len();
        super::super::cmod_double_inplace(&mut b_dbl, &v, p, ctrl_d);
        let cdouble_ccx = ccx_count(&b_dbl.ops[start..]);

        let mut b_half = super::super::B::new();
        let v = b_half.alloc_qubits(256);
        let ctrl_h = b_half.alloc_qubit();
        let start = b_half.ops.len();
        super::super::cmod_halve_inplace(&mut b_half, &v, p, ctrl_h);
        let chalve_ccx = ccx_count(&b_half.ops[start..]);
        let cshift_ccx = cdouble_ccx.max(chalve_ccx);

        let samples = 8192usize;
        let mut rng = 0x6630_0ddc_c7a5_c0deu64;
        let mut one_div = Vec::with_capacity(samples);
        let mut scratches = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            scratches.push(256 + unary);
            // Two channels: denominator-like and coefficient-like. Direction,
            // sign, cleanup, and sidecar controls are still not included.
            one_div.push(2 * (steps * cadd_ccx + unary * cshift_ccx));
        }
        one_div.sort_unstable();
        scratches.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let scaffold_after_div = 642_716usize;
        let projected_p99 = scaffold_after_div + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        let scratch_max = *scratches.last().unwrap();
        eprintln!(
            "plus-minus unary controlled parser tax: cadd={cadd_ccx}, cdouble={cdouble_ccx}, chalve={chalve_ccx}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}, scratch_max={scratch_max}"
        );
        println!("METRIC plusminus_unary_cadd_ccx={cadd_ccx}");
        println!("METRIC plusminus_unary_cdouble_ccx={cdouble_ccx}");
        println!("METRIC plusminus_unary_chalve_ccx={chalve_ccx}");
        println!("METRIC plusminus_unary_controlled_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_unary_controlled_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_unary_controlled_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_unary_controlled_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_unary_controlled_gap_p99_to_2700k={gap_p99}");
        println!("METRIC plusminus_unary_controlled_scratch_max={scratch_max}");
        assert!(cadd_ccx > 0 && cshift_ccx > 0, "controlled primitive accounting should be nonzero");
    }

    #[test]
    fn plusminus_scaled_integer_coefficients_make_shifts_cheap_but_width_kills() {
        // Escape hatch after controlled modular shifts cost 1280 CCX: keep a
        // globally 2^S-scaled integer coefficient pair so each unary shift is a
        // cheap left shift/relabel, not a controlled modular halve.  This test
        // charges the missing state width.  The coefficients grow with the same
        // accumulated shift count that made the unary stream fit, so the cheap
        // shift representation is not automatically a 663-scratch DIV.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x51f7_6635_ca1e_d123u64;
        let mut widths = Vec::with_capacity(samples);
        let mut scales = Vec::with_capacity(samples);
        let mut steps_v = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let (w, scale, steps, _initial_twos, _final_coeff) = plusminus_scaled_coeff_width_for_divisor(x, p);
            widths.push(w);
            scales.push(scale);
            steps_v.push(steps);
        }
        widths.sort_unstable();
        scales.sort_unstable();
        steps_v.sort_unstable();
        let p99 = samples * 99 / 100;
        let p999 = samples * 999 / 1000;
        let width_p99 = widths[p99];
        let width_p999 = widths[p999];
        let width_max = *widths.last().unwrap();
        let scale_p99 = scales[p99];
        let scale_max = *scales.last().unwrap();
        let steps_p99 = steps_v[p99];
        let two_coeff_scratch_max = 2 * width_max;
        let one_den_one_coeff_scratch_max = 256 + width_max;
        eprintln!(
            "plus-minus scaled integer coeffs: width_p99={width_p99}, width_p999={width_p999}, width_max={width_max}, scale_p99={scale_p99}, scale_max={scale_max}, steps_p99={steps_p99}, two_coeff_scratch_max={two_coeff_scratch_max}"
        );
        println!("METRIC plusminus_scaled_coeff_width_p99={width_p99}");
        println!("METRIC plusminus_scaled_coeff_width_p999={width_p999}");
        println!("METRIC plusminus_scaled_coeff_width_max={width_max}");
        println!("METRIC plusminus_scaled_coeff_scale_p99={scale_p99}");
        println!("METRIC plusminus_scaled_coeff_scale_max={scale_max}");
        println!("METRIC plusminus_scaled_coeff_steps_p99={steps_p99}");
        println!("METRIC plusminus_scaled_coeff_two_coeff_scratch_max={two_coeff_scratch_max}");
        println!("METRIC plusminus_scaled_coeff_one_den_one_coeff_scratch_max={one_den_one_coeff_scratch_max}");
        assert!(width_max > 0, "scaled coefficient accounting should be nonzero");
    }

    fn smag_mod_u256_for_plusminus_test(x: SignedMagU512ForHalfGcdTest, p: U256) -> U256 {
        let p512 = u512_from_u256_for_halfgcd_test(p);
        let r512 = x.mag % p512;
        let limbs = r512.as_limbs();
        let r = U256::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3]]);
        if x.neg && !r.is_zero() { p - r } else { r }
    }

    fn two_inv_pow_u256_for_plusminus_test(p: U256, iters: usize) -> U256 {
        let two_inv = (p.wrapping_add(U256::from(1u64))) >> 1;
        let mut acc = U256::from(1u64);
        let mut base = two_inv;
        let mut e = iters as u64;
        while e > 0 {
            if (e & 1) != 0 { acc = acc.mul_mod(base, p); }
            e >>= 1;
            if e != 0 { base = base.mul_mod(base, p); }
        }
        acc
    }

    fn local_count_ccx_for_plusminus_cost(ops: &[crate::circuit::Op]) -> usize {
        ops.iter()
            .filter(|o| matches!(o.kind, crate::circuit::OperationType::CCX | crate::circuit::OperationType::CCZ))
            .count()
    }

    fn local_cswap_for_plusminus_cost(
        b: &mut super::super::B,
        ctrl: super::super::QubitId,
        a: super::super::QubitId,
        t: super::super::QubitId,
    ) {
        b.cx(t, a);
        b.ccx(ctrl, a, t);
        b.cx(t, a);
    }

    fn emit_controlled_integer_add_for_plusminus(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        subtract: bool,
    ) {
        assert_eq!(acc.len(), a.len());
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        if subtract {
            super::super::sub_nbit_qq_fast(b, &f, acc);
        } else {
            super::super::add_nbit_qq_fast(b, &f, acc);
        }
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn controlled_integer_add_cost_for_plusminus(width: usize) -> usize {
        let mut b = super::super::B::new();
        let acc = b.alloc_qubits(width);
        let a = b.alloc_qubits(width);
        let ctrl = b.alloc_qubit();
        let start = b.ops.len();
        emit_controlled_integer_add_for_plusminus(&mut b, &acc, &a, ctrl, false);
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn compare_cost_for_plusminus(width: usize) -> usize {
        let mut b = super::super::B::new();
        let a = b.alloc_qubits(width);
        let c = b.alloc_qubits(width);
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        super::super::with_lt(&mut b, &a, &c, flag, |_b| {});
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn cswap_lanes_cost_for_plusminus(widths: &[usize]) -> usize {
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let start = b.ops.len();
        for &w in widths {
            let a = b.alloc_qubits(w);
            let c = b.alloc_qubits(w);
            for i in 0..w {
                local_cswap_for_plusminus_cost(&mut b, ctrl, a[i], c[i]);
            }
        }
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn trailing_zero_unary_generator_cost_for_plusminus(width: usize) -> usize {
        // Reversible-ish prefix-zero generator floor.  active[j] means all
        // lower bits were zero; unary[j] = active[j] & !d[j].  We count forward
        // prefix plus a symmetric cleanup pass after unary bits are copied/used.
        let mut b = super::super::B::new();
        let d = b.alloc_qubits(width);
        let active = b.alloc_qubits(width + 1);
        let unary = b.alloc_qubits(width);
        b.x(active[0]);
        let start = b.ops.len();
        for j in 0..width {
            b.x(d[j]);
            b.ccx(active[j], d[j], unary[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.x(d[j]);
        }
        for j in (0..width).rev() {
            b.x(d[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.ccx(active[j], d[j], unary[j]);
            b.x(d[j]);
        }
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn emit_trailing_zero_active_chain_history_for_plusminus(
        b: &mut super::super::B,
        d: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
    ) {
        assert_eq!(active.len(), d.len() + 1);
        assert_eq!(hist.len(), d.len());
        b.x(active[0]);
        for j in 0..d.len() {
            b.x(d[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.x(d[j]);
        }
        for j in 0..d.len() {
            b.cx(active[j + 1], hist[j]);
        }
        for j in (0..d.len()).rev() {
            b.x(d[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.x(d[j]);
        }
        b.x(active[0]);
    }

    fn emit_trailing_zero_active_chain_history_controlled_for_plusminus(
        b: &mut super::super::B,
        d: &[super::super::QubitId],
        start: super::super::QubitId,
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
    ) {
        // Controlled variant for fixed-iteration loops: if start=0, produce an
        // all-zero history and leave the active chain clean.  If start=1, this
        // is exactly the ordinary trailing-zero unary generator.
        assert_eq!(active.len(), d.len() + 1);
        assert_eq!(hist.len(), d.len());
        b.cx(start, active[0]);
        for j in 0..d.len() {
            b.x(d[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.x(d[j]);
        }
        for j in 0..d.len() {
            b.cx(active[j + 1], hist[j]);
        }
        for j in (0..d.len()).rev() {
            b.x(d[j]);
            b.ccx(active[j], d[j], active[j + 1]);
            b.x(d[j]);
        }
        b.cx(start, active[0]);
    }

    fn trailing_zero_active_chain_cost_for_plusminus(width: usize) -> usize {
        // Improved observation: active[j+1] itself is the unary-one bit for
        // position j.  If the active chain is the history payload, no separate
        // unary[j] AND is needed; copy/use active[1..] by CNOT/controls, then
        // reverse the chain.  This halves the prefix Toffoli floor.
        let mut b = super::super::B::new();
        let d = b.alloc_qubits(width);
        let active = b.alloc_qubits(width + 1);
        let hist = b.alloc_qubits(width);
        let start = b.ops.len();
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn emit_controlled_left_shift_nooverflow_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        spill: super::super::QubitId,
    ) {
        // Controlled arithmetic left shift on a promised no-overflow signed
        // two's-complement value.  The swap cascade leaves the old top bit in
        // `spill`; no-overflow means old_top == new_top, so one controlled CNOT
        // clears spill.  This is a real reversible gate sequence globally; on
        // invalid inputs it leaves nonzero spill instead of silently erasing.
        for i in (0..v.len()).rev() {
            let lo = if i == 0 { spill } else { v[i - 1] };
            local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
        }
        b.ccx(ctrl, v[v.len() - 1], spill);
    }

    fn emit_controlled_left_shift_nooverflow_inverse_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        spill: super::super::QubitId,
    ) {
        b.ccx(ctrl, v[v.len() - 1], spill);
        for i in 0..v.len() {
            let lo = if i == 0 { spill } else { v[i - 1] };
            local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
        }
    }

    fn emit_controlled_right_shift_exact_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        spill: super::super::QubitId,
    ) {
        // Controlled logical right shift on a promised even nonnegative value.
        // The old low bit lands in spill; valid plus-minus d values have that
        // bit zero for every active shift, so spill remains clean.
        local_cswap_for_plusminus_cost(b, ctrl, spill, v[0]);
        for i in 1..v.len() {
            local_cswap_for_plusminus_cost(b, ctrl, v[i - 1], v[i]);
        }
    }

    fn emit_controlled_left_shift_unsigned_exact_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        spill: super::super::QubitId,
    ) {
        // Controlled logical left shift for unsigned denominator lanes.  The
        // shifted-out top bit lands in spill; valid inverse traces satisfy
        // (v << 1) < 2^W on every active shift, so that bit is zero.  This is
        // deliberately different from the signed coefficient no-overflow shift,
        // whose cleanup checks sign-bit preservation.
        for i in (0..v.len()).rev() {
            let lo = if i == 0 { spill } else { v[i - 1] };
            local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
        }
    }

    fn emit_barrel_left_shift_unsigned_exact_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
        spill: &[super::super::QubitId],
    ) {
        assert!(spill.len() >= v.len());
        for (j, &ctrl) in kbits.iter().enumerate() {
            let s = 1usize << j;
            if s >= v.len() { break; }
            for i in (0..v.len()).rev() {
                let lo = if i < s { spill[i] } else { v[i - s] };
                local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
            }
        }
    }

    fn emit_barrel_left_shift_unsigned_exact_inverse_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
        spill: &[super::super::QubitId],
    ) {
        assert!(spill.len() >= v.len());
        for (j, &ctrl) in kbits.iter().enumerate().rev() {
            let s = 1usize << j;
            if s >= v.len() { continue; }
            for i in 0..v.len() {
                let lo = if i < s { spill[i] } else { v[i - s] };
                local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
            }
        }
    }

    fn emit_barrel_left_shift_signed_nooverflow_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
        spill: &[super::super::QubitId],
    ) {
        assert!(spill.len() >= v.len());
        for (j, &ctrl) in kbits.iter().enumerate() {
            let s = 1usize << j;
            if s >= v.len() { break; }
            for i in (0..v.len()).rev() {
                let lo = if i < s { spill[i] } else { v[i - s] };
                local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
            }
            let sign = v[v.len() - 1];
            for &sp in &spill[..s] {
                b.ccx(ctrl, sign, sp);
            }
        }
    }

    fn emit_barrel_left_shift_signed_nooverflow_inverse_for_plusminus(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
        spill: &[super::super::QubitId],
    ) {
        assert!(spill.len() >= v.len());
        for (j, &ctrl) in kbits.iter().enumerate().rev() {
            let s = 1usize << j;
            if s >= v.len() { continue; }
            let sign = v[v.len() - 1];
            for &sp in spill[..s].iter().rev() {
                b.ccx(ctrl, sign, sp);
            }
            for i in 0..v.len() {
                let lo = if i < s { spill[i] } else { v[i - s] };
                local_cswap_for_plusminus_cost(b, ctrl, lo, v[i]);
            }
        }
    }

    fn emit_unary_history_to_binary_count_for_plusminus(
        b: &mut super::super::B,
        hist: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
    ) {
        for &h in hist {
            super::super::cadd_nbit_const_fast(b, kbits, U256::from(1u64), h);
        }
    }

    fn emit_unary_history_to_binary_count_inverse_for_plusminus(
        b: &mut super::super::B,
        hist: &[super::super::QubitId],
        kbits: &[super::super::QubitId],
    ) {
        for &h in hist.iter().rev() {
            super::super::csub_nbit_const_fast(b, kbits, U256::from(1u64), h);
        }
    }

    fn controlled_left_shift_cost_for_plusminus(width: usize) -> usize {
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(width);
        let spill = b.alloc_qubit();
        let ctrl = b.alloc_qubit();
        let start = b.ops.len();
        emit_controlled_left_shift_nooverflow_for_plusminus(&mut b, &v, ctrl, spill);
        local_count_ccx_for_plusminus_cost(&b.ops[start..])
    }

    fn set_slice_u512_pm<R: sha3::digest::XofReader>(
        sim: &mut crate::sim::Simulator<R>,
        qs: &[super::super::QubitId],
        val: U512,
    ) {
        for (i, &q) in qs.iter().enumerate() {
            if val.bit(i) {
                *sim.qubit_mut(q) |= 1;
            } else {
                *sim.qubit_mut(q) &= !1;
            }
        }
    }

    fn get_slice_u512_pm<R: sha3::digest::XofReader>(
        sim: &crate::sim::Simulator<R>,
        qs: &[super::super::QubitId],
    ) -> U512 {
        let mut bytes = [0u8; 64];
        for (i, &q) in qs.iter().enumerate() {
            if (sim.qubit(q) & 1) != 0 {
                bytes[i / 8] |= 1u8 << (i % 8);
            }
        }
        U512::from_le_slice(&bytes)
    }

    #[test]
    fn plusminus_controlled_integer_addsub_circuit_is_phase_clean() {
        // Second actual step primitive: controlled two's-complement integer
        // add/sub with measurement-based cleanup of the ctrl&a temporary.  This
        // is the cd=cu-cv workhorse for the scaled plus-minus coefficient lanes.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let acc = b.alloc_qubits(W);
        let a = b.alloc_qubits(W);
        let ctrl = b.alloc_qubit();
        let start_add = b.ops.len();
        emit_controlled_integer_add_for_plusminus(&mut b, &acc, &a, ctrl, false);
        let add_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_add..]);
        let start_sub = b.ops.len();
        emit_controlled_integer_add_for_plusminus(&mut b, &acc, &a, ctrl, true);
        let sub_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_sub..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        for &ctrl_val in &[false, true] {
            for x in 0u64..64u64 {
                for y in 0u64..64u64 {
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-controlled-int-addsub-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &acc, U512::from(x));
                    set_slice_u512_pm(&mut sim, &a, U512::from(y));
                    if ctrl_val { *sim.qubit_mut(ctrl) |= 1; }
                    sim.apply(&ops);
                    // We emitted add followed by sub.  If both are controlled by
                    // the same ctrl, the net logical action is identity; this
                    // checks data, ancilla cleanup, and MBU phase together.
                    assert_eq!(get_slice_u512_pm(&sim, &acc).as_limbs()[0] & mask, x, "acc changed ctrl={ctrl_val} x={x} y={y}");
                    assert_eq!(get_slice_u512_pm(&sim, &a).as_limbs()[0] & mask, y, "a changed");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase ctrl={ctrl_val} x={x} y={y}");
                }
            }
        }
        eprintln!("plus-minus controlled integer add/sub circuit: width={W}, add_ccx={add_ccx}, sub_ccx={sub_ccx}, peak={peak}");
        println!("METRIC plusminus_cint_addsub_width={W}");
        println!("METRIC plusminus_cint_add_ccx={add_ccx}");
        println!("METRIC plusminus_cint_sub_ccx={sub_ccx}");
        println!("METRIC plusminus_cint_addsub_peak_q={peak}");
        assert_eq!(add_ccx, 2 * W - 1, "controlled integer add cost drifted");
        assert_eq!(sub_ccx, 2 * W - 1, "controlled integer sub cost drifted");
    }

    #[test]
    fn plusminus_controlled_active_chain_unary_generator_circuit_is_clean() {
        // Terminal-loop building block: inactive fixed iterations must write no
        // k-history, while active iterations must produce the ordinary unary
        // trailing-zero chain.  This lets a future fixed-bound DIV loop encode
        // active=0 as an all-zero history word instead of a separate flag.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let d = b.alloc_qubits(W);
        let start_flag = b.alloc_qubit();
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let start = b.ops.len();
        emit_trailing_zero_active_chain_history_controlled_for_plusminus(&mut b, &d, start_flag, &active, &hist);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for start_val in [false, true] {
            for val in 0u64..1024u64 {
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"plusminus-controlled-active-chain-generator-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                set_slice_u512_pm(&mut sim, &d, U512::from(val));
                if start_val { *sim.qubit_mut(start_flag) |= 1; }
                sim.apply(&ops);
                let tz = if val == 0 { W } else { (val.trailing_zeros() as usize).min(W) };
                let expected = if start_val {
                    if tz == W { (1u64 << W) - 1 } else { (1u64 << tz) - 1 }
                } else {
                    0
                };
                assert_eq!(get_slice_u512_pm(&sim, &hist).as_limbs()[0] & ((1u64 << W) - 1), expected, "controlled history mismatch start={start_val} val={val}");
                assert_eq!(get_slice_u512_pm(&sim, &d).as_limbs()[0] & ((1u64 << W) - 1), val, "d changed");
                assert_eq!((sim.qubit(start_flag) & 1) != 0, start_val, "start flag changed");
                assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active chain not clean start={start_val} val={val}");
                assert_eq!(sim.global_phase() & 1, 0, "unexpected phase start={start_val} val={val}");
            }
        }
        eprintln!("plus-minus controlled active-chain unary generator: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_controlled_active_chain_width={W}");
        println!("METRIC plusminus_controlled_active_chain_ccx={ccx}");
        println!("METRIC plusminus_controlled_active_chain_peak_q={peak}");
        assert_eq!(ccx, 2 * W, "controlled active-chain generator should keep 2W CCX cost");
    }

    #[test]
    fn plusminus_active_chain_unary_generator_circuit_is_clean() {
        // Actual reversible generator for the unary shift history used by the
        // plus-minus parser.  It computes the active prefix chain, copies it to
        // history, and reverses the chain, so all work bits are clean and there
        // is no measurement phase.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let d = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let start = b.ops.len();
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for val in 0u64..1024u64 {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-active-chain-generator-circuit-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &d, U512::from(val));
            sim.apply(&ops);
            let tz = if val == 0 { W } else { (val.trailing_zeros() as usize).min(W) };
            let expected = if tz == W { (1u64 << W) - 1 } else { (1u64 << tz) - 1 };
            assert_eq!(get_slice_u512_pm(&sim, &hist).as_limbs()[0] & ((1u64 << W) - 1), expected, "unary history mismatch val={val}");
            assert_eq!(get_slice_u512_pm(&sim, &d).as_limbs()[0] & ((1u64 << W) - 1), val, "d changed");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active chain not clean val={val}");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase val={val}");
        }
        eprintln!("plus-minus active-chain unary generator circuit: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_active_chain_circuit_width={W}");
        println!("METRIC plusminus_active_chain_circuit_ccx={ccx}");
        println!("METRIC plusminus_active_chain_circuit_peak_q={peak}");
        assert_eq!(ccx, 2 * W, "active-chain generator should cost 2W CCX");
    }

    fn emit_plusminus_inplace_step_forward_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        // In-place productive plus-minus step on the promised domain u>=v and
        // d=u-v nonzero. Leaves unary k history in `hist` and ordering history
        // in `flag`; all other work clean.
        b.x(one);
        super::super::sub_nbit_qq_fast(b, v, u); // u=d
        emit_trailing_zero_active_chain_history_for_plusminus(b, u, active, hist);
        for &h in hist {
            emit_controlled_right_shift_exact_for_plusminus(b, u, h, spill); // u=d>>k
        }
        emit_controlled_integer_add_for_plusminus(b, cu, cv, one, true); // cu=cu-cv
        for &h in hist {
            emit_controlled_left_shift_nooverflow_for_plusminus(b, cv, h, spill); // cv=cv<<k
        }
        super::super::cmp_lt_into(b, u, v, flag);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        b.x(one);
    }

    fn emit_plusminus_inplace_step_inverse_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        b.x(one);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        // After unswapping, flag is exactly (u < v) for u=d>>k and v=old_v.
        super::super::cmp_lt_into(b, u, v, flag);
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_nooverflow_inverse_for_plusminus(b, cv, h, spill);
        }
        emit_controlled_integer_add_for_plusminus(b, cu, cv, one, false); // cu=cd+cv=old cu
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_unsigned_exact_for_plusminus(b, u, h, spill); // u=d
        }
        emit_trailing_zero_active_chain_history_for_plusminus(b, u, active, hist); // clear hist
        super::super::add_nbit_qq_fast(b, v, u); // u=old u
        b.x(one);
    }

    fn emit_plusminus_low_unary_any_one_into_for_test(
        b: &mut super::super::B,
        x: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        flag: super::super::QubitId,
    ) {
        assert_eq!(x.len(), hist.len());
        let hits = b.alloc_qubits(x.len());
        for i in 0..x.len() {
            b.ccx(hist[i], x[i], hits[i]);
        }
        super::super::cmp_neq_zero_into(b, &hits, flag);
        for i in (0..x.len()).rev() {
            b.ccx(hist[i], x[i], hits[i]);
        }
        b.free_vec(&hits);
    }

    fn emit_plusminus_recover_direction_from_coeff_divisibility_for_test(
        b: &mut super::super::B,
        cu: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        flag: super::super::QubitId,
    ) {
        // Toggle by the direction recovered from the ordered coefficient lanes.
        // On the odd-GCD scaled-coefficient trace, k>=1 and exactly the lane
        // descended from old cv is divisible by 2^k.  Therefore the first
        // ordered coefficient lane is divisible iff the order-swap happened.
        let bad = b.alloc_qubit();
        emit_plusminus_low_unary_any_one_into_for_test(b, cu, hist, bad);
        b.x(flag);
        b.cx(bad, flag);
        emit_plusminus_low_unary_any_one_into_for_test(b, cu, hist, bad);
        b.free(bad);
    }

    fn emit_plusminus_inplace_step_forward_konly_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        emit_plusminus_inplace_step_forward_for_test(b, u, v, cu, cv, active, hist, spill, flag, one);
        emit_plusminus_recover_direction_from_coeff_divisibility_for_test(b, cu, hist, flag);
    }

    fn emit_plusminus_inplace_step_inverse_konly_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        emit_plusminus_recover_direction_from_coeff_divisibility_for_test(b, cu, hist, flag);
        emit_plusminus_inplace_step_inverse_for_test(b, u, v, cu, cv, active, hist, spill, flag, one);
    }

    fn emit_plusminus_recover_direction_from_coeff_divisibility_controlled_for_test(
        b: &mut super::super::B,
        cu: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        flag: super::super::QubitId,
    ) {
        let bad = b.alloc_qubit();
        let div = b.alloc_qubit();
        emit_plusminus_low_unary_any_one_into_for_test(b, cu, hist, bad);
        b.x(div);
        b.cx(bad, div); // div = low k bits of cu are all zero.
        b.ccx(ctrl, div, flag);
        b.cx(bad, div);
        b.x(div);
        emit_plusminus_low_unary_any_one_into_for_test(b, cu, hist, bad);
        b.free(div);
        b.free(bad);
    }

    fn emit_plusminus_inplace_step_forward_konly_active_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active_chain: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
    ) {
        let act = b.alloc_qubit();
        super::super::sub_nbit_qq_fast(b, v, u); // u=d, possibly zero.
        super::super::cmp_neq_zero_into(b, u, act);
        emit_trailing_zero_active_chain_history_controlled_for_plusminus(b, u, act, active_chain, hist);
        for &h in hist {
            emit_controlled_right_shift_exact_for_plusminus(b, u, h, spill);
        }
        emit_controlled_integer_add_for_plusminus(b, cu, cv, act, true);
        for &h in hist {
            emit_controlled_left_shift_nooverflow_for_plusminus(b, cv, h, spill);
        }
        let cmp = b.alloc_qubit();
        super::super::cmp_lt_into(b, u, v, cmp);
        b.ccx(act, cmp, flag);
        super::super::cmp_lt_into(b, u, v, cmp);
        b.free(cmp);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        emit_plusminus_recover_direction_from_coeff_divisibility_controlled_for_test(b, cu, hist, act, flag);
        // Inactive case had u=d=0; restore u=v. Active case is left alone.
        b.x(act);
        emit_controlled_integer_add_for_plusminus(b, u, v, act, false);
        b.x(act);
        super::super::cmp_neq_zero_into(b, hist, act); // hist!=0 clears active; hist=0 leaves it clear.
        b.free(act);
    }

    fn emit_plusminus_inplace_step_inverse_konly_active_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active_chain: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: super::super::QubitId,
        flag: super::super::QubitId,
    ) {
        let act = b.alloc_qubit();
        super::super::cmp_neq_zero_into(b, hist, act);
        emit_plusminus_recover_direction_from_coeff_divisibility_controlled_for_test(b, cu, hist, act, flag);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        super::super::cmp_lt_into(b, u, v, flag);
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_nooverflow_inverse_for_plusminus(b, cv, h, spill);
        }
        emit_controlled_integer_add_for_plusminus(b, cu, cv, act, false);
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_unsigned_exact_for_plusminus(b, u, h, spill);
        }
        emit_trailing_zero_active_chain_history_controlled_for_plusminus(b, u, act, active_chain, hist);
        emit_controlled_integer_add_for_plusminus(b, u, v, act, false);
        // Clear act from the restored pre-step equality/inequality.
        super::super::sub_nbit_qq_fast(b, v, u);
        super::super::cmp_neq_zero_into(b, u, act);
        super::super::add_nbit_qq_fast(b, v, u);
        b.free(act);
    }

    fn plusminus_kbit_width_for_history_len(len: usize) -> usize {
        let mut bits = 0usize;
        let mut x = len;
        while x > 0 {
            bits += 1;
            x >>= 1;
        }
        bits.max(1)
    }

    fn emit_plusminus_inplace_step_forward_konly_barrel_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: &[super::super::QubitId],
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        b.x(one);
        super::super::sub_nbit_qq_fast(b, v, u); // u=d
        emit_trailing_zero_active_chain_history_for_plusminus(b, u, active, hist);
        let kbits = b.alloc_qubits(plusminus_kbit_width_for_history_len(hist.len()));
        emit_unary_history_to_binary_count_for_plusminus(b, hist, &kbits);
        emit_barrel_left_shift_unsigned_exact_inverse_for_plusminus(b, u, &kbits, spill); // u=d>>k
        emit_controlled_integer_add_for_plusminus(b, cu, cv, one, true);
        emit_barrel_left_shift_signed_nooverflow_for_plusminus(b, cv, &kbits, spill); // cv=cv<<k
        emit_unary_history_to_binary_count_inverse_for_plusminus(b, hist, &kbits);
        b.free_vec(&kbits);
        super::super::cmp_lt_into(b, u, v, flag);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        b.x(one);
        emit_plusminus_recover_direction_from_coeff_divisibility_for_test(b, cu, hist, flag);
    }

    fn emit_plusminus_inplace_step_inverse_konly_barrel_for_test(
        b: &mut super::super::B,
        u: &[super::super::QubitId],
        v: &[super::super::QubitId],
        cu: &[super::super::QubitId],
        cv: &[super::super::QubitId],
        active: &[super::super::QubitId],
        hist: &[super::super::QubitId],
        spill: &[super::super::QubitId],
        flag: super::super::QubitId,
        one: super::super::QubitId,
    ) {
        emit_plusminus_recover_direction_from_coeff_divisibility_for_test(b, cu, hist, flag);
        b.x(one);
        for i in 0..u.len() {
            local_cswap_for_plusminus_cost(b, flag, u[i], v[i]);
            local_cswap_for_plusminus_cost(b, flag, cu[i], cv[i]);
        }
        super::super::cmp_lt_into(b, u, v, flag);
        let kbits = b.alloc_qubits(plusminus_kbit_width_for_history_len(hist.len()));
        emit_unary_history_to_binary_count_for_plusminus(b, hist, &kbits);
        emit_barrel_left_shift_signed_nooverflow_inverse_for_plusminus(b, cv, &kbits, spill);
        emit_controlled_integer_add_for_plusminus(b, cu, cv, one, false);
        emit_barrel_left_shift_unsigned_exact_for_plusminus(b, u, &kbits, spill);
        emit_unary_history_to_binary_count_inverse_for_plusminus(b, hist, &kbits);
        b.free_vec(&kbits);
        emit_trailing_zero_active_chain_history_for_plusminus(b, u, active, hist);
        super::super::add_nbit_qq_fast(b, v, u);
        b.x(one);
    }

    fn plusminus_classical_step_mod_width_for_test(
        u: &mut u64,
        v: &mut u64,
        cu: &mut u64,
        cv: &mut u64,
        width: usize,
    ) {
        let mask = (1u64 << width) - 1;
        let diff = *u - *v;
        let k = diff.trailing_zeros() as usize;
        assert!(k > 0, "odd-GCD test domain requires k>0");
        let d = diff >> k;
        let cd = cu.wrapping_sub(*cv) & mask;
        let cvs = (*cv << k) & mask;
        if d < *v {
            *u = *v;
            *v = d;
            *cu = cvs;
            *cv = cd;
        } else {
            *u = d;
            *cu = cd;
            *cv = cvs;
        }
    }

    #[test]
    fn plusminus_barrel_shift_roundtrip_is_logwidth_not_quadratic() {
        // Candidate replacement for the W^2 repeated-shift skeleton: store or
        // compute k in binary and use a reversible barrel shift.  This validates
        // the exact unsigned denominator shift primitive; binary-k generation
        // and signed coefficient barrel shifts remain separate blockers.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        const KB: usize = 4;
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(W);
        let kbits = b.alloc_qubits(KB);
        let spill = b.alloc_qubits(W);
        let start_left = b.ops.len();
        emit_barrel_left_shift_unsigned_exact_for_plusminus(&mut b, &v, &kbits, &spill);
        let left_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_left..]);
        let start_inv = b.ops.len();
        emit_barrel_left_shift_unsigned_exact_inverse_for_plusminus(&mut b, &v, &kbits, &spill);
        let inv_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_inv..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for k in 0u64..8u64 {
            for x in 0u64..(1u64 << (W - k as usize)).min(256) {
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"plusminus-barrel-unsigned-shift-roundtrip-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                set_slice_u512_pm(&mut sim, &v, U512::from(x));
                set_slice_u512_pm(&mut sim, &kbits, U512::from(k));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & ((1u64 << W) - 1), x, "v changed k={k} x={x}");
                assert_eq!(get_slice_u512_pm(&sim, &kbits).as_limbs()[0], k, "k changed");
                assert_eq!(get_slice_u512_pm(&sim, &spill), U512::ZERO, "spill dirty k={k} x={x}");
                assert_eq!(sim.global_phase() & 1, 0, "unexpected phase k={k} x={x}");
            }
        }
        let mut b64 = super::super::B::new();
        let v64 = b64.alloc_qubits(64);
        let k64 = b64.alloc_qubits(6);
        let spill64 = b64.alloc_qubits(64);
        let s64 = b64.ops.len();
        emit_barrel_left_shift_unsigned_exact_for_plusminus(&mut b64, &v64, &k64, &spill64);
        let w64_left_ccx = local_count_ccx_for_plusminus_cost(&b64.ops[s64..]);
        let extrap257_left = ((w64_left_ccx as f64) * (257.0 / 64.0) * (8.0 / 6.0)).round() as usize;
        eprintln!("plus-minus barrel unsigned shift: w16 left={left_ccx}, inv={inv_ccx}, peak={peak}, w64_left={w64_left_ccx}, extrap257_left={extrap257_left}");
        println!("METRIC plusminus_barrel_shift_w16_left_ccx={left_ccx}");
        println!("METRIC plusminus_barrel_shift_w16_inverse_ccx={inv_ccx}");
        println!("METRIC plusminus_barrel_shift_w16_peak_q={peak}");
        println!("METRIC plusminus_barrel_shift_w64_left_ccx={w64_left_ccx}");
        println!("METRIC plusminus_barrel_shift_extrap257_left_ccx={extrap257_left}");
        assert_eq!(left_ccx, W * KB, "barrel shift should be W*logW cswaps");
    }

    #[test]
    fn plusminus_signed_barrel_shift_roundtrip_is_phase_clean() {
        // Signed coefficient-lane counterpart to the unsigned denominator
        // barrel shift.  Valid inputs satisfy signed no-overflow for every
        // active stage; shifted-out sign bits are cleaned using the post-shift
        // sign bit.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        const KB: usize = 4;
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(W);
        let kbits = b.alloc_qubits(KB);
        let spill = b.alloc_qubits(W);
        let start_left = b.ops.len();
        emit_barrel_left_shift_signed_nooverflow_for_plusminus(&mut b, &v, &kbits, &spill);
        let left_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_left..]);
        let start_inv = b.ops.len();
        emit_barrel_left_shift_signed_nooverflow_inverse_for_plusminus(&mut b, &v, &kbits, &spill);
        let inv_ccx = local_count_ccx_for_plusminus_cost(&b.ops[start_inv..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1i64 << W) - 1;
        let lo = -(1i64 << (W - 1));
        let hi = (1i64 << (W - 1)) - 1;
        for k in 0u64..8u64 {
            for x in -128i64..128i64 {
                let shifted = x << k;
                if shifted < lo || shifted > hi { continue; }
                let raw = (x & mask) as u64;
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"plusminus-signed-barrel-shift-roundtrip-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                set_slice_u512_pm(&mut sim, &v, U512::from(raw));
                set_slice_u512_pm(&mut sim, &kbits, U512::from(k));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & ((1u64 << W) - 1), raw, "v changed k={k} x={x}");
                assert_eq!(get_slice_u512_pm(&sim, &kbits).as_limbs()[0], k, "k changed");
                assert_eq!(get_slice_u512_pm(&sim, &spill), U512::ZERO, "spill dirty k={k} x={x}");
                assert_eq!(sim.global_phase() & 1, 0, "unexpected phase k={k} x={x}");
            }
        }
        let mut b64 = super::super::B::new();
        let v64 = b64.alloc_qubits(64);
        let k64 = b64.alloc_qubits(6);
        let spill64 = b64.alloc_qubits(64);
        let s64 = b64.ops.len();
        emit_barrel_left_shift_signed_nooverflow_for_plusminus(&mut b64, &v64, &k64, &spill64);
        let w64_left_ccx = local_count_ccx_for_plusminus_cost(&b64.ops[s64..]);
        let extrap257_left = ((w64_left_ccx as f64) * (257.0 / 64.0) * (8.0 / 6.0)).round() as usize;
        eprintln!("plus-minus signed barrel shift: w16 left={left_ccx}, inv={inv_ccx}, peak={peak}, w64_left={w64_left_ccx}, extrap257_left={extrap257_left}");
        println!("METRIC plusminus_signed_barrel_w16_left_ccx={left_ccx}");
        println!("METRIC plusminus_signed_barrel_w16_inverse_ccx={inv_ccx}");
        println!("METRIC plusminus_signed_barrel_w16_peak_q={peak}");
        println!("METRIC plusminus_signed_barrel_w64_left_ccx={w64_left_ccx}");
        println!("METRIC plusminus_signed_barrel_extrap257_left_ccx={extrap257_left}");
        assert!(left_ccx > W * KB && left_ccx < 2 * W * KB, "signed cleanup should be small over W log W");
    }

    #[test]
    fn plusminus_unary_history_to_binary_k_roundtrip_cost() {
        // Barrel shifts need binary k.  This simple converter counts unary
        // prefix bits with controlled +1 additions, then uncomputes the counter
        // from the still-stored unary history after the shifts.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        const KB: usize = 5;
        let mut bc = super::super::B::new();
        let hist = bc.alloc_qubits(W);
        let kbits = bc.alloc_qubits(KB);
        let start_c = bc.ops.len();
        emit_unary_history_to_binary_count_for_plusminus(&mut bc, &hist, &kbits);
        let compute_ccx = local_count_ccx_for_plusminus_cost(&bc.ops[start_c..]);
        let compute_peak = bc.peak_qubits;
        let c_num_qubits = bc.next_qubit as usize;
        let c_num_bits = bc.next_bit as usize;
        let c_ops = bc.ops;
        for k in 0u64..=W as u64 {
            let hist_val = if k == W as u64 { (1u64 << W) - 1 } else { (1u64 << k) - 1 };
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-unary-to-binary-compute-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(c_num_qubits, c_num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &hist, U512::from(hist_val));
            sim.apply(&c_ops);
            assert_eq!(get_slice_u512_pm(&sim, &hist).as_limbs()[0] & ((1u64 << W) - 1), hist_val, "hist changed k={k}");
            assert_eq!(get_slice_u512_pm(&sim, &kbits).as_limbs()[0] & ((1u64 << KB) - 1), k, "kbits mismatch k={k}");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase k={k}");
        }

        let mut br = super::super::B::new();
        let hist_r = br.alloc_qubits(W);
        let kbits_r = br.alloc_qubits(KB);
        let start_r = br.ops.len();
        emit_unary_history_to_binary_count_for_plusminus(&mut br, &hist_r, &kbits_r);
        emit_unary_history_to_binary_count_inverse_for_plusminus(&mut br, &hist_r, &kbits_r);
        let roundtrip_ccx = local_count_ccx_for_plusminus_cost(&br.ops[start_r..]);
        let roundtrip_peak = br.peak_qubits;
        let r_num_qubits = br.next_qubit as usize;
        let r_num_bits = br.next_bit as usize;
        let r_ops = br.ops;
        for k in 0u64..=W as u64 {
            let hist_val = if k == W as u64 { (1u64 << W) - 1 } else { (1u64 << k) - 1 };
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-unary-to-binary-roundtrip-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(r_num_qubits, r_num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &hist_r, U512::from(hist_val));
            sim.apply(&r_ops);
            assert_eq!(get_slice_u512_pm(&sim, &hist_r).as_limbs()[0] & ((1u64 << W) - 1), hist_val, "hist changed roundtrip k={k}");
            assert_eq!(get_slice_u512_pm(&sim, &kbits_r), U512::ZERO, "kbits dirty roundtrip k={k}");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase roundtrip k={k}");
        }

        let mut b64 = super::super::B::new();
        let h64 = b64.alloc_qubits(64);
        let k64 = b64.alloc_qubits(7);
        let s64 = b64.ops.len();
        emit_unary_history_to_binary_count_for_plusminus(&mut b64, &h64, &k64);
        emit_unary_history_to_binary_count_inverse_for_plusminus(&mut b64, &h64, &k64);
        let w64_roundtrip_ccx = local_count_ccx_for_plusminus_cost(&b64.ops[s64..]);
        let extrap257_roundtrip = ((w64_roundtrip_ccx as f64) * (257.0 / 64.0) * (8.0 / 7.0)).round() as usize;
        eprintln!("plus-minus unary-to-binary k: w16 compute={compute_ccx}, compute_peak={compute_peak}, roundtrip={roundtrip_ccx}, roundtrip_peak={roundtrip_peak}, w64_roundtrip={w64_roundtrip_ccx}, extrap257_roundtrip={extrap257_roundtrip}");
        println!("METRIC plusminus_unary_to_binary_w16_compute_ccx={compute_ccx}");
        println!("METRIC plusminus_unary_to_binary_w16_compute_peak_q={compute_peak}");
        println!("METRIC plusminus_unary_to_binary_w16_roundtrip_ccx={roundtrip_ccx}");
        println!("METRIC plusminus_unary_to_binary_w16_roundtrip_peak_q={roundtrip_peak}");
        println!("METRIC plusminus_unary_to_binary_w64_roundtrip_ccx={w64_roundtrip_ccx}");
        println!("METRIC plusminus_unary_to_binary_extrap257_roundtrip_ccx={extrap257_roundtrip}");
        assert!(roundtrip_ccx == 2 * compute_ccx, "roundtrip should be compute+uncompute");
    }

    #[test]
    fn plusminus_barrel_shift_konly_step_roundtrip_and_scaling() {
        // Integrated replacement for the physically repeated shifts inside a
        // k-only step: unary history -> binary k, unsigned/signed barrel shifts,
        // then uncompute binary k from the stored unary history.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubits(W);
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        emit_plusminus_inplace_step_forward_konly_barrel_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, &spill, flag, one);
        emit_plusminus_inplace_step_inverse_konly_barrel_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, &spill, flag, one);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(91u64, 27u64), (201, 77), (255, 127), (987, 31), (65535, 1)];
        for &(uval, vval) in &cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-barrel-konly-step-roundtrip-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u, U512::from(uval));
            set_slice_u512_pm(&mut sim, &v, U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv, U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &hist), U512::ZERO, "hist not clean case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &spill), U512::ZERO, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "flag not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(one) & 1, 0, "one not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }

        let mut b64 = super::super::B::new();
        let u64 = b64.alloc_qubits(64);
        let v64 = b64.alloc_qubits(64);
        let cu64 = b64.alloc_qubits(64);
        let cv64 = b64.alloc_qubits(64);
        let active64 = b64.alloc_qubits(65);
        let hist64 = b64.alloc_qubits(64);
        let spill64 = b64.alloc_qubits(64);
        let flag64 = b64.alloc_qubit();
        let one64 = b64.alloc_qubit();
        let s64 = b64.ops.len();
        emit_plusminus_inplace_step_forward_konly_barrel_for_test(&mut b64, &u64, &v64, &cu64, &cv64, &active64, &hist64, &spill64, flag64, one64);
        let w64_forward_ccx = local_count_ccx_for_plusminus_cost(&b64.ops[s64..]);
        let extrap257_forward = ((w64_forward_ccx as f64) * (257.0 / 64.0) * (8.0 / 6.0)).round() as usize;
        eprintln!("plus-minus barrel k-only step: w16_roundtrip={ccx}, peak={peak}, w64_forward={w64_forward_ccx}, extrap257_forward={extrap257_forward}");
        println!("METRIC plusminus_barrel_step_w16_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_barrel_step_w16_peak_q={peak}");
        println!("METRIC plusminus_barrel_step_w64_forward_ccx={w64_forward_ccx}");
        println!("METRIC plusminus_barrel_step_extrap257_forward_ccx={extrap257_forward}");
        assert!(w64_forward_ccx < 3 * 9338, "barrel step should beat repeated-shift w64 forward by a wide margin");
    }

    #[test]
    fn plusminus_physical_shift_step_scaling_exposes_quadratic_tax() {
        // Cost-only scaling check for the currently wired variable shifts.  The
        // optimistic budget counted shift work per unary-one bit, but this
        // physical skeleton emits one controlled single-bit shift for every
        // possible history bit.  If this grows ~W^2 per step, production needs
        // a relabel/offset representation before 257-bit integration.
        let mut rows = Vec::new();
        for &w in &[8usize, 16, 32, 64] {
            let mut b = super::super::B::new();
            let u = b.alloc_qubits(w);
            let v = b.alloc_qubits(w);
            let cu = b.alloc_qubits(w);
            let cv = b.alloc_qubits(w);
            let active = b.alloc_qubits(w + 1);
            let hist = b.alloc_qubits(w);
            let spill = b.alloc_qubit();
            let flag = b.alloc_qubit();
            let one = b.alloc_qubit();
            let start = b.ops.len();
            emit_plusminus_inplace_step_forward_konly_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag, one);
            let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
            rows.push((w, ccx, b.peak_qubits));
        }
        for (w, ccx, peak) in &rows {
            eprintln!("plus-minus physical k-only forward scaling: w={w}, ccx={ccx}, peak={peak}");
        }
        let w8 = rows[0].1;
        let w16 = rows[1].1;
        let w32 = rows[2].1;
        let w64 = rows[3].1;
        let extrap257 = ((w64 as f64) * (257.0f64 / 64.0).powi(2)).round() as usize;
        println!("METRIC plusminus_physical_shift_w8_forward_ccx={w8}");
        println!("METRIC plusminus_physical_shift_w16_forward_ccx={w16}");
        println!("METRIC plusminus_physical_shift_w32_forward_ccx={w32}");
        println!("METRIC plusminus_physical_shift_w64_forward_ccx={w64}");
        println!("METRIC plusminus_physical_shift_extrap257_forward_ccx={extrap257}");
        assert!(w64 > 3 * w32, "current physical shift unexpectedly stopped looking quadratic");
    }

    #[test]
    fn plusminus_active_aware_step_noop_and_roundtrip_is_clean() {
        // Fixed-bound loop smoke test: u==v must be a no-op with an all-zero
        // history word, while u!=v must match the ordinary k-only step.  The
        // active bit is temporary and is cleared from hist!=0 / restored d.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 12;
        let mut bf = super::super::B::new();
        let fu = bf.alloc_qubits(W);
        let fv = bf.alloc_qubits(W);
        let fcu = bf.alloc_qubits(W);
        let fcv = bf.alloc_qubits(W);
        let factive = bf.alloc_qubits(W + 1);
        let fhist = bf.alloc_qubits(W);
        let fspill = bf.alloc_qubit();
        let fflag = bf.alloc_qubit();
        let fstart = bf.ops.len();
        emit_plusminus_inplace_step_forward_konly_active_for_test(&mut bf, &fu, &fv, &fcu, &fcv, &factive, &fhist, fspill, fflag);
        let f_ccx = local_count_ccx_for_plusminus_cost(&bf.ops[fstart..]);
        let f_peak = bf.peak_qubits;
        let f_num_qubits = bf.next_qubit as usize;
        let f_num_bits = bf.next_bit as usize;
        let f_ops = bf.ops;
        let mask = (1u64 << W) - 1;
        let forward_cases = [(91u64, 27u64, true), (201, 77, true), (45, 45, false), (127, 127, false)];
        for &(uval, vval, is_active) in &forward_cases {
            let (mut eu, mut ev, mut ecu, mut ecv) = (uval, vval, 0u64, 1u64);
            let expected_hist = if is_active {
                let k = (uval - vval).trailing_zeros() as usize;
                plusminus_classical_step_mod_width_for_test(&mut eu, &mut ev, &mut ecu, &mut ecv, W);
                (1u64 << k) - 1
            } else {
                0
            };
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-active-aware-step-forward-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(f_num_qubits, f_num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &fu, U512::from(uval));
            set_slice_u512_pm(&mut sim, &fv, U512::from(vval));
            set_slice_u512_pm(&mut sim, &fcu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &fcv, U512::from(1u64));
            sim.apply(&f_ops);
            assert_eq!(get_slice_u512_pm(&sim, &fu).as_limbs()[0] & mask, eu & mask, "forward u mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fv).as_limbs()[0] & mask, ev & mask, "forward v mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcu).as_limbs()[0] & mask, ecu & mask, "forward cu mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcv).as_limbs()[0] & mask, ecv & mask, "forward cv mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fhist).as_limbs()[0] & mask, expected_hist, "forward hist mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &factive), U512::ZERO, "forward active dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fspill) & 1, 0, "forward spill dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fflag) & 1, 0, "forward flag dirty case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "forward unexpected phase case=({uval},{vval})");
        }

        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        emit_plusminus_inplace_step_forward_konly_active_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag);
        emit_plusminus_inplace_step_inverse_konly_active_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for &(uval, vval, _) in &forward_cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-active-aware-step-roundtrip-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u, U512::from(uval));
            set_slice_u512_pm(&mut sim, &v, U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv, U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &hist), U512::ZERO, "hist not clean case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "flag not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }
        eprintln!("plus-minus active-aware one-step: width={W}, forward_ccx={f_ccx}, forward_peak={f_peak}, roundtrip_ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_active_step_width={W}");
        println!("METRIC plusminus_active_step_forward_ccx={f_ccx}");
        println!("METRIC plusminus_active_step_forward_peak_q={f_peak}");
        println!("METRIC plusminus_active_step_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_active_step_roundtrip_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_fixed_bound_packed_active_loop_roundtrip_is_clean() {
        // Fixed public iteration bound with packed high-bit history slots and
        // terminal no-ops.  Some cases converge before STEPS, so later history
        // words must be all zero and the reverse pass must LIFO-clean both real
        // and inactive steps.
        use sha3::digest::{ExtendableOutput, Update};
        const LIVE: usize = 12;
        const TOTAL: usize = 36;
        const STEPS: usize = 8;
        let cases = [(37u64, 5u64), (91, 27), (201, 77), (255, 127), (187, 45), (233, 17), (171, 65), (1001, 33)];
        let mask = (1u64 << LIVE) - 1;

        let mut bf = super::super::B::new();
        let fu_lane = bf.alloc_qubits(TOTAL);
        let fv_lane = bf.alloc_qubits(TOTAL);
        let fcu_lane = bf.alloc_qubits(TOTAL);
        let fcv_lane = bf.alloc_qubits(TOTAL);
        let factive = bf.alloc_qubits(LIVE + 1);
        let mut fslots = Vec::new();
        for i in LIVE..TOTAL {
            fslots.push(fu_lane[i]);
            fslots.push(fv_lane[i]);
            fslots.push(fcu_lane[i]);
            fslots.push(fcv_lane[i]);
        }
        assert!(fslots.len() >= STEPS * LIVE);
        let fhists: Vec<Vec<super::super::QubitId>> = (0..STEPS)
            .map(|s| fslots[s * LIVE..(s + 1) * LIVE].to_vec())
            .collect();
        let fspill = bf.alloc_qubit();
        let fflag = bf.alloc_qubit();
        let fstart = bf.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_active_for_test(&mut bf, &fu_lane[..LIVE], &fv_lane[..LIVE], &fcu_lane[..LIVE], &fcv_lane[..LIVE], &factive, &fhists[step], fspill, fflag);
        }
        let f_ccx = local_count_ccx_for_plusminus_cost(&bf.ops[fstart..]);
        let f_peak = bf.peak_qubits;
        let f_num_qubits = bf.next_qubit as usize;
        let f_num_bits = bf.next_bit as usize;
        let f_ops = bf.ops;
        for &(uval, vval) in &cases {
            let (mut eu, mut ev, mut ecu, mut ecv) = (uval, vval, 0u64, 1u64);
            let mut expected_hists = [0u64; STEPS];
            for h in expected_hists.iter_mut() {
                if eu != ev {
                    let k = (eu - ev).trailing_zeros() as usize;
                    *h = (1u64 << k) - 1;
                    plusminus_classical_step_mod_width_for_test(&mut eu, &mut ev, &mut ecu, &mut ecv, LIVE);
                }
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-fixed-bound-packed-active-forward-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(f_num_qubits, f_num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &fu_lane[..LIVE], U512::from(uval));
            set_slice_u512_pm(&mut sim, &fv_lane[..LIVE], U512::from(vval));
            set_slice_u512_pm(&mut sim, &fcu_lane[..LIVE], U512::ZERO);
            set_slice_u512_pm(&mut sim, &fcv_lane[..LIVE], U512::from(1u64));
            sim.apply(&f_ops);
            assert_eq!(get_slice_u512_pm(&sim, &fu_lane[..LIVE]).as_limbs()[0] & mask, eu & mask, "forward u mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fv_lane[..LIVE]).as_limbs()[0] & mask, ev & mask, "forward v mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcu_lane[..LIVE]).as_limbs()[0] & mask, ecu & mask, "forward cu mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcv_lane[..LIVE]).as_limbs()[0] & mask, ecv & mask, "forward cv mismatch case=({uval},{vval})");
            for (i, hist) in fhists.iter().enumerate() {
                assert_eq!(get_slice_u512_pm(&sim, hist).as_limbs()[0] & mask, expected_hists[i], "forward packed hist {i} mismatch case=({uval},{vval})");
            }
            assert_eq!(get_slice_u512_pm(&sim, &factive), U512::ZERO, "forward active dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fspill) & 1, 0, "forward spill dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fflag) & 1, 0, "forward flag dirty case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "forward unexpected phase case=({uval},{vval})");
        }

        let mut b = super::super::B::new();
        let u_lane = b.alloc_qubits(TOTAL);
        let v_lane = b.alloc_qubits(TOTAL);
        let cu_lane = b.alloc_qubits(TOTAL);
        let cv_lane = b.alloc_qubits(TOTAL);
        let active = b.alloc_qubits(LIVE + 1);
        let mut slots = Vec::new();
        for i in LIVE..TOTAL {
            slots.push(u_lane[i]);
            slots.push(v_lane[i]);
            slots.push(cu_lane[i]);
            slots.push(cv_lane[i]);
        }
        let hists: Vec<Vec<super::super::QubitId>> = (0..STEPS)
            .map(|s| slots[s * LIVE..(s + 1) * LIVE].to_vec())
            .collect();
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_active_for_test(&mut b, &u_lane[..LIVE], &v_lane[..LIVE], &cu_lane[..LIVE], &cv_lane[..LIVE], &active, &hists[step], spill, flag);
        }
        for step in (0..STEPS).rev() {
            emit_plusminus_inplace_step_inverse_konly_active_for_test(&mut b, &u_lane[..LIVE], &v_lane[..LIVE], &cu_lane[..LIVE], &cv_lane[..LIVE], &active, &hists[step], spill, flag);
        }
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for &(uval, vval) in &cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-fixed-bound-packed-active-roundtrip-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u_lane[..LIVE], U512::from(uval));
            set_slice_u512_pm(&mut sim, &v_lane[..LIVE], U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu_lane[..LIVE], U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv_lane[..LIVE], U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u_lane[..LIVE]).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v_lane[..LIVE]).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu_lane[..LIVE]).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv_lane[..LIVE]).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            for (i, &slot) in slots.iter().enumerate() {
                assert_eq!(sim.qubit(slot) & 1, 0, "packed slot {i} not clean case=({uval},{vval})");
            }
            assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "flag not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }
        let packed_slots = STEPS * LIVE;
        eprintln!("plus-minus fixed-bound packed active loop: live={LIVE}, total={TOTAL}, steps={STEPS}, packed_slots={packed_slots}, forward_ccx={f_ccx}, forward_peak={f_peak}, roundtrip_ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_fixed_bound_live_width={LIVE}");
        println!("METRIC plusminus_fixed_bound_total_width={TOTAL}");
        println!("METRIC plusminus_fixed_bound_steps={STEPS}");
        println!("METRIC plusminus_fixed_bound_packed_slots={packed_slots}");
        println!("METRIC plusminus_fixed_bound_forward_ccx={f_ccx}");
        println!("METRIC plusminus_fixed_bound_forward_peak_q={f_peak}");
        println!("METRIC plusminus_fixed_bound_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_fixed_bound_roundtrip_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_inplace_three_step_konly_roundtrip_is_clean() {
        // Multi-step smoke test with no persistent direction flags.  The flag
        // is cleared after every forward step by coefficient divisibility and
        // recomputed during inverse from the same live state plus k-history.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        const STEPS: usize = 3;
        let cases = [(91u64, 27u64), (201, 77), (255, 127), (987, 31)];
        let mask = (1u64 << W) - 1;

        let mut bf = super::super::B::new();
        let fu = bf.alloc_qubits(W);
        let fv = bf.alloc_qubits(W);
        let fcu = bf.alloc_qubits(W);
        let fcv = bf.alloc_qubits(W);
        let factive = bf.alloc_qubits(W + 1);
        let fhists: Vec<Vec<super::super::QubitId>> = (0..STEPS).map(|_| bf.alloc_qubits(W)).collect();
        let fspill = bf.alloc_qubit();
        let fflag = bf.alloc_qubit();
        let fone = bf.alloc_qubit();
        let fstart = bf.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_for_test(&mut bf, &fu, &fv, &fcu, &fcv, &factive, &fhists[step], fspill, fflag, fone);
        }
        let f_ccx = local_count_ccx_for_plusminus_cost(&bf.ops[fstart..]);
        let f_peak = bf.peak_qubits;
        let f_num_qubits = bf.next_qubit as usize;
        let f_num_bits = bf.next_bit as usize;
        let f_ops = bf.ops;
        for &(uval, vval) in &cases {
            let (mut eu, mut ev, mut ecu, mut ecv) = (uval, vval, 0u64, 1u64);
            for _ in 0..STEPS {
                plusminus_classical_step_mod_width_for_test(&mut eu, &mut ev, &mut ecu, &mut ecv, W);
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-inplace-three-step-konly-forward-v2");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(f_num_qubits, f_num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &fu, U512::from(uval));
            set_slice_u512_pm(&mut sim, &fv, U512::from(vval));
            set_slice_u512_pm(&mut sim, &fcu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &fcv, U512::from(1u64));
            sim.apply(&f_ops);
            assert_eq!(get_slice_u512_pm(&sim, &fu).as_limbs()[0] & mask, eu & mask, "forward u mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fv).as_limbs()[0] & mask, ev & mask, "forward v mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcu).as_limbs()[0] & mask, ecu & mask, "forward cu mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &fcv).as_limbs()[0] & mask, ecv & mask, "forward cv mismatch case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &factive), U512::ZERO, "forward active dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fspill) & 1, 0, "forward spill dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fflag) & 1, 0, "forward recovered flag dirty case=({uval},{vval})");
            assert_eq!(sim.qubit(fone) & 1, 0, "forward one dirty case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "forward unexpected phase case=({uval},{vval})");
        }

        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hists: Vec<Vec<super::super::QubitId>> = (0..STEPS).map(|_| b.alloc_qubits(W)).collect();
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flag, one);
        }
        for step in (0..STEPS).rev() {
            emit_plusminus_inplace_step_inverse_konly_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flag, one);
        }
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for &(uval, vval) in &cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-inplace-three-step-konly-roundtrip-v2");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u, U512::from(uval));
            set_slice_u512_pm(&mut sim, &v, U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv, U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            for (i, hist) in hists.iter().enumerate() {
                assert_eq!(get_slice_u512_pm(&sim, hist), U512::ZERO, "hist {i} not clean case=({uval},{vval})");
            }
            assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "recovered direction flag not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(one) & 1, 0, "one not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }
        eprintln!("plus-minus in-place three-step k-only roundtrip: width={W}, steps={STEPS}, forward_ccx={f_ccx}, forward_peak={f_peak}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_konly_three_step_width={W}");
        println!("METRIC plusminus_konly_three_step_steps={STEPS}");
        println!("METRIC plusminus_konly_three_step_forward_ccx={f_ccx}");
        println!("METRIC plusminus_konly_three_step_forward_peak_q={f_peak}");
        println!("METRIC plusminus_konly_three_step_ccx={ccx}");
        println!("METRIC plusminus_konly_three_step_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_inplace_width32_eight_step_konly_roundtrip_is_clean() {
        // Wider/longer k-only stress smoke test.  This exercises unsigned
        // denominator left-shift restoration near the top of the word, signed
        // coefficient shifts, repeated direction recovery, and LIFO k-history
        // cleanup beyond the single-step happy path.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 32;
        const STEPS: usize = 8;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hists: Vec<Vec<super::super::QubitId>> = (0..STEPS).map(|_| b.alloc_qubits(W)).collect();
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flag, one);
        }
        for step in (0..STEPS).rev() {
            emit_plusminus_inplace_step_inverse_konly_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flag, one);
        }
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [
            (123456789u64, 98765u64),
            (400000001, 1234567),
            (0xfffffff1, 12345),
            (65537, 17),
            (987654321, 123456789),
            (2147483647, 1),
        ];
        for &(uval, vval) in &cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-inplace-width32-eight-step-konly-v2");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u, U512::from(uval));
            set_slice_u512_pm(&mut sim, &v, U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu, U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv, U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            for (i, hist) in hists.iter().enumerate() {
                assert_eq!(get_slice_u512_pm(&sim, hist), U512::ZERO, "hist {i} not clean case=({uval},{vval})");
            }
            assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "direction recovery flag not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(one) & 1, 0, "one not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }
        eprintln!("plus-minus width32 eight-step k-only roundtrip: width={W}, steps={STEPS}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_konly_w32_steps8_width={W}");
        println!("METRIC plusminus_konly_w32_steps8_steps={STEPS}");
        println!("METRIC plusminus_konly_w32_steps8_ccx={ccx}");
        println!("METRIC plusminus_konly_w32_steps8_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_slack_slot_history_roundtrip_is_clean() {
        // First toy slack-packing circuit: the per-step k-history qubits are
        // not separate registers.  They are fixed public high-bit slots inside
        // the four live lanes; arithmetic touches only the low LIVE bits.  This
        // validates the core reversible layout idea before implementing the
        // full 256-bit public envelope.
        use sha3::digest::{ExtendableOutput, Update};
        const LIVE: usize = 12;
        const TOTAL: usize = 24;
        const STEPS: usize = 3;
        let mut b = super::super::B::new();
        let u_lane = b.alloc_qubits(TOTAL);
        let v_lane = b.alloc_qubits(TOTAL);
        let cu_lane = b.alloc_qubits(TOTAL);
        let cv_lane = b.alloc_qubits(TOTAL);
        let active = b.alloc_qubits(LIVE + 1);
        let mut slots = Vec::new();
        for i in LIVE..TOTAL {
            slots.push(u_lane[i]);
            slots.push(v_lane[i]);
            slots.push(cu_lane[i]);
            slots.push(cv_lane[i]);
        }
        assert!(slots.len() >= STEPS * LIVE);
        let hists: Vec<Vec<super::super::QubitId>> = (0..STEPS)
            .map(|s| slots[s * LIVE..(s + 1) * LIVE].to_vec())
            .collect();
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_konly_for_test(&mut b, &u_lane[..LIVE], &v_lane[..LIVE], &cu_lane[..LIVE], &cv_lane[..LIVE], &active, &hists[step], spill, flag, one);
        }
        for step in (0..STEPS).rev() {
            emit_plusminus_inplace_step_inverse_konly_for_test(&mut b, &u_lane[..LIVE], &v_lane[..LIVE], &cu_lane[..LIVE], &cv_lane[..LIVE], &active, &hists[step], spill, flag, one);
        }
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << LIVE) - 1;
        let cases = [(91u64, 27u64), (201, 77), (255, 127), (187, 45), (233, 17), (171, 65)];
        for &(uval, vval) in &cases {
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"plusminus-slack-slot-history-roundtrip-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_pm(&mut sim, &u_lane[..LIVE], U512::from(uval));
            set_slice_u512_pm(&mut sim, &v_lane[..LIVE], U512::from(vval));
            set_slice_u512_pm(&mut sim, &cu_lane[..LIVE], U512::ZERO);
            set_slice_u512_pm(&mut sim, &cv_lane[..LIVE], U512::from(1u64));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_pm(&sim, &u_lane[..LIVE]).as_limbs()[0] & mask, uval & mask, "u changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &v_lane[..LIVE]).as_limbs()[0] & mask, vval & mask, "v changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cu_lane[..LIVE]).as_limbs()[0] & mask, 0, "cu changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &cv_lane[..LIVE]).as_limbs()[0] & mask, 1, "cv changed case=({uval},{vval})");
            assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean case=({uval},{vval})");
            for (i, &slot) in slots.iter().enumerate() {
                assert_eq!(sim.qubit(slot) & 1, 0, "packed slot {i} not clean case=({uval},{vval})");
            }
            assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(flag) & 1, 0, "direction flag not clean case=({uval},{vval})");
            assert_eq!(sim.qubit(one) & 1, 0, "one not clean case=({uval},{vval})");
            assert_eq!(sim.global_phase() & 1, 0, "unexpected phase case=({uval},{vval})");
        }
        let packed_slots = STEPS * LIVE;
        eprintln!("plus-minus slack-slot k-history roundtrip: live={LIVE}, total={TOTAL}, steps={STEPS}, packed_slots={packed_slots}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_slack_slot_live_width={LIVE}");
        println!("METRIC plusminus_slack_slot_total_width={TOTAL}");
        println!("METRIC plusminus_slack_slot_steps={STEPS}");
        println!("METRIC plusminus_slack_slot_packed_slots={packed_slots}");
        println!("METRIC plusminus_slack_slot_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_slack_slot_roundtrip_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_inplace_one_step_roundtrip_is_clean() {
        // Actual in-place forward followed by explicit inverse. This is the
        // critical quantum-compatibility test for a wireable step: history and
        // direction are produced, used to reverse, and then cleaned with no
        // fresh old/new lane coexistence.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        emit_plusminus_inplace_step_forward_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag, one);
        emit_plusminus_inplace_step_inverse_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag, one);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(37u64, 5u64), (91, 27), (128, 64), (201, 77), (255, 127)];
        for &(uval, vval) in &cases {
            for cuv in [0u64, 1, 7, 123] {
                for cvv in [0u64, 3, 11, 15] {
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-inplace-one-step-roundtrip-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &u, U512::from(uval));
                    set_slice_u512_pm(&mut sim, &v, U512::from(vval));
                    set_slice_u512_pm(&mut sim, &cu, U512::from(cuv));
                    set_slice_u512_pm(&mut sim, &cv, U512::from(cvv));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed");
                    assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, cuv & mask, "cu changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, cvv & mask, "cv changed");
                    assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean");
                    assert_eq!(get_slice_u512_pm(&sim, &hist), U512::ZERO, "hist not clean");
                    assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                    assert_eq!(sim.qubit(flag) & 1, 0, "flag not clean");
                    assert_eq!(sim.qubit(one) & 1, 0, "one not clean");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
                }
            }
        }
        eprintln!("plus-minus in-place one-step roundtrip: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_inplace_roundtrip_width={W}");
        println!("METRIC plusminus_inplace_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_inplace_roundtrip_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_inplace_three_step_roundtrip_is_clean() {
        // Multi-step smoke test with separate public history slots. This catches
        // whether the in-place step can actually be chained like a DIV loop,
        // not just run once in isolation.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        const STEPS: usize = 3;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hists: Vec<Vec<super::super::QubitId>> = (0..STEPS).map(|_| b.alloc_qubits(W)).collect();
        let flags = b.alloc_qubits(STEPS);
        let spill = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        for step in 0..STEPS {
            emit_plusminus_inplace_step_forward_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flags[step], one);
        }
        for step in (0..STEPS).rev() {
            emit_plusminus_inplace_step_inverse_for_test(&mut b, &u, &v, &cu, &cv, &active, &hists[step], spill, flags[step], one);
        }
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(91u64, 27u64), (128, 64), (201, 77), (255, 127), (1000, 17), (987, 31)];
        for &(uval, vval) in &cases {
            for cuv in [0u64, 1, 7, 123] {
                for cvv in [0u64, 1, 3, 5] {
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-inplace-three-step-roundtrip-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &u, U512::from(uval));
                    set_slice_u512_pm(&mut sim, &v, U512::from(vval));
                    set_slice_u512_pm(&mut sim, &cu, U512::from(cuv));
                    set_slice_u512_pm(&mut sim, &cv, U512::from(cvv));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed");
                    assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, cuv & mask, "cu changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, cvv & mask, "cv changed");
                    assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean");
                    for (i, hist) in hists.iter().enumerate() {
                        assert_eq!(get_slice_u512_pm(&sim, hist), U512::ZERO, "hist {i} not clean");
                    }
                    for (i, &flag) in flags.iter().enumerate() {
                        assert_eq!(sim.qubit(flag) & 1, 0, "flag {i} not clean");
                    }
                    assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                    assert_eq!(sim.qubit(one) & 1, 0, "one not clean");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
                }
            }
        }
        eprintln!("plus-minus in-place three-step roundtrip: width={W}, steps={STEPS}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_inplace_three_step_width={W}");
        println!("METRIC plusminus_inplace_three_step_steps={STEPS}");
        println!("METRIC plusminus_inplace_three_step_ccx={ccx}");
        println!("METRIC plusminus_inplace_three_step_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_inplace_one_step_forward_matches_classical() {
        // First genuinely low-scratch productive step: mutate (u,v,cu,cv) in
        // place and leave unary k plus one direction bit as history. This is
        // the shape that can actually be wired; no fresh output lanes coexist
        // with old state.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        emit_plusminus_inplace_step_forward_for_test(&mut b, &u, &v, &cu, &cv, &active, &hist, spill, flag, one);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(37u64, 5u64), (91, 27), (128, 64), (201, 77), (255, 127)];
        for &(uval, vval) in &cases {
            for cuv in [0u64, 1, 7, 123] {
                for cvv in [0u64, 3, 11, 15] {
                    let diff = uval - vval;
                    let k = diff.trailing_zeros() as usize;
                    let d_class = diff >> k;
                    let cd_class = cuv.wrapping_sub(cvv) & mask;
                    let cvs_class = (cvv << k) & mask;
                    let dir = d_class < vval;
                    let (enu, env, encu, encv) = if dir {
                        (vval, d_class, cvs_class, cd_class)
                    } else {
                        (d_class, vval, cd_class, cvs_class)
                    };
                    let expected_hist = if k >= W { (1u64 << W) - 1 } else { (1u64 << k) - 1 };
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-inplace-one-step-forward-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &u, U512::from(uval));
                    set_slice_u512_pm(&mut sim, &v, U512::from(vval));
                    set_slice_u512_pm(&mut sim, &cu, U512::from(cuv));
                    set_slice_u512_pm(&mut sim, &cv, U512::from(cvv));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, enu & mask, "u mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, env & mask, "v mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, encu & mask, "cu mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, encv & mask, "cv mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &hist).as_limbs()[0] & mask, expected_hist, "hist mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean");
                    assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                    assert_eq!((sim.qubit(flag) & 1) != 0, dir, "direction flag mismatch");
                    assert_eq!(sim.qubit(one) & 1, 0, "one not clean");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
                }
            }
        }
        eprintln!("plus-minus in-place one-step forward: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_inplace_step_width={W}");
        println!("METRIC plusminus_inplace_step_ccx={ccx}");
        println!("METRIC plusminus_inplace_step_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_one_step_output_skeleton_matches_classical() {
        // Productive one-step skeleton: compute the ordered next denominator
        // and scaled coefficient lanes into fresh output registers, then clean
        // all temporary work.  This is still not the final in-place low-scratch
        // update, but it validates that the real primitives compose into a
        // phase-clean, non-noop plus-minus step.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let d = b.alloc_qubits(W);
        let dsh = b.alloc_qubits(W);
        let cd = b.alloc_qubits(W);
        let cvs = b.alloc_qubits(W);
        let nu = b.alloc_qubits(W);
        let nv = b.alloc_qubits(W);
        let ncu = b.alloc_qubits(W);
        let ncv = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let one = b.alloc_qubit();
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        b.x(one);

        for i in 0..W { b.cx(u[i], d[i]); }
        super::super::sub_nbit_qq_fast(&mut b, &v, &d); // d = u-v
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);

        for i in 0..W { b.cx(d[i], dsh[i]); }
        for &h in &hist {
            emit_controlled_right_shift_exact_for_plusminus(&mut b, &dsh, h, spill);
        }
        for i in 0..W { b.cx(cv[i], cvs[i]); }
        for &h in &hist {
            emit_controlled_left_shift_nooverflow_for_plusminus(&mut b, &cvs, h, spill);
        }
        for i in 0..W { b.cx(cu[i], cd[i]); }
        emit_controlled_integer_add_for_plusminus(&mut b, &cd, &cv, one, true); // cd = cu-cv

        for i in 0..W {
            b.cx(dsh[i], nu[i]);
            b.cx(v[i], nv[i]);
            b.cx(cd[i], ncu[i]);
            b.cx(cvs[i], ncv[i]);
        }
        super::super::with_lt(&mut b, &dsh, &v, flag, |b| {
            for i in 0..W {
                local_cswap_for_plusminus_cost(b, flag, nu[i], nv[i]);
                local_cswap_for_plusminus_cost(b, flag, ncu[i], ncv[i]);
            }
        });

        emit_controlled_integer_add_for_plusminus(&mut b, &cd, &cv, one, false);
        for i in (0..W).rev() { b.cx(cu[i], cd[i]); }
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_nooverflow_inverse_for_plusminus(&mut b, &cvs, h, spill);
        }
        for i in (0..W).rev() { b.cx(cv[i], cvs[i]); }
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_unsigned_exact_for_plusminus(&mut b, &dsh, h, spill);
        }
        for i in (0..W).rev() { b.cx(d[i], dsh[i]); }
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        super::super::add_nbit_qq_fast(&mut b, &v, &d);
        for i in (0..W).rev() { b.cx(u[i], d[i]); }
        b.x(one);

        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(37u64, 5u64), (91, 27), (128, 64), (201, 77), (255, 127)];
        for &(uval, vval) in &cases {
            for cuv in [0u64, 1, 7, 123] {
                for cvv in [0u64, 3, 11, 15] {
                    let diff = uval - vval;
                    let k = diff.trailing_zeros() as usize;
                    let d_class = diff >> k;
                    let cd_class = cuv.wrapping_sub(cvv) & mask;
                    let cvs_class = (cvv << k) & mask;
                    let (enu, env, encu, encv) = if d_class < vval {
                        (vval, d_class, cvs_class, cd_class)
                    } else {
                        (d_class, vval, cd_class, cvs_class)
                    };
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-one-step-output-skeleton-v2");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &u, U512::from(uval));
                    set_slice_u512_pm(&mut sim, &v, U512::from(vval));
                    set_slice_u512_pm(&mut sim, &cu, U512::from(cuv));
                    set_slice_u512_pm(&mut sim, &cv, U512::from(cvv));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_pm(&sim, &nu).as_limbs()[0] & mask, enu & mask, "nu mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &nv).as_limbs()[0] & mask, env & mask, "nv mismatch");
                    assert_eq!(get_slice_u512_pm(&sim, &ncu).as_limbs()[0] & mask, encu & mask, "ncu mismatch u={uval} v={vval} cu={cuv} cv={cvv} k={k}");
                    assert_eq!(get_slice_u512_pm(&sim, &ncv).as_limbs()[0] & mask, encv & mask, "ncv mismatch");
                    for (name, reg) in [("d", &d), ("dsh", &dsh), ("cd", &cd), ("cvs", &cvs), ("active", &active), ("hist", &hist)] {
                        assert_eq!(get_slice_u512_pm(&sim, reg), U512::ZERO, "{name} not clean");
                    }
                    assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                    assert_eq!(sim.qubit(one) & 1, 0, "one not clean");
                    assert_eq!(sim.qubit(flag) & 1, 0, "flag not clean");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
                }
            }
        }
        eprintln!("plus-minus one-step output skeleton: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_step_output_width={W}");
        println!("METRIC plusminus_step_output_ccx={ccx}");
        println!("METRIC plusminus_step_output_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_one_step_bennett_skeleton_is_phase_clean() {
        // Integration-risk test: compose the actual pieces used in one scaled
        // plus-minus step in a Bennett compute/uncompute shell.  It computes
        // d=u-v, derives active-chain k history from d, computes cv<<k and
        // cd=cu-cv, then reverses everything.  This is not the production
        // in-place update, but it exercises the real MBU add/sub, active-chain
        // controls, controlled shifts, and inverse shifts together.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let u = b.alloc_qubits(W);
        let v = b.alloc_qubits(W);
        let cu = b.alloc_qubits(W);
        let cv = b.alloc_qubits(W);
        let d = b.alloc_qubits(W);
        let cd = b.alloc_qubits(W);
        let cvs = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let one = b.alloc_qubit();
        let start = b.ops.len();
        b.x(one);
        for i in 0..W { b.cx(u[i], d[i]); }
        super::super::sub_nbit_qq_fast(&mut b, &v, &d);
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        for i in 0..W { b.cx(cv[i], cvs[i]); }
        for &h in &hist {
            emit_controlled_left_shift_nooverflow_for_plusminus(&mut b, &cvs, h, spill);
        }
        for i in 0..W { b.cx(cu[i], cd[i]); }
        emit_controlled_integer_add_for_plusminus(&mut b, &cd, &cv, one, true);
        emit_controlled_integer_add_for_plusminus(&mut b, &cd, &cv, one, false);
        for i in (0..W).rev() { b.cx(cu[i], cd[i]); }
        for &h in hist.iter().rev() {
            emit_controlled_left_shift_nooverflow_inverse_for_plusminus(&mut b, &cvs, h, spill);
        }
        for i in (0..W).rev() { b.cx(cv[i], cvs[i]); }
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        super::super::add_nbit_qq_fast(&mut b, &v, &d);
        for i in (0..W).rev() { b.cx(u[i], d[i]); }
        b.x(one);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1u64 << W) - 1;
        let cases = [(37u64, 5u64), (91, 27), (128, 64), (255, 127), (409, 137)];
        for &(uval, vval) in &cases {
            for cuv in [0u64, 1, 7, 123, 0x7fffu64] {
                for cvv in [0u64, 3, 11, 77, 0x003fu64] {
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"plusminus-one-step-bennett-skeleton-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_pm(&mut sim, &u, U512::from(uval));
                    set_slice_u512_pm(&mut sim, &v, U512::from(vval));
                    set_slice_u512_pm(&mut sim, &cu, U512::from(cuv));
                    set_slice_u512_pm(&mut sim, &cv, U512::from(cvv));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_pm(&sim, &u).as_limbs()[0] & mask, uval & mask, "u changed");
                    assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & mask, vval & mask, "v changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cu).as_limbs()[0] & mask, cuv & mask, "cu changed");
                    assert_eq!(get_slice_u512_pm(&sim, &cv).as_limbs()[0] & mask, cvv & mask, "cv changed");
                    for (name, reg) in [("d", &d), ("cd", &cd), ("cvs", &cvs), ("active", &active), ("hist", &hist)] {
                        assert_eq!(get_slice_u512_pm(&sim, reg), U512::ZERO, "{name} not clean");
                    }
                    assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                    assert_eq!(sim.qubit(one) & 1, 0, "one not clean");
                    assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
                }
            }
        }
        eprintln!("plus-minus one-step Bennett skeleton: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_step_skeleton_width={W}");
        println!("METRIC plusminus_step_skeleton_ccx={ccx}");
        println!("METRIC plusminus_step_skeleton_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_active_chain_controls_shift_roundtrip_circuit_is_clean() {
        // Compose the two new primitives in the way the parser needs them:
        // generate unary active-chain controls from d, use them to apply a
        // variable left shift to a signed lane, then rerun the generator to
        // clear the history.  This is a small self-contained reversible circuit
        // and catches exactly the kind of phase/ancilla issue that pure cost
        // models miss.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let d = b.alloc_qubits(W);
        let lane = b.alloc_qubits(W);
        let active = b.alloc_qubits(W + 1);
        let hist = b.alloc_qubits(W);
        let spill = b.alloc_qubit();
        let start = b.ops.len();
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        for &h in &hist {
            emit_controlled_left_shift_nooverflow_for_plusminus(&mut b, &lane, h, spill);
        }
        emit_trailing_zero_active_chain_history_for_plusminus(&mut b, &d, &active, &hist);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1i64 << W) - 1;
        for dval in 1u64..128u64 {
            let k = dval.trailing_zeros() as usize;
            for x in -64i64..64i64 {
                let raw = (x & mask) as u64;
                let expected = ((x << k) & mask) as u64;
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"plusminus-active-chain-shift-roundtrip-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                set_slice_u512_pm(&mut sim, &d, U512::from(dval));
                set_slice_u512_pm(&mut sim, &lane, U512::from(raw));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_pm(&sim, &lane).as_limbs()[0] & ((1u64 << W) - 1), expected, "lane mismatch d={dval} k={k} x={x}");
                assert_eq!(get_slice_u512_pm(&sim, &d).as_limbs()[0] & ((1u64 << W) - 1), dval, "d changed");
                assert_eq!(get_slice_u512_pm(&sim, &active), U512::ZERO, "active not clean");
                assert_eq!(get_slice_u512_pm(&sim, &hist), U512::ZERO, "history not clean");
                assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean");
                assert_eq!(sim.global_phase() & 1, 0, "unexpected phase");
            }
        }
        eprintln!("plus-minus active-chain controlled shift roundtrip: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_active_shift_roundtrip_width={W}");
        println!("METRIC plusminus_active_shift_roundtrip_ccx={ccx}");
        println!("METRIC plusminus_active_shift_roundtrip_peak_q={peak}");
        assert!(ccx > 0 && peak > 0);
    }

    #[test]
    fn plusminus_nooverflow_controlled_shift_circuit_is_clean() {
        // First actual reversible circuit skeleton for the scaled-integer route.
        // It validates the promised no-overflow controlled left shift used by
        // the cost model: valid signed inputs produce 2*x, ctrl=0 is identity,
        // and the spill ancilla returns to zero with no measurement phase.
        use sha3::digest::{ExtendableOutput, Update};
        const W: usize = 16;
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(W);
        let ctrl = b.alloc_qubit();
        let spill = b.alloc_qubit();
        let start = b.ops.len();
        emit_controlled_left_shift_nooverflow_for_plusminus(&mut b, &v, ctrl, spill);
        let ccx = local_count_ccx_for_plusminus_cost(&b.ops[start..]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mask = (1i64 << W) - 1;
        for &ctrl_val in &[false, true] {
            for x in -128i64..128i64 {
                let raw = (x & mask) as u64;
                let expected = if ctrl_val { ((2 * x) & mask) as u64 } else { raw };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"plusminus-nooverflow-shift-circuit-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                set_slice_u512_pm(&mut sim, &v, U512::from(raw));
                if ctrl_val { *sim.qubit_mut(ctrl) |= 1; }
                sim.apply(&ops);
                assert_eq!(get_slice_u512_pm(&sim, &v).as_limbs()[0] & ((1u64 << W) - 1), expected, "shift value mismatch ctrl={ctrl_val} x={x}");
                assert_eq!(sim.qubit(spill) & 1, 0, "spill not clean ctrl={ctrl_val} x={x}");
                assert_eq!(sim.global_phase() & 1, 0, "unexpected phase ctrl={ctrl_val} x={x}");
            }
        }
        eprintln!("plus-minus no-overflow controlled shift circuit: width={W}, ccx={ccx}, peak={peak}");
        println!("METRIC plusminus_shift_circuit_width={W}");
        println!("METRIC plusminus_shift_circuit_ccx={ccx}");
        println!("METRIC plusminus_shift_circuit_peak_q={peak}");
        assert_eq!(ccx, W + 1, "controlled no-overflow shift should cost width+1 CCX");
    }

    #[test]
    fn plusminus_scaled_integer_coefficients_recover_inverse() {
        // Algebra check for the cheap-shift representation: final coefficient c
        // represents 1 = c*x/2^S, so c*2^-S is x^-1 mod p.  This must hold
        // before any circuit cost model using scaled integer coefficients is
        // meaningful.
        let p = SECP256K1_P;
        let samples = 2048usize;
        let mut rng = 0x1eaf_6635_1eed_c0deu64;
        let mut max_scale = 0usize;
        let mut max_width = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let (width, scale, _steps, _initial_twos, coeff) = plusminus_scaled_coeff_width_for_divisor(x, p);
            let c = smag_mod_u256_for_plusminus_test(coeff, p);
            let inv_scale = two_inv_pow_u256_for_plusminus_test(p, scale);
            let inv = c.mul_mod(inv_scale, p);
            assert_eq!(x.mul_mod(inv, p), U256::from(1u64), "scaled plus-minus coefficient did not recover inverse");
            max_scale = max_scale.max(scale);
            max_width = max_width.max(width);
        }
        eprintln!("plus-minus scaled coefficient inverse recovery: samples={samples}, max_scale={max_scale}, max_width={max_width}");
        println!("METRIC plusminus_scaled_coeff_inverse_samples={samples}");
        println!("METRIC plusminus_scaled_coeff_inverse_max_scale={max_scale}");
        println!("METRIC plusminus_scaled_coeff_inverse_max_width={max_width}");
    }

    #[test]
    fn plusminus_scaled_integer_controlled_step_cost_is_sota_shaped_if_overflow_clean() {
        // Replace the modular controlled halve/double tax (1280 CCX) with the
        // operation suggested by the scaled-integer algebra: controlled signed
        // integer add/sub for cd=cu-cv and controlled left-shift/relabel for
        // retained coefficients.  This charges the obvious cswap/add floor but
        // still leaves the hard proof obligations explicit: signed overflow/top
        // bit recovery, direction controls, and full reversible cleanup.
        const WIDTH: usize = 257; // 256 magnitude bits plus sign/guard.
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(WIDTH);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(WIDTH);
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0xadd5_1f75_6635_c0deu64;
        let mut one_div = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            // Two channels (denominator-like + coefficient-like).  Direction
            // muxing and sign/range canonicalization are not included.
            one_div.push(2 * (steps * cint_add_ccx + unary * cshift_ccx));
        }
        one_div.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let scaffold_after_div = 642_716usize;
        let projected_p99 = scaffold_after_div + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus scaled integer controlled step cost: cint_add={cint_add_ccx}, cshift={cshift_ccx}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_scaled_integer_cint_add_ccx={cint_add_ccx}");
        println!("METRIC plusminus_scaled_integer_cshift_ccx={cshift_ccx}");
        println!("METRIC plusminus_scaled_integer_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_scaled_integer_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_scaled_integer_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_scaled_integer_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_scaled_integer_gap_p99_to_2700k={gap_p99}");
        assert!(gap_p99 < 0, "scaled-integer plus-minus floor is not SOTA-shaped before cleanup tax");
    }

    #[test]
    fn plusminus_trailing_zero_unary_generator_tax_still_fits() {
        // Charge a concrete prefix-zero unary k generator for each plus-minus
        // step.  This is the cost of discovering k from d=u-v, separate from
        // the scaled coefficient add/shift floor.  Even a forward+cleanup prefix
        // scan should be much smaller than the modular-control tax if the route
        // is viable.
        const WIDTH: usize = 256;
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(WIDTH);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x7a11_1e20_6635_c0deu64;
        let mut one_div = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            one_div.push(2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * gen_ccx);
        }
        one_div.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let projected_p99 = 642_716usize + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus unary generator tax: gen_ccx={gen_ccx}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_unary_generator_ccx={gen_ccx}");
        println!("METRIC plusminus_unary_generator_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_unary_generator_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_unary_generator_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_unary_generator_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_unary_generator_gap_p99_to_2700k={gap_p99}");
        assert!(gap_p99 < 0, "plus-minus unary generator tax erases SOTA margin");
    }

    #[test]
    fn plusminus_ordering_compare_and_swap_tax_still_fits() {
        // Add the next omitted per-step costs: compare d vs old v to choose the
        // ordered output, then conditionally swap the denominator lane and the
        // paired scaled-coefficient lane.  This still omits sign canonicalization
        // and actual slack-pack moves, but it is the main non-arithmetic control
        // tax in the plus-minus step.
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x0ede_6635_c5aa_900du64;
        let mut one_div = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
            one_div.push(2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * step_tax);
        }
        one_div.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let projected_p99 = 642_716usize + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus ordering tax: cmp={cmp_ccx}, cswap={cswap_ccx}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_order_cmp_ccx={cmp_ccx}");
        println!("METRIC plusminus_order_cswap_ccx={cswap_ccx}");
        println!("METRIC plusminus_order_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_order_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_order_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_order_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_order_gap_p99_to_2700k={gap_p99}");
        assert!(gap_p99 < 0, "plus-minus compare/swap tax erases SOTA margin");
    }

    #[test]
    fn plusminus_posthoc_variable_scale_correction_kills_margin() {
        // The scaled coefficient gives c*x/2^S=1.  If we convert c to the
        // canonical quotient only after DIV, a data-dependent 2^-S correction is
        // needed.  Charging one controlled modular halve per unary shift with
        // existing primitives is the obvious post-hoc route; if it kills the
        // margin, scale must be absorbed into the affine algebra or into a
        // product-clean inverse, not corrected at the end.
        let scale_halve_ccx = {
            let mut b = super::super::B::new();
            let v = b.alloc_qubits(256);
            let ctrl = b.alloc_qubit();
            let start = b.ops.len();
            super::super::cmod_halve_inplace(&mut b, &v, SECP256K1_P, ctrl);
            local_count_ccx_for_plusminus_cost(&b.ops[start..])
        };
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x5ca1_e663_c0de_2026u64;
        let mut one_div = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
            let core = 2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * step_tax;
            one_div.push(core + unary * scale_halve_ccx);
        }
        one_div.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let projected_p99 = 642_716usize + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus posthoc scale correction: scale_halve={scale_halve_ccx}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_scale_halve_ccx={scale_halve_ccx}");
        println!("METRIC plusminus_scale_posthoc_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_scale_posthoc_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_scale_posthoc_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_scale_posthoc_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_scale_posthoc_gap_p99_to_2700k={gap_p99}");
        assert!(gap_p99 > 0, "post-hoc variable scale correction unexpectedly still fits SOTA");
    }

    fn solinas_history_carry_scale_dp_for_plusminus(max_len: usize) -> (Vec<usize>, Vec<usize>) {
        let mut cost_by_k = vec![0usize; 23];
        for k in 1..=22 {
            let (_cur, no_threshold, _exact, _threshold) =
                super::super::primitive_costs::direct_solinas_multihalve_chunk_cost_split(k);
            cost_by_k[k] = no_threshold;
        }
        let inf = usize::MAX / 4;
        let mut dp = vec![inf; max_len + 1];
        let mut chunks = vec![0usize; max_len + 1];
        dp[0] = 0;
        for i in 1..=max_len {
            for k in 1..=22.min(i) {
                let cand = dp[i - k].saturating_add(cost_by_k[k]);
                if cand < dp[i] {
                    dp[i] = cand;
                    chunks[i] = chunks[i - k] + 1;
                }
            }
        }
        (dp, chunks)
    }

    #[test]
    fn plusminus_scale_correction_with_solinas_history_carry_fits_margin() {
        // Replace post-hoc controlled halves by the direct Solinas multihalve
        // history-carry primitive costed earlier.  This is still a model: it
        // assumes the one-bit residual per chunk can be packed with the unary
        // history and that variable-S chunk selection is driven by the same
        // public/slack-packed history.  It decides whether the scale problem is
        // worth a real circuit skeleton.
        let p = SECP256K1_P;
        let samples = 4096usize;
        let mut rng = 0x5011_6635_ca11_2026u64;
        let mut traces = Vec::with_capacity(samples);
        let mut max_scale = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            max_scale = max_scale.max(unary);
            traces.push((unary, steps));
        }
        let (scale_dp, chunk_dp) = solinas_history_carry_scale_dp_for_plusminus(max_scale);
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let mut one_div = Vec::with_capacity(samples);
        let mut scale_costs = Vec::with_capacity(samples);
        let mut scale_chunks = Vec::with_capacity(samples);
        for (unary, steps) in traces {
            let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
            let core = 2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * step_tax;
            one_div.push(core + scale_dp[unary]);
            scale_costs.push(scale_dp[unary]);
            scale_chunks.push(chunk_dp[unary]);
        }
        one_div.sort_unstable();
        scale_costs.sort_unstable();
        scale_chunks.sort_unstable();
        let p99 = samples * 99 / 100;
        let one_div_p99 = one_div[p99];
        let one_div_max = *one_div.last().unwrap();
        let scale_p99 = scale_costs[p99];
        let scale_max = *scale_costs.last().unwrap();
        let chunks_max = *scale_chunks.last().unwrap();
        let two_div_p99 = 2 * one_div_p99;
        let projected_p99 = 642_716usize + two_div_p99;
        let gap_p99 = projected_p99 as isize - 2_700_000isize;
        eprintln!(
            "plus-minus Solinas history-carry scale: max_scale={max_scale}, scale_p99={scale_p99}, scale_max={scale_max}, chunks_max={chunks_max}, one_div_p99={one_div_p99}, projected_p99={projected_p99}, gap_p99={gap_p99}"
        );
        println!("METRIC plusminus_solinas_scale_max_bits={max_scale}");
        println!("METRIC plusminus_solinas_scale_cost_p99_ccx={scale_p99}");
        println!("METRIC plusminus_solinas_scale_cost_max_ccx={scale_max}");
        println!("METRIC plusminus_solinas_scale_chunks_max={chunks_max}");
        println!("METRIC plusminus_solinas_scale_one_div_p99_ccx={one_div_p99}");
        println!("METRIC plusminus_solinas_scale_one_div_max_ccx={one_div_max}");
        println!("METRIC plusminus_solinas_scale_two_div_p99_ccx={two_div_p99}");
        println!("METRIC plusminus_solinas_scale_projected_p99_toffoli={projected_p99}");
        println!("METRIC plusminus_solinas_scale_gap_p99_to_2700k={gap_p99}");
        assert!(gap_p99 < 0, "Solinas history-carry scale correction does not preserve plus-minus margin");
    }

    #[test]
    fn plusminus_active_chain_solinas_stress_max() {
        // Same all-in stress as the Solinas model, but with the active prefix
        // chain serving as unary history (512 CCX/step generator).
        let p = SECP256K1_P;
        let samples = 32_768usize;
        let mut rng = 0xac71_6635_2026_0002u64;
        let mut traces = Vec::with_capacity(samples);
        let mut max_scale = 0usize;
        let mut max_steps = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            max_scale = max_scale.max(unary);
            max_steps = max_steps.max(steps);
            traces.push((unary, steps));
        }
        let (scale_dp, chunk_dp) = solinas_history_carry_scale_dp_for_plusminus(max_scale);
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_active_chain_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let mut projected = Vec::with_capacity(samples);
        let mut chunks = Vec::with_capacity(samples);
        for (unary, steps) in traces {
            let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
            let one_div = 2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * step_tax + scale_dp[unary];
            projected.push(642_716usize + 2 * one_div);
            chunks.push(chunk_dp[unary]);
        }
        projected.sort_unstable();
        chunks.sort_unstable();
        let p99 = samples * 99 / 100;
        let p999 = samples * 999 / 1000;
        let projected_p99 = projected[p99];
        let projected_p999 = projected[p999];
        let projected_max = *projected.last().unwrap();
        let gap_max = projected_max as isize - 2_700_000isize;
        let chunks_max = *chunks.last().unwrap();
        eprintln!(
            "plus-minus active-chain Solinas stress: samples={samples}, max_scale={max_scale}, max_steps={max_steps}, projected_p99={projected_p99}, projected_p999={projected_p999}, projected_max={projected_max}, gap_max={gap_max}, chunks_max={chunks_max}"
        );
        println!("METRIC plusminus_active_solinas_stress_samples={samples}");
        println!("METRIC plusminus_active_solinas_stress_max_scale={max_scale}");
        println!("METRIC plusminus_active_solinas_stress_max_steps={max_steps}");
        println!("METRIC plusminus_active_solinas_stress_projected_p99={projected_p99}");
        println!("METRIC plusminus_active_solinas_stress_projected_p999={projected_p999}");
        println!("METRIC plusminus_active_solinas_stress_projected_max={projected_max}");
        println!("METRIC plusminus_active_solinas_stress_gap_max_to_2700k={gap_max}");
        println!("METRIC plusminus_active_solinas_stress_chunks_max={chunks_max}");
        assert!(gap_max < 0, "active-chain plus-minus Solinas sample max exceeds Google target");
    }

    #[test]
    fn plusminus_solinas_scale_stress_max_stays_below_google() {
        // Larger deterministic stress for the all-in plus-minus cost model:
        // scaled integer add/shift, unary generator, ordering, and Solinas
        // history-carry scale correction.  Still not a proof, but checks that
        // p999/max tails are not close to the 2.7M boundary.
        let p = SECP256K1_P;
        let samples = 32_768usize;
        let mut rng = 0x51a5_6635_2026_0001u64;
        let mut traces = Vec::with_capacity(samples);
        let mut max_scale = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            let unary: usize = ks.iter().sum();
            let steps = ks.len();
            max_scale = max_scale.max(unary);
            traces.push((unary, steps));
        }
        let (scale_dp, chunk_dp) = solinas_history_carry_scale_dp_for_plusminus(max_scale);
        let cmp_ccx = compare_cost_for_plusminus(256);
        let cswap_ccx = cswap_lanes_cost_for_plusminus(&[256, 257]);
        let gen_ccx = trailing_zero_unary_generator_cost_for_plusminus(256);
        let cint_add_ccx = controlled_integer_add_cost_for_plusminus(257);
        let cshift_ccx = controlled_left_shift_cost_for_plusminus(257);
        let mut projected = Vec::with_capacity(samples);
        let mut chunks = Vec::with_capacity(samples);
        for (unary, steps) in traces {
            let step_tax = gen_ccx + cmp_ccx + cswap_ccx;
            let one_div = 2 * (steps * cint_add_ccx + unary * cshift_ccx) + steps * step_tax + scale_dp[unary];
            projected.push(642_716usize + 2 * one_div);
            chunks.push(chunk_dp[unary]);
        }
        projected.sort_unstable();
        chunks.sort_unstable();
        let p99 = samples * 99 / 100;
        let p999 = samples * 999 / 1000;
        let projected_p99 = projected[p99];
        let projected_p999 = projected[p999];
        let projected_max = *projected.last().unwrap();
        let gap_max = projected_max as isize - 2_700_000isize;
        let chunks_max = *chunks.last().unwrap();
        eprintln!(
            "plus-minus Solinas all-in stress: samples={samples}, max_scale={max_scale}, projected_p99={projected_p99}, projected_p999={projected_p999}, projected_max={projected_max}, gap_max={gap_max}, chunks_max={chunks_max}"
        );
        println!("METRIC plusminus_solinas_stress_samples={samples}");
        println!("METRIC plusminus_solinas_stress_max_scale={max_scale}");
        println!("METRIC plusminus_solinas_stress_projected_p99={projected_p99}");
        println!("METRIC plusminus_solinas_stress_projected_p999={projected_p999}");
        println!("METRIC plusminus_solinas_stress_projected_max={projected_max}");
        println!("METRIC plusminus_solinas_stress_gap_max_to_2700k={gap_max}");
        println!("METRIC plusminus_solinas_stress_chunks_max={chunks_max}");
        assert!(gap_max < 0, "sample max exceeds Google target in plus-minus Solinas model");
    }

    fn smag_shl_for_plusminus_test(x: SignedMagU512ForHalfGcdTest, k: usize) -> SignedMagU512ForHalfGcdTest {
        smag_for_halfgcd_test(x.neg, x.mag << k)
    }

    fn smag_shr_exact_for_plusminus_test(x: SignedMagU512ForHalfGcdTest, k: usize) -> Option<SignedMagU512ForHalfGcdTest> {
        if k == 0 {
            return Some(x);
        }
        let low_mask = (U512::from(1u64) << k) - U512::from(1u64);
        if !(x.mag & low_mask).is_zero() {
            None
        } else {
            Some(smag_for_halfgcd_test(x.neg, x.mag >> k))
        }
    }

    fn plusminus_scaled_slack_deficit_for_divisor(x: U256, p: U256) -> (isize, usize, usize) {
        assert!(!x.is_zero());
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
        let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
        let mut history = initial_twos;
        let mut max_deficit = history as isize;
        let mut max_used = 0usize;
        let mut steps = 0usize;
        if u < v {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut cu, &mut cv);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
            let cv_scaled = smag_for_halfgcd_test(cv.neg, cv.mag << k);
            history += k;
            steps += 1;
            if v >= d {
                u = v;
                v = d;
                cu = cv_scaled;
                cv = cd;
            } else {
                u = d;
                cu = cd;
                cv = cv_scaled;
            }
            let coeff_bits = |z: SignedMagU512ForHalfGcdTest| -> usize {
                if z.mag.is_zero() { 1 } else { 1 + u512_bit_len_for_halfgcd_test(z.mag) }
            };
            let used = u512_bit_len_for_halfgcd_test(u)
                + u512_bit_len_for_halfgcd_test(v)
                + coeff_bits(cu)
                + coeff_bits(cv);
            max_used = max_used.max(used);
            let slack = 4isize * 256isize - used as isize;
            max_deficit = max_deficit.max(history as isize - slack);
        }
        (max_deficit, max_used, steps)
    }

    #[test]
    fn plusminus_scaled_integer_state_accounting_needs_history_fusion() {
        // Correct the optimistic state accounting for the scaled-integer route.
        // tx/ty can host one denominator and one coefficient, but an extended
        // plus-minus replay still appears to need the other denominator, the
        // other coefficient, and the unary k history unless history is consumed
        // online or packed into shrinking state.  This is the scratch analogue
        // of the direction/parser blocker.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0xacc0_6635_7a7e_5eedu64;
        let mut unary_payloads = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let ks = plusminus_k_sequence_for_divisor(x, p);
            unary_payloads.push(ks.iter().sum::<usize>());
        }
        unary_payloads.sort_unstable();
        let p99 = samples * 99 / 100;
        let unary_p99 = unary_payloads[p99];
        let unary_max = *unary_payloads.last().unwrap();
        let second_denominator = 256usize;
        let second_coefficient = 256usize;
        let scratch_p99 = second_denominator + second_coefficient + unary_p99;
        let scratch_max = second_denominator + second_coefficient + unary_max;
        let over_google_max = scratch_max as isize - 663isize;
        eprintln!(
            "plus-minus scaled integer state accounting: unary_p99={unary_p99}, unary_max={unary_max}, scratch_p99={scratch_p99}, scratch_max={scratch_max}, over_google_max={over_google_max}"
        );
        println!("METRIC plusminus_scaled_state_unary_p99={unary_p99}");
        println!("METRIC plusminus_scaled_state_unary_max={unary_max}");
        println!("METRIC plusminus_scaled_state_scratch_p99={scratch_p99}");
        println!("METRIC plusminus_scaled_state_scratch_max={scratch_max}");
        println!("METRIC plusminus_scaled_state_over_google_max_bits={over_google_max}");
        assert!(over_google_max > 0, "plus-minus scaled state would fit without history fusion; revisit architecture");
    }

    fn plusminus_scaled_lane_history_trace_for_divisor(x: U256, p: U256) -> Vec<([usize; 4], usize)> {
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
        let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
        let mut history = initial_twos;
        let coeff_bits = |z: SignedMagU512ForHalfGcdTest| -> usize {
            if z.mag.is_zero() { 1 } else { 1 + u512_bit_len_for_halfgcd_test(z.mag) }
        };
        let mut out = Vec::new();
        if u < v {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut cu, &mut cv);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
            let cv_scaled = smag_for_halfgcd_test(cv.neg, cv.mag << k);
            history += k;
            if v >= d {
                u = v;
                v = d;
                cu = cv_scaled;
                cv = cd;
            } else {
                u = d;
                cu = cd;
                cv = cv_scaled;
            }
            out.push(([
                u512_bit_len_for_halfgcd_test(u),
                u512_bit_len_for_halfgcd_test(v),
                coeff_bits(cu),
                coeff_bits(cv),
            ], history));
        }
        out
    }

    fn plusminus_scaled_lane_history_ambig_trace_for_divisor(x: U256, p: U256) -> Vec<([usize; 4], usize, usize)> {
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
        let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
        let mut history = initial_twos;
        let mut ambiguous_dirs = 0usize;
        let coeff_bits = |z: SignedMagU512ForHalfGcdTest| -> usize {
            if z.mag.is_zero() { 1 } else { 1 + u512_bit_len_for_halfgcd_test(z.mag) }
        };
        let div_by_pow2 = |z: SignedMagU512ForHalfGcdTest, k: usize| -> bool {
            k == 0 || z.mag.is_zero() || z.mag.trailing_zeros() as usize >= k
        };
        let mut out = Vec::new();
        if u < v {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut cu, &mut cv);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
            let cv_scaled = smag_for_halfgcd_test(cv.neg, cv.mag << k);
            // Simple local reverse rule: for k>0, cv_scaled is divisible by
            // 2^k.  If cd is not, the ordered output reveals which coefficient
            // lane was scaled; otherwise a persistent direction bit is needed.
            if k == 0 || div_by_pow2(cd, k) {
                ambiguous_dirs += 1;
            }
            history += k;
            if v >= d {
                u = v;
                v = d;
                cu = cv_scaled;
                cv = cd;
            } else {
                u = d;
                cu = cd;
                cv = cv_scaled;
            }
            out.push(([
                u512_bit_len_for_halfgcd_test(u),
                u512_bit_len_for_halfgcd_test(v),
                coeff_bits(cu),
                coeff_bits(cv),
            ], history, ambiguous_dirs));
        }
        out
    }

    fn plusminus_scaled_used_history_trace_for_divisor(x: U256, p: U256) -> Vec<(usize, usize)> {
        let mut u = u512_from_u256_for_halfgcd_test(p);
        let mut v = u512_from_u256_for_halfgcd_test(x);
        let initial_twos = x.trailing_zeros() as usize;
        v >>= initial_twos;
        let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
        let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
        let mut history = initial_twos;
        let coeff_bits = |z: SignedMagU512ForHalfGcdTest| -> usize {
            if z.mag.is_zero() { 1 } else { 1 + u512_bit_len_for_halfgcd_test(z.mag) }
        };
        let mut out = Vec::new();
        if u < v {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut cu, &mut cv);
        }
        while u != v {
            let mut d = u - v;
            let k = d.trailing_zeros() as usize;
            d >>= k;
            let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
            let cv_scaled = smag_for_halfgcd_test(cv.neg, cv.mag << k);
            history += k;
            if v >= d {
                u = v;
                v = d;
                cu = cv_scaled;
                cv = cd;
            } else {
                u = d;
                cu = cd;
                cv = cv_scaled;
            }
            let used = u512_bit_len_for_halfgcd_test(u)
                + u512_bit_len_for_halfgcd_test(v)
                + coeff_bits(cu)
                + coeff_bits(cv);
            out.push((used, history));
        }
        out
    }

    #[test]
    fn plusminus_scaled_state_slack_may_pack_unary_history() {
        // Try the obvious fusion escape for the 914q naive state: the four live
        // 256-bit lanes (two denominators + two scaled coefficients, with tx/ty
        // supplying two of them) develop high-bit slack as the GCD values shrink.
        // If the consumed unary history fits in that slack, only the two extra
        // lanes (≈512 scratch) plus a small sidecar are needed.  This is still a
        // packing model, not a reversible circuit.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x51ac_6635_9acc_5eedu64;
        let mut deficits = Vec::with_capacity(samples);
        let mut useds = Vec::with_capacity(samples);
        let mut steps_v = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let (deficit, used, steps) = plusminus_scaled_slack_deficit_for_divisor(x, p);
            deficits.push(deficit.max(0) as usize);
            useds.push(used);
            steps_v.push(steps);
        }
        deficits.sort_unstable();
        useds.sort_unstable();
        steps_v.sort_unstable();
        let p99 = samples * 99 / 100;
        let deficit_p99 = deficits[p99];
        let deficit_max = *deficits.last().unwrap();
        let used_max = *useds.last().unwrap();
        let scratch_p99 = 512 + deficit_p99;
        let scratch_max = 512 + deficit_max;
        let over_google_max = scratch_max as isize - 663isize;
        eprintln!(
            "plus-minus scaled slack packing: deficit_p99={deficit_p99}, deficit_max={deficit_max}, scratch_p99={scratch_p99}, scratch_max={scratch_max}, used_max={used_max}, over_google_max={over_google_max}"
        );
        println!("METRIC plusminus_scaled_slack_deficit_p99={deficit_p99}");
        println!("METRIC plusminus_scaled_slack_deficit_max={deficit_max}");
        println!("METRIC plusminus_scaled_slack_scratch_p99={scratch_p99}");
        println!("METRIC plusminus_scaled_slack_scratch_max={scratch_max}");
        println!("METRIC plusminus_scaled_slack_used_max={used_max}");
        println!("METRIC plusminus_scaled_slack_over_google_max_bits={over_google_max}");
        assert!(scratch_p99 <= 663, "slack-packed plus-minus misses even p99 Google scratch");
    }

    #[test]
    fn plusminus_scaled_slack_has_public_step_envelope() {
        // The previous slack model used per-sample live bitlengths.  A circuit
        // cannot cheaply know arbitrary bitlength boundaries, so try a stronger
        // public schedule: for each step index, reserve enough live width for
        // the worst sampled state at that index, and ask whether the consumed
        // unary history still fits in the remaining high bits.  If this holds,
        // slack packing may be scheduled by loop index rather than by a dense
        // bitlength oracle.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x57ee_6635_5ced_5eedu64;
        let mut max_used_by_step: Vec<usize> = Vec::new();
        let mut max_hist_by_step: Vec<usize> = Vec::new();
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            for (i, (used, hist)) in plusminus_scaled_used_history_trace_for_divisor(x, p).into_iter().enumerate() {
                if i == max_used_by_step.len() {
                    max_used_by_step.push(0);
                    max_hist_by_step.push(0);
                }
                max_used_by_step[i] = max_used_by_step[i].max(used);
                max_hist_by_step[i] = max_hist_by_step[i].max(hist);
            }
        }
        let mut max_deficit = 0isize;
        let mut worst_step = 0usize;
        for i in 0..max_used_by_step.len() {
            let slack = 4isize * 256isize - max_used_by_step[i] as isize;
            let deficit = max_hist_by_step[i] as isize - slack;
            if deficit > max_deficit {
                max_deficit = deficit;
                worst_step = i;
            }
        }
        let max_deficit_u = max_deficit.max(0) as usize;
        let scratch = 512 + max_deficit_u;
        let over_google = scratch as isize - 663isize;
        let max_steps = max_used_by_step.len();
        eprintln!(
            "plus-minus public-step slack envelope: max_steps={max_steps}, worst_step={worst_step}, deficit={max_deficit_u}, scratch={scratch}, over_google={over_google}, used_at_worst={}, hist_at_worst={} ",
            max_used_by_step[worst_step], max_hist_by_step[worst_step]
        );
        println!("METRIC plusminus_scaled_public_slack_max_steps={max_steps}");
        println!("METRIC plusminus_scaled_public_slack_worst_step={worst_step}");
        println!("METRIC plusminus_scaled_public_slack_deficit={max_deficit_u}");
        println!("METRIC plusminus_scaled_public_slack_scratch={scratch}");
        println!("METRIC plusminus_scaled_public_slack_over_google_bits={over_google}");
        assert!(scratch <= 663, "public step-index slack envelope misses Google scratch");
    }

    #[test]
    fn plusminus_scaled_public_lane_envelope_has_packable_slack() {
        // Strengthen total public slack to lane-specific public slack.  For each
        // step index, take the worst sampled width of each of the four live
        // lanes.  History must fit in the sum of the four high-bit slack bands;
        // the single-lane metric says whether a trivial one-lane stack would be
        // enough or whether a multi-lane fixed map is needed.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x1a4e_6635_51ac_5eedu64;
        let mut max_lane_by_step: Vec<[usize; 4]> = Vec::new();
        let mut max_hist_by_step: Vec<usize> = Vec::new();
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            for (i, (lanes, hist)) in plusminus_scaled_lane_history_trace_for_divisor(x, p).into_iter().enumerate() {
                if i == max_lane_by_step.len() {
                    max_lane_by_step.push([0; 4]);
                    max_hist_by_step.push(0);
                }
                for j in 0..4 {
                    max_lane_by_step[i][j] = max_lane_by_step[i][j].max(lanes[j]);
                }
                max_hist_by_step[i] = max_hist_by_step[i].max(hist);
            }
        }
        let mut total_deficit = 0isize;
        let mut single_deficit = 0isize;
        let mut worst_step = 0usize;
        for i in 0..max_lane_by_step.len() {
            let slack_sum: isize = max_lane_by_step[i].iter().map(|&w| 256isize - w as isize).sum();
            let max_lane_slack: isize = max_lane_by_step[i].iter().map(|&w| 256isize - w as isize).max().unwrap();
            let d_total = max_hist_by_step[i] as isize - slack_sum;
            let d_single = max_hist_by_step[i] as isize - max_lane_slack;
            if d_total > total_deficit {
                total_deficit = d_total;
                worst_step = i;
            }
            single_deficit = single_deficit.max(d_single);
        }
        let total_deficit_u = total_deficit.max(0) as usize;
        let single_deficit_u = single_deficit.max(0) as usize;
        let scratch = 512 + total_deficit_u;
        let over_google = scratch as isize - 663isize;
        eprintln!(
            "plus-minus lane slack envelope: steps={}, total_deficit={total_deficit_u}, single_lane_deficit={single_deficit_u}, scratch={scratch}, over_google={over_google}, worst_step={worst_step}",
            max_lane_by_step.len()
        );
        println!("METRIC plusminus_scaled_lane_slack_steps={}", max_lane_by_step.len());
        println!("METRIC plusminus_scaled_lane_slack_total_deficit={total_deficit_u}");
        println!("METRIC plusminus_scaled_lane_slack_single_deficit={single_deficit_u}");
        println!("METRIC plusminus_scaled_lane_slack_scratch={scratch}");
        println!("METRIC plusminus_scaled_lane_slack_over_google_bits={over_google}");
        assert_eq!(total_deficit_u, 0, "lane-specific public slack does not fit history");
    }

    #[test]
    fn plusminus_scaled_public_lane_envelope_with_direction_bits() {
        // The in-place circuit integration found that the ordering bit cannot
        // be uncomputed after swapping the compared lanes unless a reverse
        // local-recovery circuit is supplied.  Charge one direction bit per
        // plus-minus step as history and ask whether the public lane slack still
        // fits the Google 663 scratch envelope.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0xd1ec_6635_b17d_5eedu64;
        let mut max_lane_by_step: Vec<[usize; 4]> = Vec::new();
        let mut max_hist_by_step: Vec<usize> = Vec::new();
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            for (i, (lanes, hist)) in plusminus_scaled_lane_history_trace_for_divisor(x, p).into_iter().enumerate() {
                if i == max_lane_by_step.len() {
                    max_lane_by_step.push([0; 4]);
                    max_hist_by_step.push(0);
                }
                for j in 0..4 {
                    max_lane_by_step[i][j] = max_lane_by_step[i][j].max(lanes[j]);
                }
                max_hist_by_step[i] = max_hist_by_step[i].max(hist + i + 1);
            }
        }
        let mut total_deficit = 0isize;
        let mut worst_step = 0usize;
        for i in 0..max_lane_by_step.len() {
            let slack_sum: isize = max_lane_by_step[i].iter().map(|&w| 256isize - w as isize).sum();
            let d_total = max_hist_by_step[i] as isize - slack_sum;
            if d_total > total_deficit {
                total_deficit = d_total;
                worst_step = i;
            }
        }
        let total_deficit_u = total_deficit.max(0) as usize;
        let scratch = 512 + total_deficit_u;
        let over_google = scratch as isize - 663isize;
        eprintln!(
            "plus-minus lane slack with direction bits: steps={}, total_deficit={total_deficit_u}, scratch={scratch}, over_google={over_google}, worst_step={worst_step}, hist_dir_at_worst={}",
            max_lane_by_step.len(), max_hist_by_step[worst_step]
        );
        println!("METRIC plusminus_scaled_dir_slack_steps={}", max_lane_by_step.len());
        println!("METRIC plusminus_scaled_dir_slack_total_deficit={total_deficit_u}");
        println!("METRIC plusminus_scaled_dir_slack_scratch={scratch}");
        println!("METRIC plusminus_scaled_dir_slack_over_google_bits={over_google}");
        assert!(over_google > 0 && over_google < 32, "unexpected direction-bit slack result");
    }

    #[test]
    fn plusminus_scaled_public_lane_envelope_with_ambiguous_direction_bits() {
        // Try to avoid storing every ordering bit.  Given k, exactly one output
        // coefficient lane should be divisible by 2^k unless cd=cu-cv is also
        // divisible (or k=0).  Only those ambiguous steps need a persistent
        // direction bit; all other directions are locally recoverable from the
        // ordered coefficient lanes.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0xa6d1_6635_b17d_5eedu64;
        let mut max_lane_by_step: Vec<[usize; 4]> = Vec::new();
        let mut max_hist_by_step: Vec<usize> = Vec::new();
        let mut max_ambig_by_step: Vec<usize> = Vec::new();
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            for (i, (lanes, hist, ambig)) in plusminus_scaled_lane_history_ambig_trace_for_divisor(x, p).into_iter().enumerate() {
                if i == max_lane_by_step.len() {
                    max_lane_by_step.push([0; 4]);
                    max_hist_by_step.push(0);
                    max_ambig_by_step.push(0);
                }
                for j in 0..4 {
                    max_lane_by_step[i][j] = max_lane_by_step[i][j].max(lanes[j]);
                }
                max_hist_by_step[i] = max_hist_by_step[i].max(hist + ambig);
                max_ambig_by_step[i] = max_ambig_by_step[i].max(ambig);
            }
        }
        let mut total_deficit = 0isize;
        let mut worst_step = 0usize;
        for i in 0..max_lane_by_step.len() {
            let slack_sum: isize = max_lane_by_step[i].iter().map(|&w| 256isize - w as isize).sum();
            let d_total = max_hist_by_step[i] as isize - slack_sum;
            if d_total > total_deficit {
                total_deficit = d_total;
                worst_step = i;
            }
        }
        let total_deficit_u = total_deficit.max(0) as usize;
        let scratch = 512 + total_deficit_u;
        let over_google = scratch as isize - 663isize;
        let max_ambig = max_ambig_by_step.iter().copied().max().unwrap_or(0);
        eprintln!(
            "plus-minus lane slack with ambiguous direction bits: steps={}, max_ambig={max_ambig}, total_deficit={total_deficit_u}, scratch={scratch}, over_google={over_google}, worst_step={worst_step}, hist_ambig_at_worst={}",
            max_lane_by_step.len(), max_hist_by_step[worst_step]
        );
        println!("METRIC plusminus_scaled_ambig_dir_steps={}", max_lane_by_step.len());
        println!("METRIC plusminus_scaled_ambig_dir_max_bits={max_ambig}");
        println!("METRIC plusminus_scaled_ambig_dir_total_deficit={total_deficit_u}");
        println!("METRIC plusminus_scaled_ambig_dir_scratch={scratch}");
        println!("METRIC plusminus_scaled_ambig_dir_over_google_bits={over_google}");
        assert!(scratch <= 663, "ambiguous-only direction history misses Google scratch");
    }

    #[test]
    fn plusminus_scaled_public_packing_map_moves_are_clifford_only() {
        // Build a deterministic slot map from the lane envelope: at each public
        // step, history bit j lives in the j-th available high slack slot across
        // all lanes.  Changing this map between steps only requires SWAP/CX
        // movement of already-classical-looking history qubits (Clifford cost),
        // but the move count estimates scheduling complexity.
        let p = SECP256K1_P;
        let samples = 8192usize;
        let mut rng = 0x5a10_6635_c11f_f00du64;
        let mut max_lane_by_step: Vec<[usize; 4]> = Vec::new();
        let mut max_hist_by_step: Vec<usize> = Vec::new();
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            for (i, (lanes, hist)) in plusminus_scaled_lane_history_trace_for_divisor(x, p).into_iter().enumerate() {
                if i == max_lane_by_step.len() {
                    max_lane_by_step.push([0; 4]);
                    max_hist_by_step.push(0);
                }
                for j in 0..4 {
                    max_lane_by_step[i][j] = max_lane_by_step[i][j].max(lanes[j]);
                }
                max_hist_by_step[i] = max_hist_by_step[i].max(hist);
            }
        }
        let slots_for = |lanes: [usize; 4], hist: usize| -> Vec<(usize, usize)> {
            let mut slots = Vec::new();
            for lane in 0..4 {
                for bit in lanes[lane]..256 {
                    slots.push((lane, bit));
                }
            }
            assert!(slots.len() >= hist, "not enough slack slots for public packing map");
            slots.truncate(hist);
            slots
        };
        let mut prev: Vec<(usize, usize)> = Vec::new();
        let mut total_moves = 0usize;
        let mut max_step_moves = 0usize;
        let mut max_hist = 0usize;
        for i in 0..max_lane_by_step.len() {
            let cur = slots_for(max_lane_by_step[i], max_hist_by_step[i]);
            max_hist = max_hist.max(cur.len());
            let common = prev.len().min(cur.len());
            let moves = (0..common).filter(|&j| prev[j] != cur[j]).count();
            total_moves += moves;
            max_step_moves = max_step_moves.max(moves);
            prev = cur;
        }
        eprintln!(
            "plus-minus fixed public packing map: steps={}, max_hist={max_hist}, total_slot_moves={total_moves}, max_step_moves={max_step_moves}",
            max_lane_by_step.len()
        );
        println!("METRIC plusminus_scaled_packmap_steps={}", max_lane_by_step.len());
        println!("METRIC plusminus_scaled_packmap_max_history={max_hist}");
        println!("METRIC plusminus_scaled_packmap_total_slot_moves={total_moves}");
        println!("METRIC plusminus_scaled_packmap_max_step_moves={max_step_moves}");
        assert!(max_hist <= 512, "history map unexpectedly exceeds two-lane scratch budget");
    }

    fn smag_to_twos_for_plusminus_test(x: SignedMagU512ForHalfGcdTest, width: usize) -> U512 {
        let modulus = U512::from(1u64) << width;
        let mask = modulus - U512::from(1u64);
        if x.neg && !x.mag.is_zero() {
            (modulus - (x.mag & mask)) & mask
        } else {
            x.mag & mask
        }
    }

    fn twos_to_smag_for_plusminus_test(x: U512, width: usize) -> SignedMagU512ForHalfGcdTest {
        let modulus = U512::from(1u64) << width;
        let mask = modulus - U512::from(1u64);
        let x = x & mask;
        if x.bit(width - 1) {
            smag_for_halfgcd_test(true, (modulus - x) & mask)
        } else {
            smag_for_halfgcd_test(false, x)
        }
    }

    #[test]
    fn plusminus_scaled_full_state_recovers_direction_locally() {
        // The k-sequence alone has huge reverse preimage rank, but the actual
        // scaled DIV state also contains coefficient lanes.  Check the local
        // inverse relation with `(u,v,cu,cv,k)`: in each direction exactly one
        // of the new coefficient lanes should be divisible by 2^k in the way
        // required to reconstruct the old cv.  If this is unique on real
        // trajectories, direction bits need not be stored separately; reverse
        // can derive them from the full live state.
        let p = SECP256K1_P;
        let p512 = u512_from_u256_for_halfgcd_test(p);
        let samples = 512usize;
        let mut rng = 0x6635_d1ec_71a5_c0deu64;
        let mut ambiguous = 0usize;
        let mut total_steps = 0usize;
        let mut max_candidates = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = p512;
            let mut v = u512_from_u256_for_halfgcd_test(x) >> x.trailing_zeros();
            let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
            let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
            if u < v {
                core::mem::swap(&mut u, &mut v);
                core::mem::swap(&mut cu, &mut cv);
            }
            while u != v {
                let old = (u, v, cu, cv);
                let mut d = u - v;
                let k = d.trailing_zeros() as usize;
                d >>= k;
                let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
                let cv_scaled = smag_shl_for_plusminus_test(cv, k);
                if v >= d {
                    u = v;
                    v = d;
                    cu = cv_scaled;
                    cv = cd;
                } else {
                    u = d;
                    cu = cd;
                    cv = cv_scaled;
                }
                let new = (u, v, cu, cv);
                let mut candidates = 0usize;

                // Reverse dir=1 candidate: new=(old_v,d, old_cv<<k, old_cu-old_cv).
                if let Some(old_cv) = smag_shr_exact_for_plusminus_test(new.2, k) {
                    let old_cu = signed_add_for_halfgcd_test(new.3, old_cv);
                    let prev = (new.0 + (new.1 << k), new.0, old_cu, old_cv);
                    if prev.0 <= p512 && prev == old { candidates += 1; }
                }
                // Reverse dir=0 candidate: new=(d,old_v, old_cu-old_cv, old_cv<<k).
                if let Some(old_cv) = smag_shr_exact_for_plusminus_test(new.3, k) {
                    let old_cu = signed_add_for_halfgcd_test(new.2, old_cv);
                    let prev = (new.1 + (new.0 << k), new.1, old_cu, old_cv);
                    if prev.0 <= p512 && prev == old { candidates += 1; }
                }
                max_candidates = max_candidates.max(candidates);
                if candidates != 1 { ambiguous += 1; }
                total_steps += 1;
            }
        }
        eprintln!(
            "plus-minus scaled full-state local direction recovery: samples={samples}, steps={total_steps}, ambiguous={ambiguous}, max_candidates={max_candidates}"
        );
        println!("METRIC plusminus_scaled_direction_samples={samples}");
        println!("METRIC plusminus_scaled_direction_steps={total_steps}");
        println!("METRIC plusminus_scaled_direction_ambiguous_steps={ambiguous}");
        println!("METRIC plusminus_scaled_direction_max_candidates={max_candidates}");
        assert_eq!(ambiguous, 0, "scaled full state does not locally recover direction");
    }

    #[test]
    fn plusminus_scaled_coefficients_fit_finite_twos_complement_width() {
        // Cost probes used a 257-bit two's-complement-ish integer add/shift
        // floor.  Verify that the exact scaled coefficient recurrence matches
        // arithmetic modulo 2^W with no wrap for W=257 on sampled secp traces.
        // This is still classical, but it closes the most obvious overflow
        // objection to the controlled-shift floor.
        const WIDTH: usize = 257;
        let p = SECP256K1_P;
        let mask = (U512::from(1u64) << WIDTH) - U512::from(1u64);
        let samples = 2048usize;
        let mut rng = 0x2575_6635_7005_c0deu64;
        let mut total_steps = 0usize;
        let mut max_mag_bits = 0usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = u512_from_u256_for_halfgcd_test(p);
            let mut v = u512_from_u256_for_halfgcd_test(x) >> x.trailing_zeros();
            let mut cu = smag_for_halfgcd_test(false, U512::ZERO);
            let mut cv = smag_for_halfgcd_test(false, U512::from(1u64));
            let mut cu_w = smag_to_twos_for_plusminus_test(cu, WIDTH);
            let mut cv_w = smag_to_twos_for_plusminus_test(cv, WIDTH);
            if u < v {
                core::mem::swap(&mut u, &mut v);
                core::mem::swap(&mut cu, &mut cv);
                core::mem::swap(&mut cu_w, &mut cv_w);
            }
            while u != v {
                let mut d = u - v;
                let k = d.trailing_zeros() as usize;
                d >>= k;
                let cd = signed_add_for_halfgcd_test(cu, signed_neg_for_halfgcd_test(cv));
                let cv_scaled = smag_shl_for_plusminus_test(cv, k);
                let cd_w = cu_w.wrapping_sub(cv_w) & mask;
                let cv_scaled_w = (cv_w << k) & mask;
                if v >= d {
                    u = v;
                    v = d;
                    cu = cv_scaled;
                    cv = cd;
                    cu_w = cv_scaled_w;
                    cv_w = cd_w;
                } else {
                    u = d;
                    cu = cd;
                    cv = cv_scaled;
                    cu_w = cd_w;
                    cv_w = cv_scaled_w;
                }
                assert_eq!(twos_to_smag_for_plusminus_test(cu_w, WIDTH), cu, "cu wrapped at width {WIDTH}");
                assert_eq!(twos_to_smag_for_plusminus_test(cv_w, WIDTH), cv, "cv wrapped at width {WIDTH}");
                max_mag_bits = max_mag_bits
                    .max(u512_bit_len_for_halfgcd_test(cu.mag))
                    .max(u512_bit_len_for_halfgcd_test(cv.mag));
                total_steps += 1;
            }
        }
        eprintln!("plus-minus finite two's-complement width: samples={samples}, steps={total_steps}, width={WIDTH}, max_mag_bits={max_mag_bits}");
        println!("METRIC plusminus_scaled_twos_width={WIDTH}");
        println!("METRIC plusminus_scaled_twos_samples={samples}");
        println!("METRIC plusminus_scaled_twos_steps={total_steps}");
        println!("METRIC plusminus_scaled_twos_max_mag_bits={max_mag_bits}");
    }

    #[test]
    fn plusminus_raw_k_live_x_parser_recompute_is_gate_dead() {
        // Last live-parser objection for the plus-minus stream: if raw k bits
        // fit only without delimiters, maybe a parser recomputes the odd-GCD
        // prefix from the live denominator around each use.  This optimistic
        // bitlength-width proxy is even larger than for centered Euclid; simple
        // prefix recomputation is not the missing raw parser.
        let p = SECP256K1_P;
        let samples = 1024usize;
        let mut rng = 0x91a7_600d_11fe_c0deu64;
        let mut weights = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x.is_zero() { x = U256::from(1u64); }
            let mut u = u512_from_u256_for_halfgcd_test(p);
            let mut v = u512_from_u256_for_halfgcd_test(x);
            let initial_twos = x.trailing_zeros() as usize;
            v >>= initial_twos;
            if u < v {
                core::mem::swap(&mut u, &mut v);
            }
            let mut prefix = usize_bit_len_for_payload_test(initial_twos) * 256usize;
            let mut total = 0usize;
            while u != v {
                let mut d = u - v;
                let k = d.trailing_zeros() as usize;
                prefix += usize_bit_len_for_payload_test(k) * u512_bit_len_for_halfgcd_test(u).max(1);
                total += 2 * prefix;
                d >>= k;
                if v >= d {
                    u = v;
                    v = d;
                } else {
                    u = d;
                }
            }
            weights.push(total);
        }
        weights.sort_unstable();
        let mean = weights.iter().sum::<usize>() as f64 / samples as f64;
        let p99 = weights[samples * 99 / 100];
        eprintln!("plus-minus raw-k live-x parser recompute weight: mean={mean:.1}, p99={p99}");
        println!("METRIC plusminus_kseq_live_recompute_weight_mean={mean:.3}");
        println!("METRIC plusminus_kseq_live_recompute_weight_p99={p99}");
        assert!(mean > 8_000_000.0);
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

    fn inv_mod_u64_for_power_two_odd(a: u64, k: usize) -> u64 {
        debug_assert!(k < 64);
        debug_assert_eq!(a & 1, 1);
        let mask = (1u64 << k) - 1;
        let mut x = 1u64;
        // Newton iteration modulo 2,4,8,...; enough for k<=32 in these tests.
        for _ in 0..6 {
            x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x))) & mask;
        }
        x & mask
    }

    fn u256_from_low_u512_for_multihalve_test(x: U512) -> U256 {
        let l = x.as_limbs();
        debug_assert_eq!(l[4], 0);
        debug_assert_eq!(l[5], 0);
        debug_assert_eq!(l[6], 0);
        debug_assert_eq!(l[7], 0);
        U256::from_limbs([l[0], l[1], l[2], l[3]])
    }

    fn halve_once_mod_p_for_multihalve_test(x: U256, p: U256) -> U256 {
        let mut wide = u512_from_u256_for_halfgcd_test(x);
        if (x.as_limbs()[0] & 1) != 0 {
            wide += u512_from_u256_for_halfgcd_test(p);
        }
        u256_from_low_u512_for_multihalve_test(wide >> 1)
    }

    fn solinas_direct_multihalve_step_for_test(x: U256, k: usize) -> (U256, u64) {
        debug_assert!((1..=31).contains(&k));
        let p = SECP256K1_P;
        let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
        let mask = (1u64 << k) - 1;
        let c_low = c.as_limbs()[0] & mask;
        let c_inv = inv_mod_u64_for_power_two_odd(c_low, k);
        let x_low = x.as_limbs()[0] & mask;
        let t = ((x_low as u128 * c_inv as u128) as u64) & mask;
        let wide = u512_from_u256_for_halfgcd_test(x)
            + u512_from_u256_for_halfgcd_test(p) * U512::from(t);
        let y = u256_from_low_u512_for_multihalve_test(wide >> k);
        debug_assert!(y < p);
        (y, t)
    }

    fn multihalve_output_quotient_anf_stats(n: usize, p: u16, k: usize, mask: u16) -> (usize, usize) {
        let size = 1usize << n;
        let mut anf = vec![0u8; size];
        for y in 0..p {
            let q = (((1u32 << k) * y as u32) / p as u32) as u16;
            anf[y as usize] = ((q & mask).count_ones() & 1) as u8;
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

    fn solinas_direct_multihalve_history_step_for_test(x: U256, k: usize) -> (U256, u8) {
        let (y, t) = solinas_direct_multihalve_step_for_test(x, k);
        let h = (y >> (256usize - k)).as_limbs()[0] & ((1u64 << k) - 1);
        let e = t.wrapping_sub(h) & ((1u64 << k) - 1);
        assert!(e <= 1, "residual quotient history should be one correction bit, got {e} for k={k}");
        (y, e as u8)
    }

    fn solinas_direct_multihalve_history_reverse_step_for_test(y: U256, e: u8, k: usize) -> U256 {
        let p = SECP256K1_P;
        let h = (y >> (256usize - k)).as_limbs()[0] & ((1u64 << k) - 1);
        let t = h + e as u64;
        assert!(t < (1u64 << k));
        let wide_y = u512_from_u256_for_halfgcd_test(y) << k;
        let wide_tp = u512_from_u256_for_halfgcd_test(p) * U512::from(t);
        assert!(wide_y >= wide_tp, "reverse multihalve numerator underflow");
        let x_wide = wide_y - wide_tp;
        let x = u256_from_low_u512_for_multihalve_test(x_wide);
        assert!(x < p);
        x
    }

    #[test]
    fn solinas_chunk_multihalve_has_sparse_quotient_cleanup_shape() {
        // Direct k-bit halving over p=2^256-c is classically much nicer than a
        // generic divide: choose t from the low k input bits such that
        // x+t*p is divisible by 2^k.  The cleanup quotient is also recoverable
        // from the output as floor(2^k*y/p), which for a Solinas p is basically
        // the high k bits plus a tiny threshold correction.  This keeps the
        // chunked correction-loop idea alive; the remaining work is a fully
        // reversible in-place circuit for computing/clearing t and its small
        // product t*c.
        let p = SECP256K1_P;
        let mut rng = 0x5eed_d1ec_7a1f_0001u64;
        let samples = 1024usize;
        for _ in 0..samples {
            let mut x = rand_u256(&mut rng);
            if x >= p { x %= p; }
            let (y22, t22) = solinas_direct_multihalve_step_for_test(x, 22);
            let mut reference = x;
            for _ in 0..22 {
                reference = halve_once_mod_p_for_multihalve_test(reference, p);
            }
            assert_eq!(y22, reference);
            let numerator: U512 = u512_from_u256_for_halfgcd_test(y22) << 22usize;
            let q_wide: U512 = numerator / u512_from_u256_for_halfgcd_test(p);
            let q = q_wide.as_limbs()[0];
            assert_eq!(q, t22);
        }
        let mut chunked = rand_u256(&mut rng) % p;
        let mut reference = chunked;
        for _ in 0..404 {
            reference = halve_once_mod_p_for_multihalve_test(reference, p);
        }
        let mut remaining = 404usize;
        while remaining > 0 {
            let k = remaining.min(22);
            chunked = solinas_direct_multihalve_step_for_test(chunked, k).0;
            remaining -= k;
        }
        assert_eq!(chunked, reference);

        let (degree, density) = multihalve_output_quotient_anf_stats(14, 16381, 4, 0b1011);
        eprintln!("Solinas multihalve output quotient ANF: n=14,k=4 degree={degree}, density={density}/16384");
        println!("METRIC solinas_multihalve_chunk_bits={}", 22);
        println!("METRIC solinas_multihalve_product_bits={}", 54);
        println!("METRIC solinas_multihalve_output_quotient_degree_n14={degree}");
        println!("METRIC solinas_multihalve_output_quotient_density_n14={density}");
        assert!(degree >= 12);
        assert!(density <= 32);
    }

    #[test]
    fn solinas_chunk_multihalve_history_carry_is_exact_small_history() {
        // Exact threshold recomputation is locally too expensive, but the
        // residual q-high(y) correction is only one bit per chunk. Keeping that
        // bit through the protected body and consuming it in the explicit
        // inverse avoids dense threshold cleanup while staying reversible.
        let p = SECP256K1_P;
        let schedule404: Vec<usize> = (0..56).map(|_| 7usize).chain((0..2).map(|_| 6usize)).collect();
        let schedule401: Vec<usize> = (0..55).map(|_| 7usize).chain((0..2).map(|_| 8usize)).collect();
        assert_eq!(schedule404.iter().sum::<usize>(), 404);
        assert_eq!(schedule401.iter().sum::<usize>(), 401);
        let mut rng = 0x51de_ca77_1ed5_0001u64;
        let samples = 1024usize;
        for schedule in [&schedule404, &schedule401] {
            for _ in 0..samples {
                let x0 = rand_u256(&mut rng) % p;
                let mut y = x0;
                let mut hist = Vec::with_capacity(schedule.len());
                for &k in schedule.iter() {
                    let (next, e) = solinas_direct_multihalve_history_step_for_test(y, k);
                    y = next;
                    hist.push(e);
                }
                let mut reference = x0;
                for _ in 0..schedule.iter().sum::<usize>() {
                    reference = halve_once_mod_p_for_multihalve_test(reference, p);
                }
                assert_eq!(y, reference);
                for (&k, &e) in schedule.iter().rev().zip(hist.iter().rev()) {
                    y = solinas_direct_multihalve_history_reverse_step_for_test(y, e, k);
                }
                assert_eq!(y, x0);
            }
        }
        eprintln!(
            "Solinas multihalve history-carry exact: bits404={}, bits401={}, samples={samples}",
            schedule404.len(),
            schedule401.len()
        );
        println!("METRIC solinas_multihalve_history_bits_404={}", schedule404.len());
        println!("METRIC solinas_multihalve_history_bits_401={}", schedule401.len());
        println!("METRIC solinas_multihalve_history_exact_samples={samples}");
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
