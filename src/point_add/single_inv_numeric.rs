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

use alloy_primitives::U256;

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
        let (mut qx, mut qy) = (0u16, 0u16);
        'outer: for x in 1..p {
            for y in 1..p {
                if is_curve_point_u16_for_phase_test(x, y, p) {
                    qx = x;
                    qy = y;
                    break 'outer;
                }
            }
        }
        let (degree, density) = top_level_measured_input_phase_anf_stats(8, p, qx, qy, 0b1010_0101_0101_1010);
        eprintln!(
            "Top-level old-point MBUC phase ANF: q=({qx},{qy}), degree={degree}, density={density}/65536"
        );
        assert!(degree >= 14);
        assert!(density > 10_000);
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
