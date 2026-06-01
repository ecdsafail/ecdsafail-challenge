//! (refactor) Mechanically extracted from kaliski.rs. No logic changes.
use super::*;
/// Specialized real forward primitive for the first few guaranteed-bulk
/// Kaliski iterations where `f = 1` and `v_w != 0` are known a priori.
///
/// This keeps the same persistent-state interface as `kaliski_iteration`
/// (notably `m_i` ends in the same value that the generic step would have
/// produced), but drops STEP 0 / `f` handling entirely.
///
/// Not wired into the live inversion path yet: a direct forward-only swap-in
/// attempt did not preserve full point-add correctness, so this remains an
/// experimental helper while the history/backward compatibility conditions are
/// worked out.
pub(crate) fn kaliski_iteration_bulk_prefix3(
    b: &mut B,
    p: U256,
    u: &[QubitId],
    v_w: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m_i: QubitId,
    m_future: &[QubitId],
    iter_idx: usize,
    coeff: Option<(&[QubitId], &[QubitId])>,
    frame: &mut Option<QubitId>,
    is_last: bool,
) {
    // (r,s) cswap boundary-merge is only valid on the default coeff=None channel.
    let merge_rs = coeff.is_none() && kal_cswap_rs_merge_enabled();
    let merge_uv = merge_rs && kal_cswap_uv_merge_enabled();
    let uv_safe_iters = kal_cswap_uv_merge_safe_iters();
    let uv_merge_in = merge_uv && iter_idx < uv_safe_iters;
    let uv_merge_out = merge_uv && !is_last && iter_idx + 1 < uv_safe_iters;
    let uv_frame_in = if uv_merge_in { *frame } else { None };
    let gz = gz_step4_slow();
    let gz_dbl = gz_double_direct();
    let a_f = b.alloc_qubit();
    let b_f = b.alloc_qubit();
    let add_f = b.alloc_qubit();

    let _kal_saved_phase = b.phase;

    // STEP 0 is a no-op on the guaranteed-bulk prefix (v_w != 0 so the
    // is_zero flag is always 0). The forward measurement-uncompute phases of
    // the OR chain are self-cancelling within with_eq_zero_fast, so dropping
    // the call entirely on both forward and backward is consistent.
    let _ = iter_idx;
    b.set_phase("kal_bulk_step1");
    // Specialized STEP 1 for f=1; the generic z HMR scaffold is a self-
    // cancelling noop (alloc-0 + ccx + hmr + matching cz_if) so we skip it.
    b.x(a_f);
    b.cx(u[0], a_f); // a_f = !u0
    b.x(v_w[0]);
    b.ccx(u[0], v_w[0], m_i); // m_i = u0 & !v0
    b.x(v_w[0]);
    if let Some(frame_in) = uv_frame_in {
        // The previous iter's deferred STEP-9 (u,v_w) swap means physical
        // u/v are conditionally exchanged by frame_in. Correct STEP-1 flags
        // to canonical basis by toggling on frame_in & (u0 xor v0).
        b.cx(v_w[0], u[0]);
        if env_flag_enabled("KAL_UV_STEP1_FANOUT", true) {
            let t = b.alloc_qubit();
            b.ccx(frame_in, u[0], t);
            b.cx(t, a_f);
            b.cx(t, m_i);
            let tm = b.alloc_bit();
            b.hmr(t, tm);
            b.cz_if(frame_in, u[0], tm);
            b.free(t);
        } else {
            b.ccx(frame_in, u[0], a_f);
            b.ccx(frame_in, u[0], m_i);
        }
        b.cx(v_w[0], u[0]);
    }
    b.cx(a_f, b_f);
    b.cx(m_i, b_f); // b_f = a_f xor m_i

    b.set_phase("kal_bulk_step2");
    // Late-iter comparator truncation: bitlen(u)+bitlen(v_w) ≤ 2n-iter_idx so
    // high bits are 0 and don't affect u > v_w.
    let cmp_width = (if iter_idx < u.len() {
        u.len()
    } else {
        2 * u.len() - iter_idx
    })
    .min(kal_wtrunc_width(iter_idx, u.len()));
    let l_gt = b.alloc_qubit();
    with_gt(b, &u[..cmp_width], &v_w[..cmp_width], l_gt, |b| {
        if let Some(frame_in) = uv_frame_in {
            // `with_gt` computed physical_gt. In the equality-free early prefix,
            // canonical_gt = physical_gt xor frame_in.
            b.cx(frame_in, l_gt);
        }
        b.x(b_f);
        let t = b.alloc_qubit();
        b.ccx(l_gt, b_f, t);
        b.cx(t, a_f);
        b.cx(t, m_i);
        {
            let tm = b.alloc_bit();
            b.hmr(t, tm);
            b.cz_if(l_gt, b_f, tm);
        }
        b.free(t);
        // add_dummy scaffold (self-cancelling noop) skipped.
        b.x(b_f);
        if let Some(frame_in) = uv_frame_in {
            b.cx(frame_in, l_gt);
        }
    });
    b.free(l_gt);

    b.set_phase("kal_bulk_step3_cswap");
    // Late-iter truncation: bitlen(u)+bitlen(v_w) ≤ 2n-iter_idx (Kaliski invariant).
    let uv_width_step3 = if iter_idx < u.len() {
        u.len()
    } else {
        2 * u.len() - iter_idx
    };
    if let Some(frame_in) = uv_frame_in {
        // Merge previous STEP-9 uv swap with this STEP-3 uv swap. Control is
        // a_{k-1} xor a_k, built transiently in a_f.
        b.cx(frame_in, a_f);
        for j in 0..uv_width_step3 {
            cswap(b, a_f, u[j], v_w[j]);
        }
        b.cx(frame_in, a_f);
    } else {
        for j in 0..uv_width_step3 {
            cswap(b, a_f, u[j], v_w[j]);
        }
    }
    let rs_width_step3 = if iter_idx + 1 < u.len() {
        iter_idx + 1
    } else {
        u.len()
    };
    // (r,s) STEP 3 — merged with the deferred STEP 9 of the previous iteration
    // when merge_rs and an incoming frame parity is present.
    if let (true, Some(frame_in)) = (merge_rs, *frame) {
        // frame_in = a_{k-1} (previous iter's deferred step9 control).
        // Merged cswap control = a_{k-1} ⊕ a_k. Build into a_f (free CX),
        // emit one cswap (width = min(k+1,n) = step9(k-1) width = step3(k)
        // width), then restore a_f to a_k. After: physical = canonical-post
        // step3(k).
        b.cx(frame_in, a_f); // a_f = a_{k-1} ⊕ a_k
        for j in 0..rs_width_step3 {
            cswap(b, a_f, r[j], s[j]);
        }
        b.cx(frame_in, a_f); // a_f = a_k (restored)
        // Reset frame_in (= a_{k-1}) to |0⟩ via the step10 reroute of the
        // previous iter, evaluated on the now-canonical (r,s) with a_k (= a_f)
        // as the select bit (distinct qubit from frame_in → no self-control):
        //   a_{k-1} = NOT(a_k ? r[0] : s[0])
        // frame_in ^= NOT(a_f ? r[0] : s[0]):
        b.cx(s[0], frame_in);
        b.x(frame_in); // frame_in ^= NOT s[0]
        // CCX-ELIM: frame_in ^= a_f&(r[0]^s[0]) — fold 2 CCX -> 1 CCX + free CX
        // (distribute the AND over the XOR using s[0] as a transient pivot on r[0]).
        b.cx(s[0], r[0]);
        b.ccx(a_f, r[0], frame_in);
        b.cx(s[0], r[0]);
        b.free(frame_in); // frame_in now |0⟩
        *frame = None;
    } else {
        for j in 0..rs_width_step3 {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    if let Some((cr, cs)) = coeff {
        b.set_phase("kal_bulk_coeff_step3_cswap");
        coeff_channel_cswap(b, a_f, cr, cs);
    }

    b.set_phase("kal_bulk_step4");
    // Specialized STEP 4 with add_f = !b_f.
    b.x(add_f);
    b.cx(b_f, add_f);
    {
        let n = u.len();
        // Narrow load/sub width to the late-iter bound (same formula as sub_width).
        // Before this fix: load_width = n, sub_width = max(2n-k, n) → load too wide.
        // After: load_width = sub_width = max(2n-iter_idx, n). Saves n CCX/qubits per iter.
        // W-TRUNC: further narrow to the empirical bitlen envelope (min, never wider).
        let load_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        let tmp = b.alloc_qubits(n);
        for i in 0..load_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        // Narrow load/sub width to the late-iter bound.
        // Both tmp and v_w are 256 qubits. Use slice [0..load_width] for each.
        // 9n-floor: carry-BORROW fast Cuccaro — host the n-1 carry register on
        // clean future m_hist bits (restored to |0>), so the STEP-4 binder
        // drops by up to n-1 at FLAT Toffoli.
        if gz {
            sub_nbit_qq_fast_mfut(b, &tmp[..load_width], &v_w[..load_width], m_future);
        } else {
            sub_nbit_qq_fast(b, &tmp[..load_width], &v_w[..load_width]);
        }
        let transform_width = if iter_idx + 1 < n { iter_idx + 1 } else { n };
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        for i in 0..transform_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        let add_width = if iter_idx + 2 < n { iter_idx + 2 } else { n };
        let mut tmp_slice: Vec<QubitId> = tmp[0..transform_width].to_vec();
        let tmp_pad = if add_width > transform_width {
            let q = b.alloc_qubit();
            tmp_slice.push(q);
            Some(q)
        } else {
            None
        };
        let s_slice: Vec<QubitId> = s[0..add_width].to_vec();
        if gz {
            // Late-iter recovery: also borrow u's provably-|0> high bits so the
            // full-width s-add never falls back to slow Cuccaro (flat peak).
            if gz_late_recover() {
                let u_clean = gz_u_clean_high(u, iter_idx);
                add_nbit_qq_fast_mfut_pool(b, &tmp_slice, &s_slice, m_future, u_clean);
            } else {
                add_nbit_qq_fast_mfut(b, &tmp_slice, &s_slice, m_future);
            }
        } else {
            add_nbit_qq_fast(b, &tmp_slice, &s_slice);
        }
        if let Some(q) = tmp_pad {
            b.free(q);
        }
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            if i < transform_width {
                b.cz_if(add_f, r[i], m);
            } else if i < load_width {
                // W-TRUNC: bits in [transform_width, load_width) hold add_f&u.
                b.cz_if(add_f, u[i], m);
            }
            // W-TRUNC: bits >= load_width were never loaded (tmp[i]=0); the HMR
            // of |0⟩ needs no phase correction.
        }
        b.free_vec(&tmp);
    }
    if let Some((cr, cs)) = coeff {
        b.set_phase("kal_bulk_coeff_step4_add");
        coeff_channel_cadd(b, p, cr, cs, add_f);
    }

    b.set_phase("kal_bulk_step5");
    b.x(b_f);
    {
        let sm = b.alloc_bit();
        b.hmr(add_f, sm);
        b.z_if(b_f, sm);
    }
    b.x(b_f);
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);

    b.set_phase("kal_bulk_step6_7_8");
    for i in 0..(u.len() - 1) {
        b.swap(v_w[i], v_w[i + 1]);
    }
    if iter_idx < r_small_threshold() {
        mod_double_no_corr(b, r);
    } else if gz_dbl {
        // 9n-floor: register-free direct const-add double (drops the
        // const-register + carry-register that bind step6_7_8 at 2457).
        mod_double_inplace_direct(b, r, p);
    } else {
        mod_double_inplace_fast(b, r, p);
    }
    if let Some((cr, _cs)) = coeff {
        b.set_phase("kal_bulk_coeff_step8_double");
        coeff_channel_double(b, p, cr);
    }

    b.set_phase("kal_bulk_step9_cswap");
    // Late-iter truncation: same uv-width bound as step3.
    let uv_width_step9 = if iter_idx < u.len() {
        u.len()
    } else {
        2 * u.len() - iter_idx
    };
    if !uv_merge_out {
        for j in 0..uv_width_step9 {
            cswap(b, a_f, u[j], v_w[j]);
        }
    }
    let rs_width_step9 = if iter_idx + 2 < u.len() {
        iter_idx + 2
    } else {
        u.len()
    };
    if merge_rs && !is_last {
        // DEFER the (r,s) STEP 9 cswap: carry a_k as the outgoing frame parity
        // (allocated here, consumed by the next iter's merged step3). This
        // qubit is NOT live during STEP 4 (allocated after step6_7_8, freed at
        // the next step3 before step4) → peak-neutral. a_f (= a_k) is then
        // reset to |0⟩ for free using the frame copy as select.
        let frame_out = b.alloc_qubit();
        b.cx(a_f, frame_out); // frame_out = a_k
        b.cx(frame_out, a_f); // a_f = a_k ^ a_k = 0
        *frame = Some(frame_out);
    } else {
        // Eager (r,s) STEP 9 (edge: last iter, or merge disabled), then STEP 10.
        for j in 0..rs_width_step9 {
            cswap(b, a_f, r[j], s[j]);
        }
        if let Some((cr, cs)) = coeff {
            b.set_phase("kal_bulk_coeff_step9_cswap");
            coeff_channel_cswap(b, a_f, cr, cs);
        }
        // STEP 10: uncompute a via a ^= NOT s[0].
        b.x(s[0]);
        b.cx(s[0], a_f);
        b.x(s[0]);
    }

    b.free(add_f);
    b.free(b_f);
    b.free(a_f);
    b.set_phase(_kal_saved_phase);
}

pub(crate) fn kaliski_iteration(
    b: &mut B,
    p: U256,
    u: &[QubitId],
    v_w: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m_i: QubitId,
    m_future: &[QubitId],
    f: QubitId,
    iter_idx: usize,
    coeff: Option<(&[QubitId], &[QubitId])>,
    frame: &mut Option<QubitId>,
    is_last: bool,
) {
    let n = u.len();
    // (r,s) cswap boundary-merge is only valid on the default coeff=None channel.
    let merge_rs = coeff.is_none() && kal_cswap_rs_merge_enabled();
    let gz = gz_step4_slow();
    let gz_dbl = gz_double_direct();
    // Iter-local flags (zero at iter start and iter end): alloc fresh here
    // so they don't live during body (which sees lower peak by -3 qubits).
    let a_f = b.alloc_qubit();
    let b_f = b.alloc_qubit();
    let add_f = b.alloc_qubit();

    let _kal_saved_phase = b.phase;
    b.set_phase("kal_step0_eqzero");
    // ─── STEP 0: is_zero = (v_w == 0);  m[i] ^= (f AND is_zero);  f ^= m[i] ───
    // Truncated OR chain for late iter: v_w's bits [2n-iter..n-1] are 0
    // (Kaliski invariant), so OR only of low 2n-iter bits suffices.
    // W-TRUNC: further narrow to the empirical bitlen envelope.
    let or_width =
        (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
    with_eq_zero_fast(b, &v_w[0..or_width], add_f, |b| {
        b.ccx(f, add_f, m_i);
    });
    b.cx(m_i, f);

    b.set_phase("kal_step1");
    // ─── STEP 1 ───
    //   a ^= (f=1 AND u[0]=0)
    //   m[i] ^= (f=1 AND a=0 AND v_w[0]=0)  [= f AND u[0] AND NOT v_w[0]]
    //   b ^= a; b ^= m[i]
    //
    // Shared-intermediate trick: compute z = f AND u[0] once into b_f
    // (known 0 here), then derive a_f = f XOR z = f AND NOT u[0] via CX,
    // and update m_i via ccx(z, NOT v_w[0], m_i). Uncompute z, then set
    // b_f to a_f XOR m_i as before. Saves 1 CCX per iter vs mcx2+mcx3.
    b.ccx(f, u[0], b_f); // b_f = f AND u[0] (z)
    b.cx(f, a_f);
    b.cx(b_f, a_f); // a_f = f XOR z = f AND NOT u[0]
    b.x(v_w[0]);
    b.ccx(b_f, v_w[0], m_i); // m_i ^= z AND NOT v_w[0]
    b.x(v_w[0]);
    // Measurement-uncompute z (= f AND u[0]) from b_f: 0 CCX.
    {
        let zm = b.alloc_bit();
        b.hmr(b_f, zm);
        b.cz_if(f, u[0], zm);
    }
    b.cx(a_f, b_f);
    b.cx(m_i, b_f); // b_f = a_f XOR m_i

    // ─── STEP 2: with l = u > v_w: a ^= (f AND l AND ¬b); m_i ^= same.
    // Late-iter: u and v_w have bitlen ≤ 2n-iter, so only compare low 2n-iter bits.
    // W-TRUNC: further narrow to the empirical bitlen envelope.
    let cmp_width =
        (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
    let l_gt = b.alloc_qubit();
    with_gt(b, &u[0..cmp_width], &v_w[0..cmp_width], l_gt, |b| {
        b.x(b_f); // negate polarity of b_f
        b.ccx(f, l_gt, add_f); // add_f = f AND l_gt
                               // Fuse two CCX with same (add_f, b_f) controls: compute once into
                               // a fresh ancilla, fan out via CX, measurement-uncompute. Saves 1 CCX.
        let t = b.alloc_qubit();
        b.ccx(add_f, b_f, t); // t = add_f AND ¬b_f_orig
        b.cx(t, a_f); // a_f ^= t
        b.cx(t, m_i); // m_i ^= t
        {
            let tm = b.alloc_bit();
            b.hmr(t, tm);
            b.cz_if(add_f, b_f, tm);
        }
        b.free(t);
        // Measurement-uncompute add_f (= f AND l_gt): 0 CCX.
        {
            let am = b.alloc_bit();
            b.hmr(add_f, am);
            b.cz_if(f, l_gt, am);
        }
        b.x(b_f);
    });
    b.free(l_gt);

    b.set_phase("kal_step3_cswap");
    // ─── STEP 3: with control(a): swap(u, v_w); swap(r, s) ───
    // Late-iter truncation: Kaliski invariant: bitlen(u) + bitlen(v_w) ≤ 2n-iter,
    // so u[j]=v_w[j]=0 for j >= 2n-iter_idx. Truncate (u,v_w) cswap.
    // Small-iter truncation: max(r,s) ≤ 2^iter_idx, so r[j]=s[j]=0 for j >= iter_idx+1.
    let uv_width = if iter_idx < n { n } else { 2 * n - iter_idx };
    for j in 0..uv_width {
        cswap(b, a_f, u[j], v_w[j]);
    }
    let rs_width_step3 = if iter_idx + 1 < n { iter_idx + 1 } else { n };
    // (r,s) STEP 3 — merged with the deferred STEP 9 of the previous iter when
    // merge_rs and an incoming frame parity is present. (See bulk variant.)
    if let (true, Some(frame_in)) = (merge_rs, *frame) {
        b.cx(frame_in, a_f); // a_f = a_{k-1} ⊕ a_k
        for j in 0..rs_width_step3 {
            cswap(b, a_f, r[j], s[j]);
        }
        b.cx(frame_in, a_f); // a_f = a_k (restored)
        // Reset frame_in (= a_{k-1}) to |0⟩ via prev iter's step10 reroute,
        // a_k (= a_f) as select: frame_in ^= NOT(a_f ? r[0] : s[0]).
        b.cx(s[0], frame_in);
        b.x(frame_in);
        // CCX-ELIM: frame_in ^= a_f&(r[0]^s[0]) — fold 2 CCX -> 1 CCX + free CX.
        b.cx(s[0], r[0]);
        b.ccx(a_f, r[0], frame_in);
        b.cx(s[0], r[0]);
        b.free(frame_in);
        *frame = None;
    } else {
        for j in 0..rs_width_step3 {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    if let Some((cr, cs)) = coeff {
        b.set_phase("kal_coeff_step3_cswap");
        coeff_channel_cswap(b, a_f, cr, cs);
    }

    b.set_phase("kal_step4");
    // ─── STEP 4 ───
    //   add ^= (f=1 AND b=0)
    //   with control(add): v_w -= u; s += r
    //
    // Fused dual controlled sub+add: reuse one tmp register across both ops.
    // Load tmp = add_f AND u, do sub on v_w, then transform tmp to
    // add_f AND r in place (without unloading + reloading) by temporarily
    // XOR'ing r into u and re-applying ccx(add_f, u, tmp), then add tmp to
    // s and unload. Saves n CCX/iter.
    mcx2_polar(b, f, true, b_f, false, add_f);
    {
        let tmp = b.alloc_qubits(n);
        // Load tmp = add_f AND u. Late-iter bound: u[i]=0 for i >= 2n-iter.
        // W-TRUNC: further narrow to the empirical bitlen envelope.
        let load_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        for i in 0..load_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        // Sub v_w -= tmp. Late-iter: both high bits 0, truncate to load_width.
        let tmp_sub_slice: Vec<QubitId> = tmp[0..load_width].to_vec();
        let v_w_sub_slice: Vec<QubitId> = v_w[0..load_width].to_vec();
        if gz {
            sub_nbit_qq_fast_mfut(b, &tmp_sub_slice, &v_w_sub_slice, m_future);
        } else {
            sub_nbit_qq_fast(b, &tmp_sub_slice, &v_w_sub_slice);
        }
        // Transform tmp from "add_f AND u" to "add_f AND r".
        // Small-iter: only the low iter+1 bits of r can be nonzero; the
        // carry slot for s += r is handled by an explicit 0 pad instead of a
        // useless extra CCX on a known-zero r bit.
        // Late-iter: full transform (r unbounded but u high bits 0 so CCX at
        // high bits effectively produces add_f AND r from tmp=0).
        let transform_width = if iter_idx + 1 < n { iter_idx + 1 } else { n };
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        for i in 0..transform_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        // Add s += tmp. Small-iter still needs one extra carry slot above the
        // live r bits, but that top input bit is known 0.
        let add_width = if iter_idx + 2 < n { iter_idx + 2 } else { n };
        let mut tmp_slice: Vec<QubitId> = tmp[0..transform_width].to_vec();
        let tmp_pad = if add_width > transform_width {
            let q = b.alloc_qubit();
            tmp_slice.push(q);
            Some(q)
        } else {
            None
        };
        let s_slice: Vec<QubitId> = s[0..add_width].to_vec();
        if gz {
            // Late-iter recovery: also borrow u's provably-|0> high bits so the
            // full-width s-add never falls back to slow Cuccaro (flat peak).
            if gz_late_recover() {
                let u_clean = gz_u_clean_high(u, iter_idx);
                add_nbit_qq_fast_mfut_pool(b, &tmp_slice, &s_slice, m_future, u_clean);
            } else {
                add_nbit_qq_fast_mfut(b, &tmp_slice, &s_slice, m_future);
            }
        } else {
            add_nbit_qq_fast(b, &tmp_slice, &s_slice);
        }
        if let Some(q) = tmp_pad {
            b.free(q);
        }
        // Unload: bits < transform_width have tmp = add_f AND r;
        // bits [transform_width..load_width) have tmp = add_f AND u (transform skipped, load done);
        // bits >= load_width have tmp = 0 (load skipped).
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            if i < transform_width {
                b.cz_if(add_f, r[i], m);
            } else if i < load_width {
                b.cz_if(add_f, u[i], m);
            }
            // else: tmp[i]=0, no phase correction needed.
        }
        b.free_vec(&tmp);
    }
    if let Some((cr, cs)) = coeff {
        b.set_phase("kal_coeff_step4_add");
        coeff_channel_cadd(b, p, cr, cs, add_f);
    }

    b.set_phase("kal_step5");
    // ─── STEP 5: uncompute add; uncompute b ───
    // Measurement-uncompute add_f = f AND (NOT b_f): 0 CCX.
    b.x(b_f);
    {
        let sm = b.alloc_bit();
        b.hmr(add_f, sm);
        b.cz_if(f, b_f, sm);
    }
    b.x(b_f);
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);

    b.set_phase("kal_step6_7_8");
    // ─── STEP 6: v_w := v_w / 2 (shift right by 1). Unconditional swap chain.
    // Invariant: v_w[0]=0 before this step whether f=1 (STEP 4 made v_w even)
    // or f=0 (algorithm terminated with v_w=0). Unconditional shift of 0 is 0.
    // Saves 255 CCX per iter vs cswap-controlled version.
    let _ = f;
    for i in 0..(n - 1) {
        b.swap(v_w[i], v_w[i + 1]);
    }

    // ─── STEP 7 + 8: r := 2*r mod p ───────────────────────────────────
    // For iter_idx < r_small_threshold(), r's top bit is guaranteed 0 (since
    // max(r,s) ≤ 2^iter_idx by induction). mod_double's Solinas correction
    // is identity; a plain shift suffices. Saves ~255 CCX per small iter.
    if iter_idx < r_small_threshold() {
        mod_double_no_corr(b, r);
    } else if gz_dbl {
        mod_double_inplace_direct(b, r, p);
    } else {
        mod_double_inplace_fast(b, r, p);
    }
    if let Some((cr, _cs)) = coeff {
        b.set_phase("kal_coeff_step8_double");
        coeff_channel_double(b, p, cr);
    }

    b.set_phase("kal_step9_cswap");
    // ─── STEP 9: with control(a): swap(u, v_w); swap(r, s) (again) ───
    // Late-iter (u,v_w) truncation per Kaliski invariant (same as STEP 3).
    // Small-iter (r,s) truncation: after STEP 4 s ≤ 2^{iter+1}, after STEP 7+8 r ≤ 2^{iter+1}.
    let uv_width = if iter_idx < n { n } else { 2 * n - iter_idx };
    for j in 0..uv_width {
        cswap(b, a_f, u[j], v_w[j]);
    }
    let rs_width_step9 = if iter_idx + 2 < n { iter_idx + 2 } else { n };
    if merge_rs && !is_last {
        // DEFER the (r,s) STEP 9: carry a_k as the outgoing frame parity.
        let frame_out = b.alloc_qubit();
        b.cx(a_f, frame_out); // frame_out = a_k
        b.cx(frame_out, a_f); // a_f = 0 (reset via frame copy)
        *frame = Some(frame_out);
    } else {
        for j in 0..rs_width_step9 {
            cswap(b, a_f, r[j], s[j]);
        }
        if let Some((cr, cs)) = coeff {
            b.set_phase("kal_coeff_step9_cswap");
            coeff_channel_cswap(b, a_f, cr, cs);
        }
        // ─── STEP 10: uncompute a via `a ^= NOT s[0]` ───
        b.x(s[0]);
        b.cx(s[0], a_f);
        b.x(s[0]);
    }

    // Free iter-local flags (all at 0 now).
    b.free(add_f);
    b.free(b_f);
    b.free(a_f);
    b.set_phase(_kal_saved_phase);
}

/// Like `with_eq_zero` but uses measurement-based uncomputation for the
/// backward OR chain (0 Toffoli instead of n-1 CCX). NOT safe inside
/// emit_inverse blocks (uses HMR ops).
pub(crate) fn with_eq_zero_fast<F: FnOnce(&mut B)>(b: &mut B, v: &[QubitId], flag: QubitId, body: F) {
    let n = v.len();
    assert!(n > 0);
    if n == 1 {
        b.x(v[0]);
        b.cx(v[0], flag);
        body(b);
        b.cx(v[0], flag);
        b.x(v[0]);
        return;
    }
    let or_chain: Vec<QubitId> = b.alloc_qubits(n - 1);
    // Forward OR chain (n-1 CCX)
    or_step(b, v[0], v[1], or_chain[0]);
    for i in 1..n - 1 {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    body(b);
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    // Measurement-based uncompute (0 Toffoli)
    for i in (1..n - 1).rev() {
        or_step_uncompute(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    or_step_uncompute(b, v[0], v[1], or_chain[0]);
    b.free_vec(&or_chain);
}

/// Measurement-based uncompute of one or_step: uncomputes
/// `out = x OR y` using HMR + CZ (0 Toffoli).
/// Precondition: out = x OR y (was computed by or_step(x, y, out)).
/// After this: out = 0.
pub(crate) fn or_step_uncompute(b: &mut B, x: QubitId, y: QubitId, out: QubitId) {
    // out currently holds NOT((NOT x) AND (NOT y)) = x OR y.
    // Flip to get the AND value: (NOT x) AND (NOT y).
    b.x(out);
    // Now match the AND controls: flip x and y.
    b.x(x);
    b.x(y);
    let m = b.alloc_bit();
    b.hmr(out, m); // measure; out → 0
    b.cz_if(x, y, m); // phase correction with (NOT x_orig, NOT y_orig) controls
    b.x(y);
    b.x(x);
}

/// Reverse of the specialized `kaliski_iteration_bulk_prefix3` used for the
/// first few guaranteed-bulk nonterminal iterations.
pub(crate) fn kaliski_iteration_bulk_prefix3_backward(
    b: &mut B,
    p: U256,
    u: &[QubitId],
    v_w: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m_i: QubitId,
    m_future: &[QubitId],
    iter_idx: usize,
    frame: &mut Option<QubitId>,
    is_last: bool,
) {
    let n = u.len();
    // (r,s) cswap boundary-merge — bulk backward is always coeff=None.
    let merge_rs = kal_cswap_rs_merge_enabled();
    let merge_uv = merge_rs && kal_cswap_uv_merge_enabled();
    let uv_safe_iters = kal_cswap_uv_merge_safe_iters();
    let uv_merge_in = merge_uv && iter_idx < uv_safe_iters;
    let uv_merge_out = merge_uv && !is_last && iter_idx + 1 < uv_safe_iters;
    let gz = gz_step4_slow();
    let a_f = b.alloc_qubit();
    let b_f = b.alloc_qubit();
    let add_f = b.alloc_qubit();

    let _kal_saved_phase = b.phase;

    // Reverse STEP 10 + STEP 9 (r,s).
    let rs_width_step9 = if iter_idx + 2 < n { iter_idx + 2 } else { n };
    if merge_rs && !is_last {
        // Reverse of forward step9-defer: recreate a_f = a_k from the incoming
        // frame parity, then zero+free the frame.
        b.set_phase("bk_bulk_step9_cswap");
        let frame_in = frame.expect("merged backward expects an incoming frame");
        b.cx(frame_in, a_f); // a_f = a_k
        b.cx(a_f, frame_in); // frame = 0
        b.free(frame_in);
        *frame = None;
    } else {
        // Eager reverse STEP 10 then STEP 9 (r,s) — edge (last iter) / merge off.
        b.set_phase("bk_bulk_step10");
        b.x(s[0]);
        b.cx(s[0], a_f);
        b.x(s[0]);
        b.set_phase("bk_bulk_step9_cswap");
        for j in (0..rs_width_step9).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    // Reverse STEP 9 (u,v) — always eager.
    b.set_phase("bk_bulk_step9_cswap");
    let uv_width_step9 = if iter_idx < n { n } else { 2 * n - iter_idx };
    if !uv_merge_out {
        for j in (0..uv_width_step9).rev() {
            cswap(b, a_f, u[j], v_w[j]);
        }
    }

    // Reverse STEP 8+7 and STEP 6.
    // Bug fix: forward uses mod_double_inplace_fast (with Solinas correction)
    // for iter_idx >= R_SMALL_THRESHOLD, so backward must mirror with
    // mod_halve_inplace_fast to cover the case where r[255]=1 pre-double.
    // Previously unconditional mod_halve_no_corr was a latent bug that
    // happened not to manifest in tested seeds.
    b.set_phase("bk_bulk_step6_7_8");
    if iter_idx < r_small_threshold() {
        mod_halve_no_corr(b, r);
    } else {
        let mut dirty: Vec<QubitId> = u.to_vec();
        dirty.extend_from_slice(v_w);
        mod_halve_inplace_fast_with_dirty(b, r, p, Some(&dirty));
    }
    for i in (0..(n - 1)).rev() {
        b.swap(v_w[i], v_w[i + 1]);
    }

    // Reverse STEP 5.
    b.set_phase("bk_bulk_step5");
    b.cx(a_f, b_f);
    b.cx(m_i, b_f);
    b.x(add_f);
    b.cx(b_f, add_f);

    // Reverse STEP 4.
    b.set_phase("bk_bulk_step4");
    {
        let tmp = b.alloc_qubits(n);
        let load_width = if iter_idx + 1 < n { iter_idx + 1 } else { n };
        for i in 0..load_width {
            b.ccx(add_f, r[i], tmp[i]);
        }
        let sub_width = if iter_idx + 2 < n { iter_idx + 2 } else { n };
        let tmp_sub_slice: Vec<QubitId> = tmp[0..sub_width].to_vec();
        let s_slice: Vec<QubitId> = s[0..sub_width].to_vec();
        if gz {
            // Late-iter recovery: also borrow u's provably-|0> high bits so the
            // full-width s-sub never falls back to slow Cuccaro (flat peak).
            if gz_late_recover() {
                let u_clean = gz_u_clean_high(u, iter_idx);
                sub_nbit_qq_fast_mfut_pool(b, &tmp_sub_slice, &s_slice, m_future, u_clean);
            } else {
                sub_nbit_qq_fast_mfut(b, &tmp_sub_slice, &s_slice, m_future);
            }
        } else if std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1") {
            sub_nbit_qq(b, &tmp_sub_slice, &s_slice);
        } else {
            sub_nbit_qq_fast(b, &tmp_sub_slice, &s_slice);
        }
        // Late-iter denominator bits above 2n-iter_idx are known zero.  The
        // high tmp bits loaded from r only participate in the s-subtraction;
        // they do not need to be transformed into add_f&u or added back into
        // v_w.  This mirrors `kaliski_iteration_backward` and saves one CCX
        // plus two CX per skipped high bit in the bulk reverse tail.
        // W-TRUNC: must match forward load_width exactly (this undoes the
        // forward `v_w -= u` over load_width) → same min envelope.
        let transform_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        for i in 0..transform_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        // After transforming tmp from r to u, high bits of tmp above the
        // late-iter denominator width are known zero.  Truncate the reverse
        // add into v_w just like the generic backward iteration does.
        // W-TRUNC: same envelope as transform_width / forward load_width.
        let add_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        let tmp_add_slice: Vec<QubitId> = tmp[0..add_width].to_vec();
        let v_w_slice: Vec<QubitId> = v_w[0..add_width].to_vec();
        if gz {
            add_nbit_qq_fast_mfut(b, &tmp_add_slice, &v_w_slice, m_future);
        } else if std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1") {
            add_nbit_qq(b, &tmp_add_slice, &v_w_slice);
        } else {
            add_nbit_qq_fast(b, &tmp_add_slice, &v_w_slice);
        }
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            if i < transform_width {
                b.cz_if(add_f, u[i], m);
            } else if i < load_width {
                b.cz_if(add_f, r[i], m);
            }
        }
        b.free_vec(&tmp);
    }
    b.cx(b_f, add_f);
    b.x(add_f);

    // Reverse STEP 3.
    b.set_phase("bk_bulk_step3_cswap");
    let rs_width_step3 = if iter_idx + 1 < n { iter_idx + 1 } else { n };
    // Late-iter truncation mirrors forward step3.
    let uv_width_step3 = if iter_idx < n { n } else { 2 * n - iter_idx };
    // Reverse of the forward (r,s) STEP 3. When merged, recreate the outgoing
    // frame parity (= a_{k-1}) and hand it to the previous (backward-later) iter.
    // Iter 0's forward step3 is an explicit edge (no incoming frame), so its
    // reverse is the plain eager cswap.
    if merge_rs && iter_idx != 0 {
        let frame_out = b.alloc_qubit();
        // Reverse reroute (recreate frame_out = a_{k-1}), a_f = a_k as select.
        // CCX-ELIM: frame_out ^= a_f&(r[0]^s[0]) — fold 2 CCX -> 1 CCX + free CX.
        b.cx(s[0], r[0]);
        b.ccx(a_f, r[0], frame_out);
        b.cx(s[0], r[0]);
        b.x(frame_out);
        b.cx(s[0], frame_out);
        // Reverse the merged cswap: control a_{k-1} ⊕ a_k.
        b.cx(frame_out, a_f); // a_f = a_{k-1} ⊕ a_k
        for j in (0..rs_width_step3).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
        if uv_merge_in {
            for j in (0..uv_width_step3).rev() {
                cswap(b, a_f, u[j], v_w[j]);
            }
        }
        b.cx(frame_out, a_f); // a_f = a_k (restored)
        *frame = Some(frame_out);
    } else {
        for j in (0..rs_width_step3).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    if !(uv_merge_in && merge_rs && iter_idx != 0) {
        for j in (0..uv_width_step3).rev() {
            cswap(b, a_f, u[j], v_w[j]);
        }
    }
    let uv_frame_out = if uv_merge_in { *frame } else { None };

    // Reverse STEP 2.
    b.set_phase("bk_bulk_step2");
    // Mirror forward bulk STEP2 comparator truncation.
    // W-TRUNC: same envelope as forward bulk STEP2 cmp_width.
    let cmp_width =
        (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
    let l_gt = b.alloc_qubit();
    with_gt(b, &u[..cmp_width], &v_w[..cmp_width], l_gt, |b| {
        if let Some(frame_out) = uv_frame_out {
            b.cx(frame_out, l_gt);
        }
        b.x(b_f);
        let t = b.alloc_qubit();
        b.ccx(l_gt, b_f, t);
        b.cx(t, m_i);
        b.cx(t, a_f);
        // Measurement-uncompute t = l_gt & !b_f.  This mirrors the forward
        // bulk step and saves one CCX per reversed bulk iteration.
        let tm = b.alloc_bit();
        b.hmr(t, tm);
        b.cz_if(l_gt, b_f, tm);
        b.free(t);
        b.x(b_f);
        if let Some(frame_out) = uv_frame_out {
            b.cx(frame_out, l_gt);
        }
    });
    b.free(l_gt);

    // Reverse STEP 1.
    b.set_phase("bk_bulk_step1");
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);
    if let Some(frame_out) = uv_frame_out {
        b.cx(v_w[0], u[0]);
        if env_flag_enabled("KAL_UV_STEP1_FANOUT", true) {
            let t = b.alloc_qubit();
            b.ccx(frame_out, u[0], t);
            b.cx(t, m_i);
            b.cx(t, a_f);
            let tm = b.alloc_bit();
            b.hmr(t, tm);
            b.cz_if(frame_out, u[0], tm);
            b.free(t);
        } else {
            b.ccx(frame_out, u[0], m_i);
            b.ccx(frame_out, u[0], a_f);
        }
        b.cx(v_w[0], u[0]);
    }
    b.x(v_w[0]);
    b.ccx(u[0], v_w[0], m_i);
    b.x(v_w[0]);
    b.cx(u[0], a_f);
    b.x(a_f);

    b.free(add_f);
    b.free(b_f);
    b.free(a_f);
    b.set_phase(_kal_saved_phase);
}

/// Reverse of a single kaliski_iteration. Uses measurement-based
/// uncomputation for the OR chain (with_eq_zero) and the step-4 tmp
/// unload, saving ~511 CCX per iteration vs the gate-reversed version.
pub(crate) fn kaliski_iteration_backward(
    b: &mut B,
    p: U256,
    u: &[QubitId],
    v_w: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m_i: QubitId,
    m_future: &[QubitId],
    f: QubitId,
    iter_idx: usize,
    frame: &mut Option<QubitId>,
    is_last: bool,
) {
    let n = u.len();
    // (r,s) cswap boundary-merge — generic backward is always coeff=None here.
    let merge_rs = kal_cswap_rs_merge_enabled();
    let gz = gz_step4_slow();
    // Iter-local flags alloc'd fresh (zero at iter start in the backward
    // direction). They are zeroed and freed at iter end to match forward.
    let a_f = b.alloc_qubit();
    let b_f = b.alloc_qubit();
    let add_f = b.alloc_qubit();

    let _kal_saved_phase = b.phase;
    // ── Reverse STEP 10 + STEP 9 (r,s) ─────────────────────────────────
    let rs_width_step9 = if iter_idx + 2 < n { iter_idx + 2 } else { n };
    if merge_rs && !is_last {
        // Reverse of forward step9-defer: recreate a_f = a_k from incoming frame.
        b.set_phase("bk_step9_cswap");
        let frame_in = frame.expect("merged backward expects an incoming frame");
        b.cx(frame_in, a_f); // a_f = a_k
        b.cx(a_f, frame_in); // frame = 0
        b.free(frame_in);
        *frame = None;
    } else {
        b.set_phase("bk_step10");
        // Reverse STEP 10. Matches forward's gated update.
        b.x(s[0]);
        b.ccx(f, s[0], a_f);
        b.x(s[0]);
        b.set_phase("bk_step9_cswap");
        for j in (0..rs_width_step9).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    let uv_width = if iter_idx < n { n } else { 2 * n - iter_idx };
    b.set_phase("bk_step9_cswap");
    for j in (0..uv_width).rev() {
        cswap(b, a_f, u[j], v_w[j]);
    }

    b.set_phase("bk_step6_7_8");
    // Reverse STEP 8 + 7 ─────────────────────────────────────────────
    // For iter_idx < r_small_threshold(), forward used mod_double_no_corr —
    // r is guaranteed even (bit 0 = 0), so a plain shift-right inverts it.
    if iter_idx < r_small_threshold() {
        mod_halve_no_corr(b, r);
    } else {
        let mut dirty: Vec<QubitId> = u.to_vec();
        dirty.extend_from_slice(v_w);
        mod_halve_inplace_fast_with_dirty(b, r, p, Some(&dirty));
    }

    // ── Reverse STEP 6 (unconditional shift-left) ───────────
    let _ = f;
    for i in (0..(n - 1)).rev() {
        b.swap(v_w[i], v_w[i + 1]);
    }

    b.set_phase("bk_step5");
    // Reverse STEP 5 ─────────────────────────────────────────────────
    b.cx(a_f, b_f);
    b.cx(m_i, b_f);
    mcx2_polar(b, f, true, b_f, false, add_f);

    b.set_phase("bk_step4");
    // Reverse STEP 4 (with measurement uncompute for unload) ─────────
    {
        let tmp = b.alloc_qubits(n);
        // Load tmp = AND(add_f, r). Small-iter: r[i]=0 for i >= iter+1.
        let load_width = if iter_idx + 1 < n { iter_idx + 1 } else { n };
        for i in 0..load_width {
            b.ccx(add_f, r[i], tmp[i]);
        }
        // Reversed (F): sub tmp from s. Small-iter width iter+2.
        let sub_width = if iter_idx + 2 < n { iter_idx + 2 } else { n };
        let tmp_sub_slice: Vec<QubitId> = tmp[0..sub_width].to_vec();
        let s_slice: Vec<QubitId> = s[0..sub_width].to_vec();
        if gz {
            // Late-iter recovery: also borrow u's provably-|0> high bits so the
            // full-width s-sub never falls back to slow Cuccaro (flat peak).
            if gz_late_recover() {
                let u_clean = gz_u_clean_high(u, iter_idx);
                sub_nbit_qq_fast_mfut_pool(b, &tmp_sub_slice, &s_slice, m_future, u_clean);
            } else {
                sub_nbit_qq_fast_mfut(b, &tmp_sub_slice, &s_slice, m_future);
            }
        } else if std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1") {
            sub_nbit_qq(b, &tmp_sub_slice, &s_slice);
        } else {
            sub_nbit_qq_fast(b, &tmp_sub_slice, &s_slice);
        }
        // Reversed (E): transform tmp from AND(add_f,r) → AND(add_f,u).
        // Late-iter: u high bits 0, so transform at those bits: cx(r,u=0)→u=r,
        //   ccx(add_f, u=r, tmp) flips tmp. tmp goes 0 → add_f AND r. Not what we
        //   want (need add_f AND u=0). For late iter, truncate transform to uv_width.
        // W-TRUNC: must match forward load_width exactly (undoes `v_w -= u`).
        let transform_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        for i in 0..transform_width {
            b.ccx(add_f, u[i], tmp[i]);
        }
        for i in 0..transform_width {
            b.cx(r[i], u[i]);
        }
        // Reversed (D): add tmp to v_w. Truncated to uv_width (late iter bound).
        let add_width = transform_width;
        let tmp_add_slice: Vec<QubitId> = tmp[0..add_width].to_vec();
        let v_w_slice: Vec<QubitId> = v_w[0..add_width].to_vec();
        if gz {
            add_nbit_qq_fast_mfut(b, &tmp_add_slice, &v_w_slice, m_future);
        } else if std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1") {
            add_nbit_qq(b, &tmp_add_slice, &v_w_slice);
        } else {
            add_nbit_qq_fast(b, &tmp_add_slice, &v_w_slice);
        }
        // Unload: bits < min(load_width, transform_width) both apply (tmp = add_f AND u after transform).
        // For bits where transform was applied, tmp = add_f AND u. For bits where transform skipped
        // (i >= transform_width), tmp stays at whatever load left it (either add_f AND r or 0).
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            if i < transform_width {
                // Transform applied: tmp = add_f AND u.
                b.cz_if(add_f, u[i], m);
            } else if i < load_width {
                // Load done but transform skipped: tmp = add_f AND r.
                b.cz_if(add_f, r[i], m);
            }
            // else: tmp = 0, no phase.
        }
        b.free_vec(&tmp);
    }
    // Reversed (A): measurement-uncompute add_f = f AND (NOT b_f)
    b.x(b_f);
    {
        let sm = b.alloc_bit();
        b.hmr(add_f, sm);
        b.cz_if(f, b_f, sm);
    }
    b.x(b_f);

    b.set_phase("bk_step3_cswap");
    // Reverse STEP 3 ─────────────────────────────────────────────────
    let rs_width_step3 = if iter_idx + 1 < n { iter_idx + 1 } else { n };
    let uv_width = if iter_idx < n { n } else { 2 * n - iter_idx };
    // Reverse of forward (r,s) STEP 3: recreate the outgoing frame parity when
    // merged. Iter 0 forward step3 was an explicit edge → plain eager reverse.
    if merge_rs && iter_idx != 0 {
        let frame_out = b.alloc_qubit();
        // Reverse reroute (recreate frame_out = a_{k-1}), a_f = a_k as select.
        // CCX-ELIM: frame_out ^= a_f&(r[0]^s[0]) — fold 2 CCX -> 1 CCX + free CX.
        b.cx(s[0], r[0]);
        b.ccx(a_f, r[0], frame_out);
        b.cx(s[0], r[0]);
        b.x(frame_out);
        b.cx(s[0], frame_out);
        b.cx(frame_out, a_f); // a_f = a_{k-1} ⊕ a_k
        for j in (0..rs_width_step3).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
        b.cx(frame_out, a_f); // a_f = a_k (restored)
        *frame = Some(frame_out);
    } else {
        for j in (0..rs_width_step3).rev() {
            cswap(b, a_f, r[j], s[j]);
        }
    }
    for j in (0..uv_width).rev() {
        cswap(b, a_f, u[j], v_w[j]);
    }

    b.set_phase("bk_step2");
    // Reverse STEP 2 (with_gt body is self-inverse) ──────────────────
    // W-TRUNC: same envelope as forward generic STEP2 cmp_width.
    let cmp_width =
        (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
    let l_gt = b.alloc_qubit();
    with_gt(b, &u[0..cmp_width], &v_w[0..cmp_width], l_gt, |b| {
        b.x(b_f);
        b.ccx(f, l_gt, add_f);
        // Fuse two CCX with same (add_f, b_f) controls into one CCX + two CX
        // + measurement uncompute. Saves 1 CCX per backward iter.
        let t = b.alloc_qubit();
        b.ccx(add_f, b_f, t);
        b.cx(t, m_i);
        b.cx(t, a_f);
        {
            let tm = b.alloc_bit();
            b.hmr(t, tm);
            b.cz_if(add_f, b_f, tm);
        }
        b.free(t);
        // Measurement-uncompute add_f = f AND l_gt: 0 CCX.
        {
            let am = b.alloc_bit();
            b.hmr(add_f, am);
            b.cz_if(f, l_gt, am);
        }
        b.x(b_f);
    });
    b.free(l_gt);

    b.set_phase("bk_step1");
    // Reverse STEP 1 ─────────────────────────────────────────────────
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);
    b.ccx(f, u[0], b_f);
    b.x(v_w[0]);
    b.ccx(b_f, v_w[0], m_i);
    b.x(v_w[0]);
    b.cx(b_f, a_f);
    b.cx(f, a_f);
    // Measurement-uncompute z = f AND u[0] from b_f: 0 CCX.
    {
        let zm = b.alloc_bit();
        b.hmr(b_f, zm);
        b.cz_if(f, u[0], zm);
    }

    b.set_phase("bk_step0_eqzero");
    // Reverse STEP 0 (with measurement uncompute of OR chain) ────────
    // Truncated for late iter: only low 2n-iter bits of v_w are possibly nonzero.
    // W-TRUNC: same envelope as forward generic STEP0 or_width.
    b.cx(m_i, f);
    {
        let or_width =
            (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n));
        let nv = or_width;
        if nv == 1 {
            b.x(v_w[0]);
            b.cx(v_w[0], add_f);
            b.ccx(f, add_f, m_i);
            b.cx(v_w[0], add_f);
            b.x(v_w[0]);
        } else {
            let or_chain: Vec<QubitId> = b.alloc_qubits(nv - 1);
            or_step(b, v_w[0], v_w[1], or_chain[0]);
            for i in 1..nv - 1 {
                or_step(b, or_chain[i - 1], v_w[i + 1], or_chain[i]);
            }
            b.x(or_chain[nv - 2]);
            b.cx(or_chain[nv - 2], add_f);
            b.x(or_chain[nv - 2]);
            // Body
            b.ccx(f, add_f, m_i);
            // Uncompute flag
            b.x(or_chain[nv - 2]);
            b.cx(or_chain[nv - 2], add_f);
            b.x(or_chain[nv - 2]);
            // Measurement-based uncompute of OR chain (0 Toffoli)
            for i in (1..nv - 1).rev() {
                or_step_uncompute(b, or_chain[i - 1], v_w[i + 1], or_chain[i]);
            }
            or_step_uncompute(b, v_w[0], v_w[1], or_chain[0]);
            b.free_vec(&or_chain);
        }
    }

    // Free iter-local flags (all at 0 now after backward steps).
    b.free(add_f);
    b.free(b_f);
    b.free(a_f);
    b.set_phase(_kal_saved_phase);
}

pub(crate) fn with_eq_const_fast<F: FnOnce(&mut B)>(
    b: &mut B,
    bits: &[QubitId],
    c: usize,
    flag: QubitId,
    body: F,
) {
    for (i, &q) in bits.iter().enumerate() {
        if ((c >> i) & 1) != 0 {
            b.x(q);
        }
    }
    with_eq_zero_fast(b, bits, flag, body);
    for (i, &q) in bits.iter().enumerate() {
        if ((c >> i) & 1) != 0 {
            b.x(q);
        }
    }
}
