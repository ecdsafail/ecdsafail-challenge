use super::{BExt, B};
use crate::circuit::{BitId, QubitId};
use std::cell::Cell;

thread_local! {
    static FFG_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FFG_SHIFTED_SQUARE_PREFIX_SCOPE: Cell<usize> = const { Cell::new(0) };
    static CUCCARO_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static CONST_CHUNK_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static ADD_CONST_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_ffg_call_index() {
    FFG_CALL_INDEX.with(|index| index.set(0));
    CUCCARO_CALL_INDEX.with(|index| index.set(0));
    CONST_CHUNK_CALL_INDEX.with(|index| index.set(0));
    ADD_CONST_CALL_INDEX.with(|index| index.set(0));
}

fn next_ffg_call_index() -> usize {
    FFG_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

pub(super) fn with_shifted_square_ffg_prefix_scope<R>(body: impl FnOnce() -> R) -> R {
    FFG_SHIFTED_SQUARE_PREFIX_SCOPE.with(|scope| {
        let prior = scope.get();
        scope.set(prior + 1);
        let result = body();
        scope.set(prior);
        result
    })
}

fn shifted_square_ffg_prefix_scope_enabled() -> bool {
    std::env::var_os("TLM_SQUARE_SHIFTED_FFG_PREFIX_SKIP").is_some()
        && FFG_SHIFTED_SQUARE_PREFIX_SCOPE.with(|scope| scope.get() > 0)
}

fn next_cuccaro_call_index() -> usize {
    CUCCARO_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_const_chunk_call_index() -> usize {
    CONST_CHUNK_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_add_const_call_index() -> usize {
    ADD_CONST_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn env_index_value(name: &str, index: usize) -> Option<usize> {
    std::env::var(name).ok().and_then(|value| {
        value
            .split(',')
            .filter_map(|item| item.trim().split_once(':'))
            .find_map(|(call, value)| {
                (call.parse::<usize>().ok()? == index)
                    .then(|| value.parse::<usize>().ok())
                    .flatten()
            })
    })
}

fn env_index_list_contains(name: &str, index: usize) -> bool {
    std::env::var(name)
        .ok()
        .map(|value| {
            value
                .split(',')
                .filter_map(|item| item.trim().parse::<usize>().ok())
                .any(|candidate| candidate == index)
        })
        .unwrap_or(false)
}

const FFG_DEAD_HYBRID_CARRY_RANGES: &[(usize, usize, usize)] = &[
    (264, 1, 46),
    (265, 1, 46),
    (266, 1, 46),
    (267, 1, 46),
    (268, 1, 46),
    (271, 1, 46),
    (272, 1, 46),
    (273, 1, 46),
    (274, 1, 46),
    (275, 1, 46),
    (277, 1, 46),
    (278, 1, 46),
    (279, 1, 46),
    (280, 1, 46),
    (281, 1, 46),
    (596, 21, 24),
    (596, 26, 31),
    (596, 42, 46),
    (2, 1, 3),
    (2, 28, 31),
    (598, 1, 3),
    (598, 27, 27),
    (598, 30, 31),
    (3, 1, 3),
    (3, 29, 29),
    (3, 31, 31),
    (340, 1, 5),
    (51, 28, 31),
    (131, 28, 31),
    (198, 28, 31),
    (201, 28, 31),
    (597, 1, 3),
    (597, 29, 29),
    (13, 29, 31),
    (37, 29, 31),
    (50, 29, 31),
    (60, 29, 31),
    (64, 28, 28),
    (64, 30, 31),
    (73, 29, 31),
    (75, 29, 31),
    (80, 29, 31),
    (105, 29, 31),
    (113, 29, 31),
    (115, 28, 28),
    (115, 30, 31),
    (116, 29, 31),
    (119, 29, 31),
    (126, 29, 31),
    (137, 29, 31),
    (139, 28, 29),
    (139, 31, 31),
    (140, 29, 31),
    (147, 29, 31),
    (178, 29, 31),
    (190, 29, 31),
    (199, 29, 31),
    (209, 28, 28),
    (209, 30, 31),
    (284, 29, 31),
    (288, 29, 31),
    (293, 29, 31),
    (295, 29, 31),
    (318, 29, 31),
    (405, 29, 31),
    (409, 28, 28),
    (409, 30, 31),
    (416, 29, 31),
    (424, 28, 28),
    (424, 30, 31),
    (433, 28, 28),
    (433, 30, 31),
    (434, 29, 31),
    (444, 29, 31),
    (464, 29, 31),
    (471, 29, 31),
    (478, 28, 28),
    (478, 30, 31),
    (487, 28, 28),
    (487, 30, 31),
    (498, 29, 31),
    (516, 29, 31),
    (518, 29, 31),
    (548, 28, 28),
    (548, 30, 31),
    (553, 29, 31),
    (559, 29, 31),
    (560, 29, 31),
    (568, 29, 31),
    (570, 29, 31),
    (575, 29, 31),
    (580, 29, 31),
    (586, 29, 31),
    (592, 29, 31),
];

fn ffg_call_has_structurally_dead_hybrid_carry(call_index: usize, bit: usize, phase: &str) -> bool {
    if shifted_square_ffg_prefix_scope_enabled() && bit > 0 {
        return true;
    }
    if std::env::var_os("TLM_FFG_SKIP_TOP_CARRY31").is_some() && bit == 31 {
        return true;
    }
    if std::env::var_os("TLM_FFG_SKIP_TOP_CARRY30").is_some() && bit == 30 {
        return true;
    }
    if std::env::var_os("TLM_FFG_SKIP_INVERSE_MOD_SUB_TOP29").is_some()
        && bit == 29
        && phase == "tlm_apply_inverse_mod_sub_fold"
        && call_index
            <= std::env::var("TLM_FFG_INVERSE_TOP29_MAX_CALL")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(usize::MAX)
    {
        return true;
    }
    if std::env::var_os("TLM_FFG_SKIP_STRUCTURAL_DEAD_CALLS").is_none() {
        return false;
    }
    if std::env::var_os("TLM_FFG_SKIP_EXACT_TOP29_REMAINDER").is_some() {
        let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if FFG_TOP29_REMAINDER_KEYS.binary_search(&key).is_ok() {
            return true;
        }
    }
    FFG_DEAD_HYBRID_CARRY_RANGES
        .iter()
        .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

const FFG_TOP29_REMAINDER_KEYS: &[u32] = &[
    1821, 2333, 3869, 6685, 7197, 7453, 13341, 15901, 19741, 19997, 20253, 22044, 25885, 26397,
    27933, 31517, 32796, 34077, 36125, 36380, 38173, 38941, 40989, 41757, 42525, 44316, 46621,
    50205, 54045, 54557, 68893, 72221, 74781, 79133, 85277, 85789, 86557, 102685, 103453, 104989,
    108061, 110365, 112669, 115741, 117789, 120861, 121117, 123165, 126493, 126749, 127260, 128797,
    129053, 130844, 131101, 137245, 144157, 144669, 147741, 149021, 149533, 151069,
];

const CONST_CHUNK_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (879, 0, 9),
    (880, 0, 8),
    (700, 0, 7),
    (881, 0, 7),
    (887, 0, 7),
    (900, 0, 7),
    (678, 0, 6),
    (686, 0, 6),
    (707, 0, 6),
    (882, 0, 6),
    (894, 0, 6),
    (906, 0, 6),
    (649, 2, 7),
    (654, 2, 7),
    (692, 1, 6),
    (715, 0, 4),
    (715, 7, 7),
    (883, 0, 5),
    (638, 0, 1),
    (638, 3, 5),
    (689, 0, 4),
    (691, 3, 7),
    (706, 2, 2),
    (706, 4, 7),
    (718, 0, 4),
    (884, 0, 4),
    (890, 0, 4),
    (897, 0, 4),
    (903, 0, 4),
    (909, 0, 4),
    (912, 1, 5),
    (919, 0, 4),
    (626, 1, 4),
    (659, 4, 7),
    (666, 0, 3),
    (671, 3, 3),
    (671, 5, 7),
    (672, 0, 3),
    (685, 5, 8),
    (703, 1, 4),
    (704, 0, 3),
    (714, 5, 8),
    (717, 2, 5),
    (891, 0, 3),
    (893, 5, 8),
    (926, 0, 3),
    (941, 0, 3),
    (949, 1, 1),
    (949, 3, 4),
    (949, 6, 6),
    (956, 2, 5),
    (481, 1, 3),
    (485, 1, 3),
    (520, 5, 7),
    (579, 1, 2),
    (579, 4, 4),
    (590, 2, 4),
    (636, 0, 2),
    (644, 4, 5),
    (644, 7, 7),
    (663, 1, 3),
    (665, 4, 4),
    (665, 6, 7),
    (669, 0, 0),
    (669, 2, 3),
    (675, 1, 3),
    (677, 4, 5),
    (677, 7, 7),
    (682, 0, 2),
    (696, 0, 2),
    (710, 0, 1),
    (710, 3, 3),
    (711, 0, 2),
    (723, 0, 2),
    (725, 0, 2),
    (727, 0, 2),
    (729, 0, 2),
    (731, 0, 2),
    (737, 0, 2),
    (739, 0, 2),
    (741, 0, 2),
    (743, 0, 2),
    (745, 0, 2),
    (749, 0, 2),
    (751, 0, 2),
    (753, 0, 2),
    (755, 0, 2),
    (757, 0, 2),
    (896, 3, 5),
    (905, 5, 5),
    (905, 7, 8),
    (911, 5, 7),
    (916, 0, 2),
    (918, 5, 7),
    (923, 0, 2),
    (931, 5, 7),
    (932, 0, 2),
    (935, 1, 3),
    (937, 5, 7),
    (944, 0, 2),
    (962, 2, 4),
    (968, 1, 3),
    (973, 1, 3),
    (978, 1, 3),
    (983, 0, 2),
    (1023, 1, 3),
    (1039, 0, 2),
    (1044, 1, 3),
    (1109, 3, 5),
    (1149, 1, 3),
    (1622, 0, 2),
];

const CONST_CHUNK_REMAINDER_KEYS: &[u32] = &[
    1281, 5376, 5377, 6913, 8449, 8960, 9473, 10496, 15616, 16641, 17153, 19201, 19713, 20736,
    22785, 23809, 27905, 29440, 38656, 40705, 43777, 49921, 51457, 57089, 58113, 59649, 64257,
    66305, 66816, 68865, 70400, 70401, 70913, 77569, 79361, 80898, 99074, 113408, 116224, 117248,
    118016, 118017, 119043, 120064, 120065, 121090, 121091, 122115, 125189, 126467, 127492, 127493,
    128772, 128774, 130308, 130309, 131589, 131590, 135936, 136711, 136960, 136961, 137990, 137991,
    139015, 139264, 140035, 140545, 141315, 141575, 141824, 142855, 143105, 144385, 144386, 145155,
    145415, 145665, 146946, 149761, 149762, 152579, 152581, 153861, 154112, 155137, 156166, 156417,
    157447, 157696, 157697, 158466, 158978, 159746, 161794, 161796, 164354, 165120, 166915, 167680,
    168960, 168961, 169475, 174339, 174848, 174849, 176132, 176133, 177922, 177923, 178432, 178433,
    178952, 179456, 179457, 181506, 181507, 182272, 182273, 185344, 185345, 185856, 185857, 186368,
    186369, 186880, 186881, 187392, 187393, 188928, 188929, 189440, 189441, 189952, 189953, 190464,
    190465, 190976, 190977, 192000, 192001, 192512, 192513, 193024, 193025, 193536, 193537, 194048,
    194049, 196609, 199680, 200192, 200193, 214017, 217089, 223233, 226822, 226824, 227328, 227329,
    230151, 232453, 234242, 234243, 236035, 236806, 236807, 237826, 237827, 240128, 240129, 241414,
    241415, 242435, 244225, 245761, 245762, 247296, 248579, 252419, 252930, 254214, 255491, 255493,
    257027, 257028, 258050, 258566, 259840, 260354, 260355, 263172, 264707, 268546, 269825, 271105,
    271106, 272385, 273154, 273155, 273415, 273664, 274689, 275968, 276743, 277767, 278016, 278789,
    278791, 279812, 279814, 281089, 281348, 281350, 282373, 282374, 285188, 285190, 286723, 286725,
    288003, 288004, 289284, 289285, 290306, 290562, 290564, 291842, 292865, 292868, 295170, 296195,
    297217, 297218, 299265, 299266, 300288, 301312, 302081, 303104, 304897, 305152, 323329, 328194,
    330499, 341761, 345857, 346369, 346881, 347905, 352512, 354049, 359168, 361217, 361729, 362241,
    367873, 368385, 369920, 371457, 374529, 375040, 375041, 376577, 378625, 380672, 380673, 381697,
    385281, 386817, 391937, 393985, 395008, 398593, 400641, 402689, 405761, 407809, 411393, 415488,
    415489, 416001, 416513,
];

fn const_chunk_call_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_CONST_CHUNK_SKIP_STRUCTURAL_DEAD_CALLS").is_none() {
        return false;
    }
    if std::env::var_os("TLM_CONST_CHUNK_SKIP_EXACT_REMAINDER").is_some() {
        let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if CONST_CHUNK_REMAINDER_KEYS.binary_search(&key).is_ok() {
            return true;
        }
    }
    CONST_CHUNK_DEAD_RANGES
        .iter()
        .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn cuccaro_call_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_CUCCARO_SKIP_STRUCTURAL_DEAD_CALLS").is_none() {
        return false;
    }
    match call_index {
        12 | 25 => (1..=127).contains(&bit),
        37 => bit <= 135,
        19 => (1..=127).contains(&bit),
        20 | 26 => bit >= 148,
        13 => bit >= 150,
        21 => matches!(bit, 147 | 148) || (150..=251).contains(&bit),
        27 => (148..=251).contains(&bit),
        22 => bit == 146 || (148..=249).contains(&bit),
        28 => (147..=249).contains(&bit),
        14 => matches!(bit, 150 | 151) || (153..=251).contains(&bit),
        15 => (151..=249).contains(&bit),
        29 => (147..=245).contains(&bit),
        23 => matches!(bit, 148 | 149) || (151..=245).contains(&bit),
        16 => (151..=245).contains(&bit),
        30 => (149..=223).contains(&bit),
        24 => (150..=223).contains(&bit),
        17 => bit == 149 || (152..=223).contains(&bit),
        _ => false,
    }
}

fn add_const_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_ADD_CONST_SKIP_STRUCTURAL_DEAD_CARRIES").is_none() {
        return false;
    }
    call_index == 0 && (bit == 55 || bit >= 57)
}

pub const F_SECP256K1: u64 = (1u64 << 32) + 977;

pub const F_BITLEN: usize = 33;

pub const PAD: usize = 19;

pub const LSBS: usize = 20 + F_BITLEN;

pub const MSBS: usize = PAD;

pub const APPLY_CHUNK: usize = 40;

#[inline]
fn cbit(c: &[u8], i: usize) -> bool {
    let byte = i / 8;
    byte < c.len() && (c[byte] >> (i % 8)) & 1 == 1
}

pub fn cuccaro_carry(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
) {
    let call_index = next_cuccaro_call_index();
    let ops_start = circ.current_ops_len();
    let s = y.len();
    assert_eq!(x.len(), s, "cuccaro_carry: x,y width mismatch");
    let fresh = if cin.is_none() {
        Some(circ.alloc_qubit())
    } else {
        None
    };
    let c: &QubitId = cin.unwrap_or_else(|| fresh.as_ref().unwrap());
    let sum = |circ: &mut B, xi: &QubitId, yi: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *xi, *yi),
        None => circ.cx(*xi, *yi),
    };
    let gated_carry = |circ: &mut B, co: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *c, *co),
        None => circ.cx(*c, *co),
    };
    if s == 0 {
        if let Some(co) = cout {
            gated_carry(circ, co);
        }
    } else {
        for i in 0..s {
            circ.cx(*c, y[i]);
            circ.cx(*c, x[i]);
            if !cuccaro_call_has_structurally_dead_carry(call_index, i) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0200_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(x[i], y[i], *c);
                crate::point_add::restore_op_trace_context(old_context);
            }
        }
        if let Some(co) = cout {
            gated_carry(circ, co);
        }

        for i in (0..s).rev() {
            if !cuccaro_call_has_structurally_dead_carry(call_index, i) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0300_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(x[i], y[i], *c);
                crate::point_add::restore_op_trace_context(old_context);
            }
            circ.cx(*c, y[i]);
            sum(circ, &x[i], &y[i]);
            circ.cx(*c, x[i]);
        }
    }
    if let Some(f) = fresh {
        circ.zero_and_free(f);
    }
    if std::env::var_os("TRACE_TLM_CUCCARO").is_some() {
        eprintln!(
            "TLM_CUCCARO call={} phase={} width={} ctrl={} cin={} cout={} ops_start={} ops_end={}",
            call_index,
            circ.phase,
            s,
            usize::from(ctrl.is_some()),
            usize::from(cin.is_some()),
            usize::from(cout.is_some()),
            ops_start,
            circ.current_ops_len(),
        );
    }
}

fn clean_add_threaded_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
) {
    let s = y.len();
    assert_eq!(x.len(), s, "vented add: x,y width mismatch");

    let gated_sum = |circ: &mut B, xi: &QubitId, yi: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *xi, *yi),
        None => circ.cx(*xi, *yi),
    };
    if s == 0 {
        if let (Some(ci), Some(co)) = (cin, cout) {
            match ctrl {
                Some(ct) => circ.ccx(*ct, *ci, *co),
                None => circ.cx(*ci, *co),
            }
        }
        return;
    }
    let n_inner = if cout.is_some() { s } else { s - 1 };
    let mut inner: Vec<Option<QubitId>> = (0..n_inner).map(|_| Some(circ.alloc_qubit())).collect();
    let produces = |i: usize| cout.is_some() || i + 1 < s;

    for i in 0..s {
        if !produces(i) {
            continue;
        }
        let co = inner[i].as_ref().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, x[i]);
            circ.cx(*ci, y[i]);
            circ.ccx(x[i], y[i], *co);
            circ.cx(*ci, *co);
        } else {
            circ.ccx(x[i], y[i], *co);
        }
    }
    if let Some(cout) = cout {
        let top = inner[s - 1].as_ref().unwrap();
        match ctrl {
            Some(ct) => circ.ccx(*ct, *top, *cout),
            None => circ.cx(*top, *cout),
        }
    }

    for i in (0..s).rev() {
        if !produces(i) {
            let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
            if let Some(ci) = ci {
                circ.cx(*ci, x[i]);
            }
            gated_sum(circ, &x[i], &y[i]);
            if let Some(ci) = ci {
                circ.cx(*ci, x[i]);
            }
            continue;
        }
        let co = inner[i].take().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, co);
        }

        let bit = circ.alloc_bit();
        circ.hmr(co, bit);
        circ.zero_and_free(co);
        circ.cz_if_bit(x[i], y[i], bit);
        if let Some(ci) = ci {
            circ.cx(*ci, y[i]);
        }
        gated_sum(circ, &x[i], &y[i]);
        if let Some(ci) = ci {
            circ.cx(*ci, x[i]);
        }
    }
}

pub(crate) fn erase_carry_gated_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    a: &[QubitId],
    b: &[QubitId],
    cin: &QubitId,
    carry: &QubitId,
    cap: Option<usize>,
) {
    let s = a.len();
    let bit = circ.alloc_bit();
    circ.hmr(*carry, bit);

    circ.loan_zero_qubit(*carry);
    circ.push_condition(bit);
    let deposit = |c: &mut B, ta: &QubitId, tb: &QubitId, c_prev: &QubitId| match ctrl {
        Some(ct) => {
            c.z(*ct);
            c.ccz(*ct, *ta, *tb);
            c.cz(*ct, *c_prev);
        }
        None => {
            c.neg();
            c.cz(*ta, *tb);
            c.z(*c_prev);
        }
    };
    match cap {
        Some(k) if k < s => {
            let lo = s - k;
            let zcin = circ.alloc_qubit();
            super::comparator::compare_geq_cin_middle(circ, &a[lo..], &b[lo..], &zcin, deposit);
            circ.zero_and_free(zcin);
        }
        _ => {
            super::comparator::compare_geq_cin_middle(circ, a, b, cin, deposit);
        }
    }
    circ.pop_condition();
}

pub(crate) fn erase_carry_gated_zero_cin_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    a: &[QubitId],
    b: &[QubitId],
    carry: &QubitId,
    cap: Option<usize>,
) {
    let s = a.len();
    let bit = circ.alloc_bit();
    circ.hmr(*carry, bit);
    circ.push_condition(bit);
    let deposit = |c: &mut B, ta: &QubitId, tb: &QubitId, c_prev: &QubitId| match ctrl {
        Some(ct) => {
            c.z(*ct);
            c.ccz(*ct, *ta, *tb);
            c.cz(*ct, *c_prev);
        }
        None => {
            c.neg();
            c.cz(*ta, *tb);
            c.z(*c_prev);
        }
    };
    match cap {
        Some(k) if k < s => {
            let lo = s - k;
            let zcin = circ.alloc_qubit();
            super::comparator::compare_geq_cin_middle(circ, &a[lo..], &b[lo..], &zcin, deposit);
            circ.zero_and_free(zcin);
        }
        _ => {
            super::comparator::compare_geq_cin_middle(circ, a, b, carry, deposit);
        }
    }
    circ.pop_condition();
}

pub fn controlled_add_vented_chunked_cout(
    circ: &mut B,
    ctrl: &QubitId,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
) {
    add_vented_chunked_opt(circ, Some(ctrl), x, y, chunk, cout, None);
}

pub const CEILING: usize = 1167;

fn emit_chunked_capped(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    bounds: &[(usize, usize)],
    plain_len: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
) {
    let n = y.len();
    let l = n - plain_len;
    let cin0 = circ.alloc_qubit();
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        clean_add_threaded_opt(circ, ctrl, &x[lo..hi], &y[lo..hi], Some(cin), Some(&cy));
        carries.push(cy);
    }
    if l < n {
        let top_cin: &QubitId = carries.last().unwrap_or(&cin0);
        clean_add_threaded_opt(circ, ctrl, &x[l..n], &y[l..n], Some(top_cin), cout);
    } else if let Some(co) = cout {
        circ.cx(*carries.last().unwrap(), *co);
    }
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        erase_carry_gated_opt(circ, ctrl, &y[lo..hi], &x[lo..hi], cin, &carry, cap);
    }
    circ.zero_and_free(cin0);
}

fn hybrid_add_plain(circ: &mut B, a: &[QubitId], b: &[QubitId], vents_budget: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "hybrid_add: a,b width mismatch");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.cx(b[0], a[0]);
        return;
    }
    let vents = vents_budget.min(n - 1);
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(b[i], b[i + 1]);
    }
    let mut vent_ancs: Vec<Option<QubitId>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            let anc = circ.alloc_qubit();
            circ.ccx(a[i], b[i], anc);
            circ.cx(anc, b[i + 1]);
            vent_ancs[i] = Some(anc);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in (0..n - 1).rev() {
        circ.cx(b[i + 1], a[i + 1]);
        if i < vents {
            let anc = vent_ancs[i].take().unwrap();
            circ.cx(anc, b[i + 1]);
            let bit = circ.alloc_bit();
            circ.hmr(anc, bit);
            circ.zero_and_free(anc);
            circ.cz_if_bit(a[i], b[i], bit);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in 1..n - 1 {
        circ.cx(b[i], b[i + 1]);
    }
    circ.cx(b[0], a[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
}

pub(crate) fn hybrid_add_adaptive(circ: &mut B, a: &[QubitId], b: &[QubitId], k: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "adaptive add: a,b width mismatch");
    if n == 0 {
        return;
    }
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    if n <= 4 || k.saturating_add(2 * c) >= n {
        hybrid_add_plain(circ, a, b, k);
        return;
    }
    if k < n.div_ceil(c) + c + super::gidney::ADAPTIVE_RES {
        let cov = (k.saturating_mul(k.saturating_sub(1)) / 2).min(n);
        if cov > 2 * k {
            unreachable!(
                "square adaptive add hit the tight chunked_then_cuccaro branch (n={n}, k={k})"
            );
        }
        hybrid_add_plain(circ, a, b, k);
        return;
    }
    let lay = super::gidney::adaptive_layout(n, k);
    let l = lay.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + lay.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }

    emit_chunked_capped(circ, None, b, a, &bounds, lay.plain_len, None, None);
}

fn add_vented_chunked_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
) {
    add_vented_chunked_opt_capped(circ, ctrl, x, y, chunk, cout, cap, usize::MAX);
}

#[allow(clippy::too_many_arguments)]
fn add_vented_chunked_opt_capped(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
    max_vents: usize,
) {
    let n = y.len();
    assert_eq!(x.len(), n, "chunked add: x,y width mismatch");
    if n == 0 {
        return;
    }

    let c = chunk.clamp(1, n);
    let live = circ.active_qubits as usize;

    let k = CEILING.saturating_sub(live).clamp(1, n).min(max_vents);
    let plain_len = if k >= n {
        n
    } else if c <= 1 {
        0
    } else {
        ((k * c).saturating_sub(n) / (c - 1)).min(n)
    };
    let l = n - plain_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    emit_chunked_capped(circ, ctrl, x, y, &bounds, plain_len, cout, cap);
}

fn ccx_cond(
    circ: &mut B,
    ctrl: &QubitId,
    c1: &QubitId,
    c2: &QubitId,
    t: &QubitId,
    b0: bool,
    b1: bool,
) {
    if b0 {
        circ.cx(*ctrl, *c1);
    }
    if b1 {
        circ.cx(*ctrl, *c2);
    }
    circ.ccx(*c1, *c2, *t);
    if b0 {
        circ.cx(*ctrl, *c1);
    }
    if b1 {
        circ.cx(*ctrl, *c2);
    }
}

fn xor_carries_off_cin(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    off: usize,
    out: &[QubitId],
    cin: &QubitId,
) {
    let n = a.len();
    for i in (1..n - 1).rev() {
        ccx_cond(
            circ,
            ctrl,
            &a[i],
            &out[i - 1],
            &out[i],
            cbit(c, off + i),
            false,
        );
    }
    for i in 0..n - 1 {
        if cbit(c, off + i) {
            circ.cx(*ctrl, out[i]);
        }
    }
    ccx_cond(circ, ctrl, cin, &a[0], &out[0], cbit(c, off), cbit(c, off));
    for i in 1..n - 1 {
        ccx_cond(
            circ,
            ctrl,
            &a[i],
            &out[i - 1],
            &out[i],
            cbit(c, off + i),
            cbit(c, off + i),
        );
    }
}

fn dirty_carryin(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    off: usize,
    dirty: &[QubitId],
    cin: &QubitId,
) {
    let n = a.len();
    debug_assert!(n >= 2 && dirty.len() >= n - 1);
    let mut bits: Vec<BitId> = Vec::with_capacity(n - 1);
    let mut cy_owned: Option<QubitId> = None;
    for i in 0..(n - 1) {
        let new = circ.alloc_qubit();
        let anc = circ.alloc_qubit();
        let on = cbit(c, off + i);
        let cyref: QubitId = match cy_owned {
            Some(q) => q,
            None => *cin,
        };
        if on {
            circ.cx(*ctrl, anc);
        }
        circ.cx(cyref, anc);
        circ.cx(cyref, a[i]);
        circ.ccx(a[i], anc, new);
        circ.cx(cyref, new);
        circ.cx(new, dirty[i]);
        circ.cx(cyref, anc);
        if on {
            circ.cx(*ctrl, anc);
            circ.cx(*ctrl, a[i]);
        }
        circ.zero_and_free(anc);
        if let Some(old) = cy_owned.take() {
            let b = circ.alloc_bit();
            circ.hmr(old, b);
            bits.push(b);
            circ.zero_and_free(old);
        }
        cy_owned = Some(new);
    }
    let cy_top = cy_owned.take().unwrap();
    if cbit(c, off + n - 1) {
        circ.cx(*ctrl, a[n - 1]);
    }
    circ.cx(cy_top, a[n - 1]);
    {
        let b = circ.alloc_bit();
        circ.hmr(cy_top, b);
        bits.push(b);
    }
    circ.zero_and_free(cy_top);
    for i in 0..(n - 1) {
        circ.z_if_bit(dirty[i], bits[i]);
    }
    for q in a {
        circ.x(*q);
    }
    xor_carries_off_cin(circ, ctrl, a, c, off, dirty, cin);
    for q in a {
        circ.x(*q);
    }
    for i in 0..(n - 1) {
        circ.z_if_bit(dirty[i], bits[i]);
    }
}

fn graduated_const_fits(n: usize, k: usize) -> bool {
    k >= 4 && (k - 3) * (k - 2) / 2 >= n
}
fn graduated_const_kmin(n: usize) -> usize {
    (4..).find(|&k| graduated_const_fits(n, k)).unwrap()
}

fn const_chunk_add_clean(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    coff: usize,
    cin: &QubitId,
    cout: &QubitId,
) {
    let call_index = next_const_chunk_call_index();
    let s = a.len();
    if std::env::var_os("TRACE_TLM_CONST_CHUNK").is_some() {
        eprintln!(
            "CONST_CHUNK call={} phase={} width={} coff={} cin={} cout={}",
            call_index, circ.phase, s, coff, cin.0, cout.0,
        );
    }
    if s == 0 {
        return;
    }
    let mut int: Vec<Option<QubitId>> = (0..s - 1).map(|_| Some(circ.alloc_qubit())).collect();
    for i in 0..s {
        let on = cbit(c, coff + i);
        let cin_ref: QubitId = if i == 0 {
            *cin
        } else {
            *int[i - 1].as_ref().unwrap()
        };
        let cout_ref: QubitId = if i == s - 1 {
            *cout
        } else {
            *int[i].as_ref().unwrap()
        };
        circ.cx(cin_ref, a[i]);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        let old_context = crate::point_add::set_op_trace_context(
            0x0800_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        if !const_chunk_call_has_structurally_dead_carry(call_index, i) {
            circ.ccx(a[i], cin_ref, cout_ref);
        }
        crate::point_add::restore_op_trace_context(old_context);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        circ.cx(cin_ref, cout_ref);
    }
    for i in 0..s {
        if cbit(c, coff + i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    for i in (0..s - 1).rev() {
        let on = cbit(c, coff + i);
        let int_i = int[i].take().unwrap();
        let cin_ref: QubitId = if i == 0 {
            *cin
        } else {
            *int[i - 1].as_ref().unwrap()
        };
        if on {
            circ.cx(*ctrl, a[i]);
        }
        circ.cx(cin_ref, int_i);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        let b = circ.alloc_bit();
        circ.hmr(int_i, b);
        circ.zero_and_free(int_i);
        circ.cz_if_bit(a[i], cin_ref, b);
        if on {
            circ.cx(*ctrl, cin_ref);
            circ.cx(*ctrl, a[i]);
        }
    }
}

fn const_chunk_add_clean_drop_cout(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    coff: usize,
    cin: &QubitId,
) {
    let s = a.len();
    if s == 0 {
        return;
    }
    if s == 1 {
        if cbit(c, coff) {
            circ.cx(*ctrl, a[0]);
        }
        circ.cx(*cin, a[0]);
        return;
    }
    let mut int: Vec<Option<QubitId>> = (0..s - 1).map(|_| Some(circ.alloc_qubit())).collect();
    for i in 0..s - 1 {
        let on = cbit(c, coff + i);
        let cin_ref: QubitId = if i == 0 {
            *cin
        } else {
            *int[i - 1].as_ref().unwrap()
        };
        let cout_ref: QubitId = *int[i].as_ref().unwrap();
        circ.cx(cin_ref, a[i]);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        circ.ccx(a[i], cin_ref, cout_ref);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        circ.cx(cin_ref, cout_ref);
    }
    for i in 0..s - 1 {
        if cbit(c, coff + i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    if cbit(c, coff + s - 1) {
        circ.cx(*ctrl, a[s - 1]);
    }
    circ.cx(*int[s - 2].as_ref().unwrap(), a[s - 1]);
    for i in (0..s - 1).rev() {
        let on = cbit(c, coff + i);
        let int_i = int[i].take().unwrap();
        let cin_ref: QubitId = if i == 0 {
            *cin
        } else {
            *int[i - 1].as_ref().unwrap()
        };
        if on {
            circ.cx(*ctrl, a[i]);
        }
        circ.cx(cin_ref, int_i);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        let b = circ.alloc_bit();
        circ.hmr(int_i, b);
        circ.zero_and_free(int_i);
        circ.cz_if_bit(a[i], cin_ref, b);
        if on {
            circ.cx(*ctrl, cin_ref);
            circ.cx(*ctrl, a[i]);
        }
    }
}

fn compare_geq_const_cin_middle<F: FnOnce(&mut B, &QubitId, &QubitId, bool)>(
    circ: &mut B,
    a: &[QubitId],
    c: &[u8],
    coff: usize,
    cin: &QubitId,
    body: F,
) {
    let s = a.len();
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(s);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0);
    cy.push(Some(c0));
    for i in 0..s - 1 {
        let on = cbit(c, coff + i);
        let next = circ.alloc_qubit();
        let ci = *cy[i].as_ref().unwrap();
        circ.ccx(a[i], ci, next);
        if !on {
            circ.cx(a[i], next);
            circ.cx(ci, next);
        }
        cy.push(Some(next));
    }
    {
        let i = s - 1;
        let on = cbit(c, coff + i);
        let ci = *cy[i].as_ref().unwrap();
        body(circ, &a[i], &ci, on);
    }
    for i in (0..s - 1).rev() {
        let on = cbit(c, coff + i);
        let next = cy[i + 1].take().unwrap();
        let ci = *cy[i].as_ref().unwrap();
        if !on {
            circ.cx(ci, next);
            circ.cx(a[i], next);
        }
        let b = circ.alloc_bit();
        circ.hmr(next, b);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, b);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0);
    circ.x(c0);
    circ.zero_and_free(c0);
}

fn controlled_erase_carry_gated_const(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    coff: usize,
    cin: &QubitId,
    carry: QubitId,
) {
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);

    circ.loan_zero_qubit(carry);
    circ.push_condition(bit);
    compare_geq_const_cin_middle(circ, a, c, coff, cin, |cc, a_top, cy_top, ctop| {
        cc.z(*ctrl);
        cc.ccz(*ctrl, *a_top, *cy_top);
        if !ctop {
            cc.cz(*ctrl, *a_top);
            cc.cz(*ctrl, *cy_top);
        }
    });
    circ.pop_condition();
}

fn controlled_add_const_chunked_graduated_off(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[QubitId],
    c: &[u8],
    coff: usize,
    cin: &QubitId,
    k: usize,
) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (mut lo, mut i) = (0usize, 0usize);
    while lo < n && k > i + 3 {
        let cc = (k - 3 - i).min(n - lo);
        bounds.push((lo, lo + cc));
        lo += cc;
        i += 1;
    }
    assert_eq!(lo, n, "graduated staircase (k={k}) covers {lo} < n={n}");
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(clo, chi)) in bounds.iter().enumerate() {
        if std::env::var("TLM_GRAD_FINAL_NO_COUT").ok().as_deref() == Some("1")
            && j + 1 == bounds.len()
        {
            let cin_ref: QubitId = if j == 0 { *cin } else { carries[j - 1] };
            const_chunk_add_clean_drop_cout(circ, ctrl, &a[clo..chi], c, coff + clo, &cin_ref);
            break;
        }
        let cout = circ.alloc_qubit();
        let cin_ref: QubitId = if j == 0 { *cin } else { carries[j - 1] };
        const_chunk_add_clean(circ, ctrl, &a[clo..chi], c, coff + clo, &cin_ref, &cout);
        carries.push(cout);
    }
    for j in (0..carries.len()).rev() {
        let (clo, chi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin_ref: QubitId = if j == 0 { *cin } else { carries[j - 1] };
        controlled_erase_carry_gated_const(
            circ,
            ctrl,
            &a[clo..chi],
            c,
            coff + clo,
            &cin_ref,
            carry,
        );
    }
}

#[allow(clippy::needless_range_loop)]
fn add_f_window_hybrid(
    circ: &mut B,
    ctrl: &QubitId,
    reg: &[QubitId],
    lsbs: usize,
    c: &[u8],
    k: usize,
    trace_call_index: usize,
) {
    let n = lsbs;
    let a: Vec<QubitId> = reg[..n].to_vec();
    let suf_dirty = n - k - 1;
    assert!(
        reg.len() >= lsbs + suf_dirty,
        "+f hybrid: not enough high bits to borrow"
    );
    let dirty: Vec<QubitId> = (lsbs..lsbs + suf_dirty).map(|i| reg[i]).collect();
    let mut cy: Vec<Option<QubitId>> = (0..k).map(|_| Some(circ.alloc_qubit())).collect();

    if cbit(c, 0) {
        circ.ccx(*ctrl, a[0], *cy[0].as_ref().unwrap());
    }
    for i in 1..k {
        let ci = *cy[i - 1].as_ref().unwrap();
        let next = *cy[i].as_ref().unwrap();
        circ.cx(ci, a[i]);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }
        if !ffg_call_has_structurally_dead_hybrid_carry(trace_call_index, i, circ.phase) {
            let old_context = crate::point_add::set_op_trace_context(
                0x0100_0000 | (((trace_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            circ.ccx(a[i], ci, next);
            crate::point_add::restore_op_trace_context(old_context);
        }
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }
        circ.cx(ci, next);
    }

    for i in 0..k {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    let release_cy0_during_suffix = std::env::var("TLM_FFG_RELEASE_CY0_DURING_SUFFIX")
        .ok()
        .as_deref()
        == Some("1")
        && (std::env::var_os("TLM_FFG_RELEASE_CY0_CALLS").is_none()
            || env_index_list_contains("TLM_FFG_RELEASE_CY0_CALLS", trace_call_index))
        && k > 1
        && cbit(c, 0);
    if release_cy0_during_suffix {
        let cy0 = *cy[0].as_ref().unwrap();

        circ.x(a[0]);
        circ.ccx(*ctrl, a[0], cy0);
        circ.x(a[0]);
        circ.loan_zero_qubit(cy0);
    }

    {
        let a_hi: Vec<QubitId> = a[k..].to_vec();
        let cin = *cy[k - 1].as_ref().unwrap();
        let sn = n - k;

        if sn >= 2 {
            controlled_add_const_chunked_graduated_off(
                circ,
                ctrl,
                &a_hi,
                c,
                k,
                &cin,
                graduated_const_kmin(sn),
            );
        } else {
            dirty_carryin(circ, ctrl, &a_hi, c, k, &dirty, &cin);
        }
    }
    if release_cy0_during_suffix {
        let cy0 = *cy[0].as_ref().unwrap();
        circ.reclaim_zero_qubit(cy0);
        circ.x(a[0]);
        circ.ccx(*ctrl, a[0], cy0);
        circ.x(a[0]);
    }

    for i in (1..k).rev() {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]);
        }
        let ci = *cy[i - 1].as_ref().unwrap();
        let next = *cy[i].as_ref().unwrap();
        circ.cx(ci, next);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }
        let nq = cy[i].take().unwrap();
        let b = circ.alloc_bit();
        circ.hmr(nq, b);
        circ.zero_and_free(nq);
        circ.cz_if_bit(a[i], ci, b);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
            circ.cx(*ctrl, a[i]);
        }
    }

    let cy0 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.cx(*ctrl, a[0]);
        let b = circ.alloc_bit();
        circ.hmr(cy0, b);
        circ.zero_and_free(cy0);
        circ.cz_if_bit(a[0], *ctrl, b);
        circ.cx(*ctrl, a[0]);
    } else {
        circ.zero_and_free(cy0);
    }
}

fn add_f_window(
    circ: &mut B,
    ctrl: &QubitId,
    reg: &[QubitId],
    lsbs: usize,
    c: &[u8],
    g_sched: Option<usize>,
) {
    let call_index = next_ffg_call_index();
    let timeline_start = circ.active_timeline.len();
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for +f window");
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.cx(*ctrl, reg[0]);
        }
        return;
    }

    let target_g = super::target_qubit_headroom(circ).map(|headroom| {
        let mut reserve = std::env::var("TLM_TARGET_FFG_RESERVE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(4);
        if let Some(call_reserve) = env_index_value("TLM_TARGET_FFG_CALL_RESERVES", call_index) {
            reserve = call_reserve;
        } else if std::env::var("TLM_TARGET_FFG_RESERVE8_CALLS")
            .ok()
            .map(|value| {
                value
                    .split(',')
                    .filter_map(|item| item.trim().parse::<usize>().ok())
                    .any(|candidate| candidate == call_index)
            })
            .unwrap_or(false)
        {
            reserve = 8;
        }
        if let Some(call_reserve) =
            env_index_value("TLM_TARGET_FFG_CALL_RESERVE_OVERRIDES", call_index)
        {
            reserve = call_reserve;
        }
        headroom.saturating_sub(reserve)
    });
    let scheduled_g = g_sched
        .map_or_else(
            || CEILING.saturating_sub(circ.active_qubits as usize),
            |g| g,
        )
        .min(target_g.unwrap_or(usize::MAX))
        .min(n - 1);
    let capped_g = std::env::var("TLM_FFG_MAX_G")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .map_or(scheduled_g, |cap| scheduled_g.min(cap));
    let g = std::env::var("TLM_FFG_FORCE_G")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .map_or(capped_g, |forced| forced.min(n - 1));
    let trace_entry_active = circ.active_qubits;
    if g >= n - 1 {
        add_f_window_clean(circ, ctrl, reg, lsbs, c);
    } else if g == 0 {
        let cin = circ.alloc_qubit();
        let a_full: Vec<QubitId> = reg[..n].to_vec();
        let dirty: Vec<QubitId> = (lsbs..lsbs + (n - 1)).map(|i| reg[i]).collect();
        dirty_carryin(circ, ctrl, &a_full, c, 0, &dirty, &cin);
        circ.zero_and_free(cin);
    } else {
        add_f_window_hybrid(circ, ctrl, reg, lsbs, c, g, call_index);
    }
    if std::env::var_os("TRACE_TLM_FFG").is_some() {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(trace_entry_active);
        eprintln!(
            "TLM_FFG call={} phase={} g={} entry_active={} local_peak={} phase_max={} ops={}",
            call_index,
            circ.phase,
            g,
            trace_entry_active,
            local_peak,
            circ.current_phase_active_max,
            circ.current_ops_len(),
        );
    }
}

fn add_f_window_clean(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for +f window");
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.cx(*ctrl, reg[0]);
        }
        return;
    }
    let a: Vec<QubitId> = reg[..n].to_vec();

    let mut cy: Vec<Option<QubitId>> = (0..n - 1).map(|_| Some(circ.alloc_qubit())).collect();

    if cbit(c, 0) {
        circ.ccx(*ctrl, a[0], *cy[0].as_ref().unwrap());
    }

    for i in 1..n - 1 {
        let ci = cy[i - 1].take().unwrap();
        let next = cy[i].take().unwrap();
        circ.cx(ci, a[i]);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }
        circ.ccx(a[i], ci, next);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }
        circ.cx(ci, next);
        cy[i - 1] = Some(ci);
        cy[i] = Some(next);
    }

    for i in 0..n - 1 {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    if cbit(c, n - 1) {
        circ.cx(*ctrl, a[n - 1]);
    }
    circ.cx(*cy[n - 2].as_ref().unwrap(), a[n - 1]);

    for i in (1..n - 1).rev() {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
        }

        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, mbit);
        if cbit(c, i) {
            circ.cx(*ctrl, ci);
            circ.cx(*ctrl, a[i]);
        }
        cy[i - 1] = Some(ci);
    }

    let cy1 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.cx(*ctrl, a[0]);
        let mbit = circ.alloc_bit();
        circ.hmr(cy1, mbit);
        circ.zero_and_free(cy1);

        circ.cz_if_bit(a[0], *ctrl, mbit);
        circ.cx(*ctrl, a[0]);
    } else {
        circ.zero_and_free(cy1);
    }
}

fn sub_f_window(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    for q in &reg[..lsbs] {
        circ.x(*q);
    }
    add_f_window(circ, ctrl, reg, lsbs, c, None);
    for q in &reg[..lsbs] {
        circ.x(*q);
    }
}

fn controlled_lt_msbs_conditional(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    a: &[QubitId],
    b: &[QubitId],
    k: usize,
    target: QubitId,
) {
    let a_top: Vec<QubitId> = a[a.len() - k..].to_vec();
    let b_top: Vec<QubitId> = b[b.len() - k..].to_vec();
    let bit = circ.alloc_bit();
    circ.hmr(target, bit);

    circ.zero_and_free(target);
    let ctrl = ctrl.copied();
    circ.push_condition(bit);

    let lt_flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(
        circ,
        &a_top,
        &b_top,
        &lt_flag,
        |c, flag| {
            c.x(*flag);
            match &ctrl {
                Some(ct) => c.cz(*ct, *flag),
                None => c.z(*flag),
            }
            c.x(*flag);
        },
        k,
    );
    circ.zero_and_free(lt_flag);
    circ.pop_condition();
}

fn controlled_add_carry_msbs_conditional(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    a: &[QubitId],
    b: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    let a_top: Vec<QubitId> = a[a.len() - k..].to_vec();
    let b_top: Vec<QubitId> = b[b.len() - k..].to_vec();
    let bit = circ.alloc_bit();
    circ.hmr(*target, bit);
    circ.push_condition(bit);
    for q in &b_top {
        circ.x(*q);
    }

    let ctrl = ctrl.copied();
    let lt_flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(
        circ,
        &b_top,
        &a_top,
        &lt_flag,
        |c, flag| {
            c.x(*flag);
            match &ctrl {
                Some(ct) => c.cz(*ct, *flag),
                None => c.z(*flag),
            }
            c.x(*flag);
        },
        k,
    );
    circ.zero_and_free(lt_flag);
    for q in &b_top {
        circ.x(*q);
    }
    circ.pop_condition();
}

pub fn controlled_mod_add_k(
    circ: &mut B,
    ctrl: &QubitId,
    x: &[QubitId],
    y: &[QubitId],
    sched_k: Option<usize>,
    ffg_g: Option<usize>,
) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 controlled_mod_add expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();

    circ.set_phase("tlm_apply_forward_mod_add_register");
    match sched_k {
        Some(k) => {
            let yr: Vec<&QubitId> = y.iter().collect();
            let xr: Vec<&QubitId> = x.iter().collect();
            super::gidney::controlled_hybrid_add_cout_refs(circ, ctrl, &yr, &xr, &anc, k);
        }
        None => controlled_add_vented_chunked_cout(circ, ctrl, x, y, APPLY_CHUNK, Some(&anc)),
    }

    circ.set_phase("tlm_apply_forward_mod_add_fold");
    add_f_window(circ, &anc, y, LSBS, &f_bytes, ffg_g);

    debug_assert_eq!(MSBS, PAD);

    circ.set_phase("tlm_apply_forward_mod_add_clean");
    controlled_lt_msbs_conditional(circ, Some(ctrl), &y[..n], &x[..n], MSBS, anc);
}

pub fn mod_sub(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_sub expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();

    for q in y {
        circ.x(*q);
    }

    if std::env::var("TLM_SQUARE_NO_VENT_REDUCE").ok().as_deref() == Some("1") {
        cuccaro_carry(circ, None, x, y, None, Some(&anc));
    } else {
        let ci = next_cuccaro_call_index();
        add_cout_vented_skip_dead(circ, x, y, &anc, ci);
    }
    for q in y {
        circ.x(*q);
    }

    sub_f_window(circ, &anc, y, LSBS, &f_bytes);

    controlled_add_carry_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, &anc);
    circ.zero_and_free(anc);
}

fn add_cout_vented_unctrl(circ: &mut B, x: &[QubitId], y: &[QubitId], cout: &QubitId) {
    let n = y.len();
    assert_eq!(x.len(), n, "add_cout_vented_unctrl: x,y width mismatch");
    let zpad = circ.alloc_qubit();
    let mut a: Vec<QubitId> = y.to_vec();
    a.push(*cout);
    let mut b: Vec<QubitId> = x.to_vec();
    b.push(zpad);
    hybrid_add_plain(circ, &a, &b, n);
    circ.zero_and_free(zpad);
}

pub fn mod_rsub_vented_loaded(circ: &mut B, t1: &[QubitId], y: &[QubitId]) {
    let n = y.len();
    assert_eq!(
        t1.len(),
        n,
        "mod_rsub_vented_loaded: t1,y must both be n=256 bits"
    );
    assert_eq!(n, 256, "secp256k1 mod_rsub_vented_loaded expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    for q in y {
        circ.x(*q);
    }
    add_cout_vented_unctrl(circ, t1, y, &anc);
    circ.x(anc);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    circ.x(anc);
    controlled_lt_msbs_conditional(circ, None, &y[..n], &t1[..n], MSBS, anc);
}

fn add_cout_vented_skip_dead(
    circ: &mut B,
    x: &[QubitId],
    y: &[QubitId],
    cout: &QubitId,
    call_index: usize,
) {
    let n = y.len();
    assert_eq!(x.len(), n, "add_cout_vented_skip_dead: x,y width mismatch");
    let dead = |i: usize| cuccaro_call_has_structurally_dead_carry(call_index, i);
    let zpad = circ.alloc_qubit();
    let live = circ.active_qubits as usize;
    let margin = std::env::var("TLM_SQUARE_VENT_MARGIN")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(SQUARE_VENT_MARGIN);
    let vents_budget = square_peak_hard_cap()
        .saturating_sub(live)
        .saturating_sub(margin);

    let mut a: Vec<QubitId> = y.to_vec();
    a.push(*cout);
    let mut b: Vec<QubitId> = x.to_vec();
    b.push(zpad);
    let m = a.len();

    let mut vents_left = vents_budget;
    for i in 1..m {
        circ.cx(b[i], a[i]);
    }
    for i in (1..m - 1).rev() {
        circ.cx(b[i], b[i + 1]);
    }
    let mut vent_ancs: Vec<Option<QubitId>> = (0..m - 1).map(|_| None).collect();
    for i in 0..m - 1 {
        if dead(i) {
            continue;
        }
        if vents_left > 0 {
            let anc = circ.alloc_qubit();
            circ.ccx(a[i], b[i], anc);
            circ.cx(anc, b[i + 1]);
            vent_ancs[i] = Some(anc);
            vents_left -= 1;
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in (0..m - 1).rev() {
        circ.cx(b[i + 1], a[i + 1]);
        if dead(i) {
            continue;
        }
        if let Some(anc) = vent_ancs[i].take() {
            circ.cx(anc, b[i + 1]);
            let bit = circ.alloc_bit();
            circ.hmr(anc, bit);
            circ.zero_and_free(anc);
            circ.cz_if_bit(a[i], b[i], bit);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in 1..m - 1 {
        circ.cx(b[i], b[i + 1]);
    }
    circ.cx(b[0], a[0]);
    for i in 1..m {
        circ.cx(b[i], a[i]);
    }
    circ.zero_and_free(zpad);
}

fn add_cout_vented_unctrl_bounded(circ: &mut B, x: &[QubitId], y: &[QubitId], cout: &QubitId) {
    let n = y.len();
    assert_eq!(
        x.len(),
        n,
        "add_cout_vented_unctrl_bounded: x,y width mismatch"
    );
    let zpad = circ.alloc_qubit();

    let live = circ.active_qubits as usize;
    let margin = std::env::var("TLM_SQUARE_VENT_MARGIN")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(SQUARE_VENT_MARGIN);
    let headroom = square_peak_hard_cap()
        .saturating_sub(live)
        .saturating_sub(margin);
    let mut a: Vec<QubitId> = y.to_vec();
    a.push(*cout);
    let mut b: Vec<QubitId> = x.to_vec();
    b.push(zpad);
    hybrid_add_plain(circ, &a, &b, headroom);
    circ.zero_and_free(zpad);
}

pub const SQUARE_PEAK_HARD_CAP: usize = 1153;

pub const SQUARE_VENT_MARGIN: usize = 30;

fn square_peak_hard_cap() -> usize {
    std::env::var("TLM_SQUARE_PEAK_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(SQUARE_PEAK_HARD_CAP)
}

pub fn mod_add(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_add: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_add expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    add_cout_vented_unctrl(circ, x, y, &anc);

    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));

    controlled_lt_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, anc);
}

pub fn mod_add_exact(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_add_exact: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_add_exact expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    add_cout_vented_unctrl(circ, x, y, &anc);
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));

    controlled_lt_msbs_conditional(circ, None, &y[..n], &x[..n], n, anc);
}

pub fn mod_add_lowpeak(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_add_lowpeak: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_add_lowpeak expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();

    if std::env::var("TLM_SQUARE_NO_VENT_REDUCE").ok().as_deref() == Some("1") {
        cuccaro_carry(circ, None, x, y, None, Some(&anc));
    } else {
        let ci = next_cuccaro_call_index();
        add_cout_vented_skip_dead(circ, x, y, &anc, ci);
    }
    add_f_window(circ, &anc, y, LSBS, &f_bytes, None);
    controlled_lt_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, anc);
}

pub fn mod_add_shifted_low(circ: &mut B, x: &[QubitId], y: &[QubitId], shift: usize) {
    let n = y.len();
    assert_eq!(n, 256, "mod_add_shifted_low expects 256-bit y");
    assert!(shift < n, "shift must be less than 256");
    assert_eq!(x.len(), n - shift, "x must be the low shifted limb");
    if shift == 0 {
        mod_add(circ, x, y);
        return;
    }
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();

    if std::env::var("TLM_SQUARE_VENT_SHIFTED").ok().as_deref() == Some("1") {
        let ci = next_cuccaro_call_index();
        add_cout_vented_skip_dead(circ, x, &y[shift..], &anc, ci);
    } else {
        cuccaro_carry(circ, None, x, &y[shift..], None, Some(&anc));
    }
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));
    controlled_lt_msbs_conditional(circ, None, &y[n - MSBS..], &x[x.len() - MSBS..], MSBS, anc);
}

pub fn mod_sub_vented(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_sub_vented: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_sub_vented expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    for q in y {
        circ.x(*q);
    }
    add_cout_vented_unctrl(circ, x, y, &anc);
    for q in y {
        circ.x(*q);
    }

    for q in &y[..LSBS] {
        circ.x(*q);
    }
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    controlled_add_carry_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, &anc);
    circ.zero_and_free(anc);
}

pub fn mod_sub_shifted_low(circ: &mut B, x: &[QubitId], y: &[QubitId], shift: usize) {
    let n = y.len();
    assert_eq!(n, 256, "mod_sub_shifted_low expects 256-bit y");
    assert!(shift < n, "shift must be less than 256");
    assert_eq!(x.len(), n - shift, "x must be the low shifted limb");
    if shift == 0 {
        mod_sub(circ, x, y);
        return;
    }
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    for q in &y[shift..] {
        circ.x(*q);
    }
    if std::env::var("TLM_SQUARE_VENT_SHIFTED").ok().as_deref() == Some("1") {
        let ci = next_cuccaro_call_index();
        add_cout_vented_skip_dead(circ, x, &y[shift..], &anc, ci);
    } else {
        cuccaro_carry(circ, None, x, &y[shift..], None, Some(&anc));
    }
    for q in &y[shift..] {
        circ.x(*q);
    }
    sub_f_window(circ, &anc, y, LSBS, &f_bytes);
    controlled_add_carry_msbs_conditional(
        circ,
        None,
        &y[n - MSBS..],
        &x[x.len() - MSBS..],
        MSBS,
        &anc,
    );
    circ.zero_and_free(anc);
}

fn toggle_pattern_mcx(circ: &mut B, pattern: &[(QubitId, bool)], target: &QubitId) {
    for &(q, expected) in pattern {
        if !expected {
            circ.x(q);
        }
    }
    let ctrls: Vec<&QubitId> = pattern.iter().map(|(q, _)| q).collect();
    super::mcx::mcx_clean_k(circ, &ctrls, target);
    for &(q, expected) in pattern.iter().rev() {
        if !expected {
            circ.x(q);
        }
    }
}

fn toggle_geq_small_const(circ: &mut B, a: &[QubitId], threshold: usize, target: &QubitId) {
    assert!(threshold < (1usize << a.len()));
    for j in (0..a.len()).rev() {
        if (threshold >> j) & 1 != 0 {
            continue;
        }
        let mut pattern = Vec::with_capacity(a.len() - j);
        for k in (j + 1)..a.len() {
            pattern.push((a[k], (threshold >> k) & 1 != 0));
        }
        pattern.push((a[j], true));
        toggle_pattern_mcx(circ, &pattern, target);
    }
    let equality: Vec<(QubitId, bool)> = a
        .iter()
        .enumerate()
        .map(|(i, &q)| (q, (threshold >> i) & 1 != 0))
        .collect();
    toggle_pattern_mcx(circ, &equality, target);
}

fn toggle_geq_p_minus_low3(circ: &mut B, y: &[QubitId], c: &[QubitId], target: &QubitId) {
    debug_assert_eq!(y.len(), 256);
    debug_assert_eq!(c.len(), 3);

    let sum: Vec<QubitId> = (0..11).map(|_| circ.alloc_qubit()).collect();
    for i in 0..10 {
        circ.cx(y[i], sum[i]);
    }
    let zeros: Vec<QubitId> = (0..8).map(|_| circ.alloc_qubit()).collect();
    let mut c11 = c.to_vec();
    c11.extend(zeros.iter().copied());
    cuccaro_carry(circ, None, &c11, &sum, None, None);

    let low_ge = circ.alloc_qubit();
    toggle_geq_small_const(circ, &sum, 47, &low_ge);
    let lower = circ.alloc_qubit();
    circ.cx(y[32], lower);
    let mut lower_pattern = Vec::with_capacity(24);
    lower_pattern.push((y[32], false));
    lower_pattern.extend(y[10..32].iter().map(|&q| (q, true)));
    lower_pattern.push((low_ge, true));
    toggle_pattern_mcx(circ, &lower_pattern, &lower);

    let mut full_pattern = Vec::with_capacity(224);
    full_pattern.push((lower, true));
    full_pattern.extend(y[33..].iter().map(|&q| (q, true)));
    toggle_pattern_mcx(circ, &full_pattern, target);

    toggle_pattern_mcx(circ, &lower_pattern, &lower);
    circ.cx(y[32], lower);
    circ.zero_and_free(lower);
    toggle_geq_small_const(circ, &sum, 47, &low_ge);
    circ.zero_and_free(low_ge);

    for q in &sum {
        circ.x(*q);
    }
    cuccaro_carry(circ, None, &c11, &sum, None, None);
    for q in &sum {
        circ.x(*q);
    }
    for i in 0..10 {
        circ.cx(y[i], sum[i]);
    }
    for q in sum {
        circ.zero_and_free(q);
    }
    for q in zeros {
        circ.zero_and_free(q);
    }
}

pub fn mod_sub_classical_low3(circ: &mut B, y: &[QubitId], c: &[BitId]) {
    assert_eq!(y.len(), 256, "mod_sub_classical_low3 expects 256-bit y");
    assert_eq!(
        c.len(),
        3,
        "mod_sub_classical_low3 expects three classical bits"
    );

    let cq: Vec<QubitId> = (0..3).map(|_| circ.alloc_qubit()).collect();
    for i in 0..3 {
        circ.x_if_bit(cq[i], c[i]);
    }

    let low_borrow = circ.alloc_qubit();
    for q in &y[..3] {
        circ.x(*q);
    }
    cuccaro_carry(circ, None, &cq, &y[..3], None, Some(&low_borrow));
    for q in &y[..3] {
        circ.x(*q);
    }

    let full_borrow = circ.alloc_qubit();
    let mut borrow_pattern = Vec::with_capacity(254);
    borrow_pattern.push((low_borrow, true));
    borrow_pattern.extend(y[3..].iter().map(|&q| (q, false)));
    toggle_pattern_mcx(circ, &borrow_pattern, &full_borrow);

    for q in &y[3..] {
        circ.x(*q);
    }
    super::mcx::cinc_khattar_gidney(circ, &y[3..], &low_borrow);
    for q in &y[3..] {
        circ.x(*q);
    }

    let low_copy: Vec<QubitId> = (0..3).map(|_| circ.alloc_qubit()).collect();
    for i in 0..3 {
        circ.cx(y[i], low_copy[i]);
    }
    cuccaro_carry(circ, None, &cq, &low_copy, None, Some(&low_borrow));
    for q in &low_copy {
        circ.x(*q);
    }
    cuccaro_carry(circ, None, &cq, &low_copy, None, None);
    for q in &low_copy {
        circ.x(*q);
    }
    for i in 0..3 {
        circ.cx(y[i], low_copy[i]);
    }
    for q in low_copy {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(low_borrow);

    let f_bytes = F_SECP256K1.to_le_bytes();
    sub_f_window(circ, &full_borrow, y, LSBS, &f_bytes);
    toggle_geq_p_minus_low3(circ, y, &cq, &full_borrow);
    circ.zero_and_free(full_borrow);

    for i in 0..3 {
        circ.x_if_bit(cq[i], c[i]);
    }
    for q in cq {
        circ.zero_and_free(q);
    }
}

pub fn mod_neg(circ: &mut B, x: &[QubitId]) {
    let n = x.len();
    assert_eq!(n, 256, "secp256k1 mod_neg expects n=256");
    let f_minus_1 = (F_SECP256K1 - 1).to_le_bytes();
    add_const_window_clean(circ, x, n, &f_minus_1);
    for q in x {
        circ.x(*q);
    }
}

fn add_const_window_clean(circ: &mut B, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    let add_const_call_index = next_add_const_call_index();
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for const window");
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.x(reg[0]);
        }
        return;
    }
    let a: Vec<QubitId> = reg[..n].to_vec();
    let mut cy: Vec<Option<QubitId>> = (0..n - 1).map(|_| Some(circ.alloc_qubit())).collect();

    if cbit(c, 0) {
        circ.cx(a[0], *cy[0].as_ref().unwrap());
    }

    for i in 1..n - 1 {
        let ci = cy[i - 1].take().unwrap();
        let next = cy[i].take().unwrap();
        circ.cx(ci, a[i]);
        if cbit(c, i) {
            circ.x(ci);
        }
        if !add_const_has_structurally_dead_carry(add_const_call_index, i) {
            let old_context = crate::point_add::set_op_trace_context(
                0x1100_0000 | (((add_const_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            circ.ccx(a[i], ci, next);
            crate::point_add::restore_op_trace_context(old_context);
        }
        if cbit(c, i) {
            circ.x(ci);
        }
        circ.cx(ci, next);
        cy[i - 1] = Some(ci);
        cy[i] = Some(next);
    }

    for i in 0..n - 1 {
        if cbit(c, i) {
            circ.x(a[i]);
        }
    }
    if cbit(c, n - 1) {
        circ.x(a[n - 1]);
    }
    circ.cx(*cy[n - 2].as_ref().unwrap(), a[n - 1]);

    for i in (1..n - 1).rev() {
        if cbit(c, i) {
            circ.x(a[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if cbit(c, i) {
            circ.x(ci);
        }
        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, mbit);
        if cbit(c, i) {
            circ.x(ci);
            circ.x(a[i]);
        }
        cy[i - 1] = Some(ci);
    }

    let cy1 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.x(a[0]);
        let mbit = circ.alloc_bit();
        circ.hmr(cy1, mbit);
        circ.zero_and_free(cy1);
        circ.z_if_bit(a[0], mbit);
        circ.x(a[0]);
    } else {
        circ.zero_and_free(cy1);
    }
}

pub fn mod_double(circ: &mut B, a: &[QubitId]) {
    let n = a.len() - 1;
    assert_eq!(n, 256, "secp256k1 mod_double expects 257-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();

    for i in (0..n).rev() {
        circ.swap(a[i], a[i + 1]);
    }

    add_f_window(circ, &a[n], a, LSBS, &f_bytes, None);

    circ.cx(a[0], a[n]);
}

pub fn mod_double_reverse(circ: &mut B, a: &[QubitId]) {
    let n = a.len() - 1;
    assert_eq!(n, 256, "secp256k1 mod_double_reverse expects 257-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();

    circ.cx(a[0], a[n]);

    sub_f_window(circ, &a[n], a, LSBS, &f_bytes);

    for i in 0..n {
        circ.swap(a[i], a[i + 1]);
    }
}

pub fn add_f_window_pub(
    circ: &mut B,
    ctrl: &QubitId,
    reg: &[QubitId],
    lsbs: usize,
    c: &[u8],
    g_sched: Option<usize>,
) {
    add_f_window(circ, ctrl, reg, lsbs, c, g_sched);
}
