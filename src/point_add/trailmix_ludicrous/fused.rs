
use super::arith::{F_SECP256K1, LSBS};
use super::{B, BExt};
use crate::circuit::{BitId, QubitId};
use std::cell::Cell;

thread_local! {
    static FOLD_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static ACTIVE_FOLD_CALL_INDEX: Cell<usize> = const { Cell::new(usize::MAX) };
    static FOLD_CHUNK_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FOLD_DIRTY_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FOLD_CLEAN_WINDOW_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FOLD_BOUNDARY_ZERO_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FUSED_CDOUBLE_FWD_SHIFT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static FUSED_CDOUBLE_REV_SHIFT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_fold_call_index() {
    FOLD_CALL_INDEX.with(|index| index.set(0));
    ACTIVE_FOLD_CALL_INDEX.with(|index| index.set(usize::MAX));
    FOLD_CHUNK_CALL_INDEX.with(|index| index.set(0));
    FOLD_DIRTY_CALL_INDEX.with(|index| index.set(0));
    FOLD_CLEAN_WINDOW_CALL_INDEX.with(|index| index.set(0));
    FOLD_BOUNDARY_ZERO_CALL_INDEX.with(|index| index.set(0));
    FUSED_CDOUBLE_FWD_SHIFT_CALL_INDEX.with(|index| index.set(0));
    FUSED_CDOUBLE_REV_SHIFT_CALL_INDEX.with(|index| index.set(0));
}

fn next_fold_call_index() -> usize {
    FOLD_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn active_fold_call_index() -> usize {
    ACTIVE_FOLD_CALL_INDEX.with(|index| index.get())
}

fn enter_fold_call_index(index: usize) -> usize {
    ACTIVE_FOLD_CALL_INDEX.with(|slot| {
        let prior = slot.get();
        slot.set(index);
        prior
    })
}

fn restore_fold_call_index(prior: usize) {
    ACTIVE_FOLD_CALL_INDEX.with(|slot| slot.set(prior));
}

fn next_fold_chunk_call_index() -> usize {
    FOLD_CHUNK_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_fold_dirty_call_index() -> usize {
    FOLD_DIRTY_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_fold_clean_window_call_index() -> usize {
    FOLD_CLEAN_WINDOW_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_fold_boundary_zero_call_index() -> usize {
    FOLD_BOUNDARY_ZERO_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_fused_cdouble_fwd_shift_call_index() -> usize {
    FUSED_CDOUBLE_FWD_SHIFT_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_fused_cdouble_rev_shift_call_index() -> usize {
    FUSED_CDOUBLE_REV_SHIFT_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn env_index_value(name: &str, index: usize) -> Option<usize> {
    std::env::var(name)
        .ok()
        .and_then(|value| {
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

fn skip_structural_dead_fused_carries() -> bool {
    std::env::var_os("TLM_FUSED_SKIP_STRUCTURAL_DEAD_CARRIES").is_some()
}

fn skip_structural_dead_fused_cdouble_shift0() -> bool {
    std::env::var_os("TLM_FUSED_SKIP_STRUCTURAL_DEAD_SHIFT0").is_some()
}

fn skip_structural_dead_fused_dirty_fold() -> bool {
    skip_structural_dead_fused_carries()
        && std::env::var_os("TLM_FUSED_SKIP_STRUCTURAL_DEAD_DIRTY_FOLD").is_some()
}

fn skip_structural_dead_fused_clean_window() -> bool {
    skip_structural_dead_fused_carries()
        && std::env::var_os("TLM_FUSED_SKIP_STRUCTURAL_DEAD_CLEAN_WINDOW").is_some()
}

fn skip_exact_fused_clean_fold() -> bool {
    skip_structural_dead_fused_carries()
        && (std::env::var_os("TLM_FUSED_SKIP_EXACT_FOLD_REMAINDER").is_some()
            || std::env::var_os("TLM_FUSED_SKIP_EXACT_CLEAN_FOLD").is_some())
}

fn skip_exact_fused_chunk_fold() -> bool {
    skip_structural_dead_fused_carries()
        && (std::env::var_os("TLM_FUSED_SKIP_EXACT_FOLD_REMAINDER").is_some()
            || std::env::var_os("TLM_FUSED_SKIP_EXACT_CHUNK_FOLD").is_some())
}

fn skip_fused_clean_fold_top31() -> bool {
    skip_structural_dead_fused_carries()
        && std::env::var_os("TLM_FUSED_CLEAN_FOLD_SKIP_TOP31").is_some()
}

const FUSED_CLEAN_FOLD_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (345, 0, 3),
    (345, 29, 31),
    (327, 0, 3),
    (327, 30, 31),
    (353, 0, 3),
    (353, 30, 31),
    (355, 0, 3),
    (355, 30, 31),
    (377, 0, 3),
    (377, 30, 31),
    (419, 0, 3),
    (419, 29, 30),
    (431, 0, 3),
    (431, 30, 31),
    (432, 0, 3),
    (432, 29, 29),
    (432, 31, 31),
    (438, 0, 3),
    (438, 30, 31),
    (440, 0, 3),
    (440, 30, 31),
    (446, 0, 3),
    (446, 30, 31),
    (456, 0, 3),
    (456, 30, 31),
    (484, 0, 3),
    (484, 30, 31),
    (491, 0, 3),
    (491, 30, 31),
    (511, 0, 3),
    (511, 30, 31),
    (0, 1, 3),
    (0, 29, 30),
    (101, 1, 3),
    (101, 30, 31),
    (12, 1, 3),
    (12, 30, 31),
    (121, 1, 3),
    (121, 29, 29),
    (121, 31, 31),
    (148, 1, 3),
    (148, 30, 31),
    (156, 1, 3),
    (156, 30, 31),
    (167, 1, 3),
    (167, 29, 29),
    (167, 31, 31),
    (191, 1, 3),
    (191, 30, 31),
    (192, 1, 3),
    (192, 30, 31),
    (318, 0, 3),
    (318, 31, 31),
    (320, 0, 3),
    (320, 31, 31),
    (321, 0, 3),
    (321, 31, 31),
    (323, 0, 3),
    (323, 31, 31),
    (324, 0, 3),
    (324, 31, 31),
    (325, 0, 3),
    (325, 31, 31),
    (329, 0, 3),
    (329, 30, 30),
    (333, 0, 3),
    (333, 31, 31),
    (337, 0, 3),
    (337, 30, 30),
    (338, 0, 3),
    (338, 31, 31),
    (340, 0, 3),
    (340, 31, 31),
    (343, 0, 3),
    (343, 31, 31),
    (344, 0, 3),
    (344, 31, 31),
    (346, 0, 3),
    (346, 31, 31),
    (347, 0, 3),
    (347, 30, 30),
    (348, 0, 3),
    (348, 31, 31),
    (349, 0, 3),
    (349, 31, 31),
    (350, 0, 3),
    (350, 31, 31),
    (354, 0, 3),
    (354, 30, 30),
    (360, 0, 3),
    (360, 31, 31),
    (362, 0, 3),
    (362, 31, 31),
    (363, 0, 3),
    (363, 31, 31),
    (365, 0, 3),
    (365, 31, 31),
    (366, 0, 3),
    (366, 31, 31),
    (368, 0, 3),
    (368, 31, 31),
    (370, 0, 3),
    (370, 31, 31),
    (371, 0, 3),
    (371, 30, 30),
    (372, 0, 3),
    (372, 31, 31),
    (373, 0, 3),
    (373, 31, 31),
    (376, 0, 3),
    (376, 31, 31),
    (378, 0, 3),
    (378, 31, 31),
    (379, 0, 3),
    (379, 31, 31),
    (380, 0, 3),
    (380, 31, 31),
    (381, 0, 3),
    (381, 31, 31),
    (385, 0, 3),
    (385, 31, 31),
    (387, 0, 3),
    (387, 31, 31),
    (393, 0, 3),
    (393, 30, 30),
    (394, 0, 3),
    (394, 31, 31),
    (399, 0, 3),
    (399, 31, 31),
    (400, 0, 3),
    (400, 31, 31),
    (402, 0, 3),
    (402, 31, 31),
    (403, 0, 3),
    (403, 31, 31),
    (406, 0, 3),
    (406, 30, 30),
    (410, 0, 3),
    (410, 31, 31),
    (411, 0, 3),
    (411, 31, 31),
    (414, 0, 3),
    (414, 31, 31),
    (415, 0, 3),
    (415, 31, 31),
    (420, 0, 3),
    (420, 31, 31),
    (425, 0, 3),
    (425, 31, 31),
    (429, 0, 3),
    (429, 31, 31),
    (434, 0, 3),
    (434, 31, 31),
    (435, 0, 3),
    (435, 31, 31),
    (436, 0, 3),
    (436, 31, 31),
    (437, 0, 3),
    (437, 31, 31),
    (447, 0, 3),
    (447, 31, 31),
    (449, 0, 3),
    (449, 31, 31),
    (452, 0, 3),
    (452, 31, 31),
    (454, 0, 3),
    (454, 31, 31),
    (455, 0, 3),
    (455, 30, 30),
    (457, 0, 3),
    (457, 31, 31),
    (461, 0, 3),
    (461, 31, 31),
    (463, 0, 3),
    (463, 31, 31),
    (465, 0, 3),
    (465, 31, 31),
    (466, 0, 3),
    (466, 31, 31),
    (475, 0, 3),
    (475, 30, 30),
    (479, 0, 3),
    (479, 31, 31),
    (480, 0, 3),
    (480, 31, 31),
    (485, 0, 3),
    (485, 31, 31),
    (486, 0, 3),
    (486, 30, 30),
    (493, 0, 3),
    (493, 31, 31),
    (496, 0, 3),
    (496, 31, 31),
    (497, 0, 3),
    (497, 31, 31),
    (500, 0, 3),
    (500, 31, 31),
    (501, 0, 3),
    (501, 31, 31),
    (502, 0, 3),
    (502, 31, 31),
    (503, 0, 3),
    (503, 30, 30),
    (507, 0, 3),
    (507, 31, 31),
    (508, 0, 3),
    (508, 30, 30),
    (509, 0, 3),
    (509, 31, 31),
    (51, 1, 3),
    (51, 30, 31),
    (510, 0, 3),
    (510, 31, 31),
    (513, 0, 3),
    (513, 31, 31),
    (54, 1, 3),
    (54, 30, 31),
    (64, 1, 3),
    (64, 30, 31),
    (10, 1, 3),
    (10, 31, 31),
    (100, 1, 3),
    (100, 31, 31),
    (104, 1, 3),
    (104, 30, 30),
    (105, 1, 3),
    (105, 31, 31),
    (107, 1, 3),
    (107, 31, 31),
    (108, 1, 3),
    (108, 31, 31),
    (111, 1, 3),
    (111, 31, 31),
    (112, 1, 3),
    (112, 31, 31),
    (115, 1, 3),
    (115, 31, 31),
    (116, 1, 3),
    (116, 31, 31),
    (117, 1, 3),
    (117, 31, 31),
    (124, 1, 3),
    (124, 31, 31),
    (125, 1, 3),
    (125, 30, 30),
    (126, 1, 3),
    (126, 31, 31),
    (127, 1, 3),
    (127, 30, 30),
    (128, 1, 3),
    (128, 31, 31),
    (13, 1, 3),
    (13, 30, 30),
    (131, 1, 3),
    (131, 31, 31),
    (132, 1, 3),
    (132, 31, 31),
    (134, 1, 3),
    (134, 31, 31),
    (135, 1, 3),
    (135, 31, 31),
    (137, 1, 3),
    (137, 31, 31),
    (14, 1, 3),
    (14, 31, 31),
    (142, 1, 3),
    (142, 31, 31),
    (144, 1, 3),
    (144, 31, 31),
    (15, 1, 3),
    (15, 31, 31),
    (154, 1, 3),
    (154, 30, 30),
    (157, 1, 3),
    (157, 31, 31),
    (158, 1, 3),
    (158, 31, 31),
    (160, 1, 3),
    (160, 31, 31),
    (161, 1, 3),
    (161, 31, 31),
    (162, 1, 3),
    (162, 31, 31),
    (164, 1, 3),
    (164, 31, 31),
    (166, 1, 3),
    (166, 31, 31),
    (169, 1, 3),
    (169, 31, 31),
    (171, 1, 3),
    (171, 31, 31),
    (172, 1, 3),
    (172, 31, 31),
    (174, 1, 3),
    (174, 31, 31),
    (175, 1, 3),
    (175, 31, 31),
    (178, 1, 3),
    (178, 30, 30),
    (179, 1, 3),
    (179, 31, 31),
    (18, 1, 3),
    (18, 31, 31),
    (180, 1, 3),
    (180, 31, 31),
    (181, 1, 3),
    (181, 31, 31),
    (182, 1, 3),
    (182, 31, 31),
    (189, 1, 3),
    (189, 31, 31),
    (2, 1, 3),
    (2, 31, 31),
    (23, 1, 3),
    (23, 31, 31),
    (24, 1, 3),
    (24, 31, 31),
    (26, 1, 3),
    (26, 31, 31),
    (29, 1, 3),
    (29, 31, 31),
    (31, 1, 3),
    (31, 31, 31),
    (322, 0, 3),
    (326, 0, 3),
    (328, 0, 3),
    (330, 0, 3),
    (331, 0, 3),
    (332, 0, 3),
    (334, 0, 3),
    (335, 0, 3),
    (336, 0, 3),
    (339, 0, 3),
    (341, 0, 3),
    (342, 0, 3),
    (351, 0, 3),
    (352, 0, 3),
    (356, 0, 3),
    (357, 0, 3),
    (358, 0, 3),
    (359, 0, 3),
    (361, 0, 3),
    (364, 0, 3),
    (367, 0, 3),
    (369, 0, 3),
    (37, 1, 3),
    (37, 31, 31),
    (374, 0, 3),
    (375, 0, 3),
    (38, 1, 3),
    (38, 31, 31),
    (382, 0, 3),
    (383, 0, 3),
    (384, 0, 3),
    (386, 0, 3),
    (388, 0, 3),
    (389, 0, 3),
    (390, 0, 3),
    (391, 0, 3),
    (392, 0, 3),
    (395, 0, 3),
    (396, 0, 3),
    (397, 0, 3),
    (398, 0, 3),
    (4, 1, 3),
    (4, 31, 31),
    (40, 1, 3),
    (40, 31, 31),
    (401, 0, 3),
    (404, 0, 3),
    (405, 0, 3),
    (407, 0, 3),
    (408, 0, 3),
    (409, 0, 3),
    (412, 0, 3),
    (413, 0, 3),
    (416, 0, 3),
    (417, 0, 3),
    (418, 0, 3),
    (421, 0, 3),
    (422, 0, 3),
    (423, 0, 3),
    (424, 0, 3),
    (426, 0, 3),
    (427, 0, 3),
    (428, 0, 3),
    (43, 1, 3),
    (43, 31, 31),
    (430, 0, 3),
    (433, 0, 3),
    (439, 0, 3),
    (44, 1, 3),
    (44, 31, 31),
    (441, 0, 3),
    (442, 0, 3),
    (443, 0, 3),
    (444, 0, 3),
    (445, 0, 3),
    (448, 0, 3),
    (45, 1, 3),
    (45, 31, 31),
    (450, 0, 3),
    (451, 0, 3),
    (453, 0, 3),
    (458, 0, 3),
    (459, 0, 3),
    (46, 1, 3),
    (46, 31, 31),
    (460, 0, 3),
    (462, 0, 3),
    (464, 0, 3),
    (467, 0, 3),
    (468, 0, 3),
    (469, 0, 3),
    (470, 0, 3),
    (471, 0, 3),
    (472, 0, 3),
    (473, 0, 3),
    (474, 0, 3),
    (476, 0, 3),
    (477, 0, 3),
    (478, 0, 3),
    (481, 0, 3),
    (482, 0, 3),
    (483, 0, 3),
    (487, 0, 3),
    (488, 0, 3),
    (489, 0, 3),
    (490, 0, 3),
    (492, 0, 3),
    (494, 0, 3),
    (495, 0, 3),
    (498, 0, 3),
    (499, 0, 3),
    (504, 0, 3),
    (505, 0, 3),
    (506, 0, 3),
    (512, 0, 3),
    (53, 1, 3),
    (53, 31, 31),
    (56, 1, 3),
    (56, 31, 31),
    (57, 1, 3),
    (57, 31, 31),
    (59, 1, 3),
    (59, 31, 31),
    (6, 1, 3),
    (6, 31, 31),
    (61, 1, 3),
    (61, 31, 31),
    (65, 1, 3),
    (65, 31, 31),
    (67, 1, 3),
    (67, 31, 31),
    (7, 1, 3),
    (7, 31, 31),
    (71, 1, 3),
    (71, 31, 31),
    (74, 1, 3),
    (74, 31, 31),
    (75, 1, 3),
    (75, 31, 31),
    (78, 1, 3),
    (78, 31, 31),
    (79, 1, 3),
    (79, 31, 31),
    (85, 1, 3),
    (85, 30, 30),
    (86, 1, 3),
    (86, 31, 31),
    (87, 1, 3),
    (87, 31, 31),
    (88, 1, 3),
    (88, 31, 31),
    (89, 1, 3),
    (89, 31, 31),
    (9, 1, 3),
    (9, 31, 31),
    (90, 1, 3),
    (90, 31, 31),
    (91, 1, 3),
    (91, 30, 30),
    (94, 1, 3),
    (94, 29, 29),
    (95, 1, 3),
    (95, 31, 31),
    (96, 1, 3),
    (96, 31, 31),
    (98, 1, 3),
    (98, 31, 31),
    (99, 1, 3),
    (99, 31, 31),
];

const FUSED_CHUNK_FOLD_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (1008, 0, 3),
    (1022, 0, 3),
    (1036, 0, 3),
    (1050, 0, 3),
    (1064, 0, 3),
    (1078, 0, 3),
    (1092, 0, 3),
    (1106, 0, 3),
    (112, 0, 3),
    (1120, 0, 3),
    (1134, 0, 3),
    (1148, 0, 3),
    (1162, 0, 3),
    (1176, 0, 3),
    (1190, 0, 3),
    (1204, 0, 3),
    (1218, 0, 3),
    (1232, 0, 3),
    (1246, 0, 3),
    (126, 0, 3),
    (1260, 0, 3),
    (1274, 0, 3),
    (1288, 0, 3),
    (1302, 0, 3),
    (1316, 0, 3),
    (1330, 0, 3),
    (1344, 0, 3),
    (1358, 0, 3),
    (1372, 0, 3),
    (1386, 0, 3),
    (140, 0, 3),
    (1400, 0, 3),
    (1414, 0, 3),
    (1428, 0, 3),
    (1442, 0, 3),
    (1456, 0, 3),
    (1470, 0, 3),
    (1484, 0, 3),
    (1498, 0, 3),
    (1512, 0, 3),
    (1526, 0, 3),
    (154, 0, 3),
    (1540, 0, 3),
    (1554, 0, 3),
    (168, 0, 3),
    (182, 0, 3),
    (196, 0, 3),
    (210, 0, 3),
    (224, 0, 3),
    (238, 0, 3),
    (252, 0, 3),
    (266, 0, 3),
    (280, 0, 3),
    (294, 0, 3),
    (308, 0, 3),
    (322, 0, 3),
    (336, 0, 3),
    (350, 0, 3),
    (364, 0, 3),
    (378, 0, 3),
    (392, 0, 3),
    (406, 0, 3),
    (420, 0, 3),
    (434, 0, 3),
    (448, 0, 3),
    (462, 0, 3),
    (476, 0, 3),
    (490, 0, 3),
    (504, 0, 3),
    (518, 0, 3),
    (532, 0, 3),
    (546, 0, 3),
    (560, 0, 3),
    (574, 0, 3),
    (763, 0, 3),
    (777, 0, 3),
    (784, 0, 3),
    (791, 0, 3),
    (796, 0, 3),
    (798, 0, 3),
    (805, 0, 3),
    (812, 0, 3),
    (819, 0, 3),
    (826, 0, 3),
    (833, 0, 3),
    (84, 0, 3),
    (840, 0, 3),
    (847, 0, 3),
    (854, 0, 3),
    (868, 0, 3),
    (882, 0, 3),
    (896, 0, 3),
    (910, 0, 3),
    (924, 0, 3),
    (938, 0, 3),
    (952, 0, 3),
    (966, 0, 3),
    (98, 0, 3),
    (980, 0, 3),
    (994, 0, 3),
];

const FUSED_DIRTY_FOLD_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (6, 0, 3),
    (6, 21, 24),
    (6, 26, 31),
    (6, 43, 45),
    (6, 47, 51),
    (5, 7, 18),
    (5, 32, 38),
    (4, 9, 18),
    (4, 30, 30),
    (4, 32, 38),
    (3, 15, 23),
    (3, 36, 42),
    (1, 11, 19),
    (1, 34, 39),
    (7, 10, 18),
    (7, 32, 32),
    (7, 34, 38),
    (2, 10, 18),
    (2, 33, 37),
    (8, 11, 19),
    (8, 35, 36),
    (8, 38, 39),
    (0, 11, 17),
    (0, 34, 34),
    (0, 36, 38),
    (9, 12, 18),
    (9, 35, 38),
];

const FUSED_CLEAN_WINDOW_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (0, 1, 3),
    (1, 1, 3),
    (2, 1, 3),
    (3, 1, 3),
    (4, 1, 3),
    (5, 0, 3),
    (6, 0, 3),
    (7, 0, 3),
    (8, 0, 3),
];

const FUSED_CLEAN_FOLD_REMAINDER_KEYS: &[u32] = &[
    257, 258, 259, 769, 770, 771, 1281, 1282, 1283, 2049, 2050, 2051, 2817, 2818,
    2819, 4097, 4098, 4099, 4353, 4354, 4355, 4865, 4866, 4867, 5121, 5122, 5123,
    5377, 5378, 5379, 5633, 5634, 5635, 6401, 6402, 6403, 6913, 6914, 6915, 7169,
    7170, 7171, 7681, 7682, 7683, 8193, 8194, 8195, 8449, 8450, 8451, 8705, 8706,
    8707, 8961, 8962, 8963, 9217, 9218, 9219, 9985, 9986, 9987, 10497, 10498,
    10499, 10753, 10754, 10755, 12033, 12034, 12035, 12289, 12290, 12291, 12545,
    12546, 12547, 12801, 12802, 12803, 13313, 13314, 13315, 14081, 14082, 14083,
    14849, 14850, 14851, 15361, 15362, 15363, 15873, 15874, 15875, 16129, 16130,
    16131, 16897, 16898, 16899, 17409, 17410, 17411, 17665, 17666, 17667, 17921,
    17922, 17923, 18433, 18434, 18435, 18689, 18690, 18691, 19457, 19458, 19459,
    19713, 19714, 19715, 20481, 20482, 20483, 20737, 20738, 20739, 20993, 20994,
    20995, 21249, 21250, 21251, 21505, 21506, 21507, 23553, 23554, 23555, 23809,
    23810, 23811, 24833, 24834, 24835, 26113, 26114, 26115, 26369, 26370, 26371,
    27137, 27138, 27139, 27905, 27906, 27907, 28161, 28162, 28163, 28929, 28930,
    28931, 29185, 29186, 29187, 30209, 30210, 30211, 30465, 30466, 30467, 30721,
    30722, 30723, 31233, 31234, 31235, 31489, 31490, 31491, 33025, 33026, 33027,
    33281, 33282, 33283, 34049, 34050, 34051, 34817, 34818, 34819, 35329, 35330,
    35331, 35585, 35586, 35587, 35841, 35842, 35843, 36097, 36098, 36099, 36609,
    36610, 36611, 37121, 37122, 37123, 37377, 37378, 37379, 37633, 37634, 37635,
    38145, 38146, 38147, 38401, 38402, 38403, 38657, 38658, 38659, 38913, 38914,
    38915, 39169, 39170, 39171, 39681, 39682, 39683, 40705, 40706, 40707, 41729,
    41730, 41731, 42241, 42242, 42243, 43009, 43010, 43011, 43521, 43522, 43523,
    44289, 44290, 44291, 45057, 45058, 45059, 45313, 45314, 45315, 46849, 46850,
    46851, 47105, 47106, 47107, 47361, 47362, 47363, 47617, 47618, 47619, 47873,
    47874, 47875, 48129, 48130, 48131, 48641, 48642, 48643, 49409, 49410, 49411,
    49921, 49922, 49923,
];

const FUSED_CHUNK_FOLD_REMAINDER_KEYS: &[u32] = &[
    1, 2, 3, 1795, 3585, 3586, 3587, 5378, 7169, 7170, 7171, 8962, 10753, 10754,
    10755, 12547, 14337, 14338, 14339, 17921, 17922, 17923, 23299, 30466, 34051,
    37635, 41219, 44801, 44803, 48387, 55555, 62723, 66306, 77059, 84227, 87810,
    87811, 89344, 92928, 112899, 116482, 116483, 120066, 123651, 134403, 141571,
    145155, 148737, 148739, 150529, 150530, 150531, 152322, 154113, 154114, 154115,
    155907, 157697, 157698, 157699, 159489, 159490, 159491, 161024, 161281, 161282,
    161283, 163075, 164608, 164865, 164866, 164867, 166658, 166659, 167939, 168449,
    168450, 168451, 170240, 170242, 170243, 171523, 171776, 172033, 172034, 172035,
    173571, 173825, 173827, 175360, 175617, 175618, 175619, 177409, 177410, 177411,
    178691, 178944, 179201, 179202, 179203, 180993, 180994, 180995, 182275, 182528,
    182785, 182786, 182787, 184578, 184579, 185859, 186112, 186369, 186370, 186371,
    187907, 188161, 188162, 188163, 189441, 189443, 189696, 189953, 189954, 189955,
    191491, 191744, 191746, 191747, 193025, 193027, 193537, 193538, 193539, 195075,
    195328, 196609, 196610, 196611, 196864, 197121, 197122, 197123, 198657, 198658,
    198659, 198912, 200192, 200194, 200195, 200448, 202241, 202243, 202496, 203776,
    204032, 205826, 206080, 207361, 207362, 207363, 207616, 209410, 209411, 209664,
    210946, 210947, 213248, 214529, 214530, 214531, 214784, 216832, 218115, 218368,
    220419, 221699, 221952, 224001, 224002, 225283, 225536, 227586, 227587, 231170,
    231171, 232704, 234754, 236288, 238339, 239872, 241922, 245506, 245507, 247040,
    263426, 263427, 264960, 267011, 277762, 288515, 290048, 292098, 295683, 299267,
    304384, 306435, 313603, 320771, 322304, 324354, 324355, 329472, 331522, 345859,
    349442, 349443, 353027, 356611, 367362, 367363, 378115, 379648, 381699, 386816,
    392450, 396035, 399619,
];

fn fused_range_contains(ranges: &[(usize, usize, usize)], call_index: usize, bit: usize) -> bool {
    skip_structural_dead_fused_carries()
        && ranges
            .iter()
            .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn fused_key_contains(keys: &[u32], call_index: usize, bit: usize) -> bool {
    let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
    keys.binary_search(&key).is_ok()
}

const FUSED_BOUNDARY_ZERO_REMAINDER_KEYS: &[u32] = &[
    1282, 3842, 6402, 8962, 14082, 16642, 21762, 24322, 32001, 39682, 42242,
    47361, 47362, 49922, 51458, 52482, 55040, 55042, 57601, 62722, 70402,
    75521, 83202, 85761, 85762, 90882, 93442, 96002, 106241, 108802, 111362,
    113921, 113922, 116482, 118018, 119041, 119042, 121600, 121602, 124161,
    126721, 126722, 128256, 128258, 129281, 129282, 130818, 131840, 131841,
    131842, 132098, 133376, 133377, 133378, 134400, 134401, 134402, 135937,
    135938, 136960, 136962, 138497, 138498, 139520, 139521, 139522, 139777,
    141056, 141058, 141314, 142080, 142081, 142082, 142337, 142338, 143616,
    143617, 143618, 144640, 144641, 144642, 144897, 146177, 146178, 147200,
    147201, 147202, 147458, 148738, 149760, 149761, 149762, 150016, 150017,
    150018, 151297, 151298, 152320, 152321, 152322, 152578, 153857, 154880,
    154881, 154882, 156418, 158976, 158978, 160001, 160002, 160258, 162560,
    162561, 164096, 164098, 165122, 167682, 170241, 172801, 172802, 175362,
    177922, 180482, 183041, 183042, 185602, 188162, 190722, 193281, 193282,
    198402, 203521, 206082, 208641, 208642, 211201, 211202, 216322, 218882,
    226561, 229122, 231682, 234241, 236802, 244481, 249601, 249602, 252161,
    253698, 259842, 262402, 267521, 270082, 275202, 277762, 280322, 285442,
];

fn fused_boundary_zero_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    std::env::var_os("TLM_FUSED_SKIP_EXACT_BOUNDARY_ZERO").is_some()
        && fused_key_contains(FUSED_BOUNDARY_ZERO_REMAINDER_KEYS, call_index, bit)
}

fn fused_clean_fold_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    fused_range_contains(FUSED_CLEAN_FOLD_DEAD_RANGES, call_index, bit)
        || (skip_fused_clean_fold_top31() && bit == 31)
        || (skip_exact_fused_clean_fold()
            && fused_key_contains(FUSED_CLEAN_FOLD_REMAINDER_KEYS, call_index, bit))
}

fn fused_chunk_fold_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    fused_range_contains(FUSED_CHUNK_FOLD_DEAD_RANGES, call_index, bit)
        || (skip_exact_fused_chunk_fold()
            && fused_key_contains(FUSED_CHUNK_FOLD_REMAINDER_KEYS, call_index, bit))
}

fn fused_dirty_fold_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    skip_structural_dead_fused_dirty_fold()
        && FUSED_DIRTY_FOLD_DEAD_RANGES
            .iter()
            .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn fused_clean_window_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    skip_structural_dead_fused_clean_window()
        && FUSED_CLEAN_WINDOW_DEAD_RANGES
            .iter()
            .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn fold_call_reserve(index: usize, default: usize) -> usize {
    let base = env_index_value("TLM_TARGET_FOLD_CALL_RESERVES", index).unwrap_or(default);
    env_index_value("TLM_TARGET_FOLD_CALL_RESERVE_OVERRIDES", index).unwrap_or(base)
}

fn fold_ctl(p: usize) -> u8 {
    match p {
        0 | 4 | 6 | 32 => 1,
        1 | 5 | 33 => 2,
        7 => 3,
        8 | 9 => 4,
        10 => 5,
        11 => 6,
        _ => 0,
    }
}

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

fn toggle_dnot_e_from_intersection(
    circ: &mut B,
    d: &QubitId,
    cc: &QubitId,
    dne: &QubitId,
) {
    circ.cx(*d, *dne);
    circ.cx(*cc, *dne);
}

fn add_carry_into_tail_prefix(circ: &mut B, y: &[QubitId], c: &QubitId) {
    if std::env::var("TLM_FOLD_TAIL_CINC").ok().as_deref() == Some("1") {

        let yv: Vec<QubitId> = y.to_vec();
        super::mcx::cinc_khattar_gidney(circ, &yv, c);
        return;
    }

    let t = y.len();
    for k in (1..t).rev() {
        let mut ctrls: Vec<&QubitId> = Vec::with_capacity(k + 1);
        ctrls.push(c);
        ctrls.extend(y[..k].iter());
        super::mcx::mcx_clean_k(circ, &ctrls, &y[k]);
    }
    circ.cx(*c, y[0]);
}

fn add_mf_fold_clean(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId]) {
    add_mf_fold_clean_tail(circ, e, d, y, None);
}

fn add_mf_fold_clean_tail(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], tail_from: Option<usize>) {
    let l = y.len();
    assert!(l >= 2, "fold needs L >= 2");
    let loop_end = tail_from.unwrap_or(l - 1);
    const LAST_DERIVED: usize = 9;
    const LAST_AND: usize = 11;

    let mut cc = Some(circ.alloc_qubit());
    circ.ccx(*e, *d, *cc.as_ref().unwrap());
    let mut dne = Some(circ.alloc_qubit());
    toggle_dnot_e_from_intersection(
        circ,
        d,
        cc.as_ref().unwrap(),
        dne.as_ref().unwrap(),
    );
    let mut sxor = Some(circ.alloc_qubit());
    circ.cx(*e, *sxor.as_ref().unwrap());
    circ.cx(*d, *sxor.as_ref().unwrap());
    let mut sor = Some(circ.alloc_qubit());
    circ.cx(*sxor.as_ref().unwrap(), *sor.as_ref().unwrap());
    circ.cx(*cc.as_ref().unwrap(), *sor.as_ref().unwrap());

    fn fc<'a>(p: usize, e: &'a QubitId, d: &'a QubitId, cc: Option<&'a QubitId>, dne: Option<&'a QubitId>, sx: Option<&'a QubitId>, so: Option<&'a QubitId>) -> Option<&'a QubitId> {
        match fold_ctl(p) {
            1 => Some(e),
            2 => Some(d),
            3 => sx,
            4 => so,
            5 => dne,
            6 => cc,
            _ => None,
        }
    }

    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(l - 1);
    let c1 = circ.alloc_qubit();
    if let Some(a0) = fc(0, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
        if !skip_structural_dead_fused_carries() {
            let old_context = crate::point_add::set_op_trace_context(
                0x0d00_0000 | (((active_fold_call_index() as u32) & 0xffff) << 8),
            );
            circ.ccx(*a0, y[0], c1);
            crate::point_add::restore_op_trace_context(old_context);
        }
        circ.cx(*a0, y[0]);
    }
    cy.push(Some(c1));
    for i in 1..loop_end {
        let next = circ.alloc_qubit();
        {
            let ci = cy[i - 1].as_ref().unwrap();
            circ.cx(*ci, y[i]);
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, *ci);
            }
            if !fused_clean_fold_has_structurally_dead_carry(active_fold_call_index(), i) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0d00_0000 | (((active_fold_call_index() as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(y[i], *ci, next);
                crate::point_add::restore_op_trace_context(old_context);
            }
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, *ci);
            }
            circ.cx(*ci, next);
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, y[i]);
            }
        }
        cy.push(Some(next));
        if i == LAST_DERIVED {
            let so = sor.take().unwrap();
            circ.cx(*sxor.as_ref().unwrap(), so);
            circ.cx(*cc.as_ref().unwrap(), so);
            circ.zero_and_free(so);
            let sx = sxor.take().unwrap();
            circ.cx(*e, sx);
            circ.cx(*d, sx);
            circ.zero_and_free(sx);
        }
        if i == LAST_AND {
            let dn = dne.take().unwrap();
            toggle_dnot_e_from_intersection(circ, d, cc.as_ref().unwrap(), &dn);
            circ.zero_and_free(dn);
            let c = cc.take().unwrap();
            clear_and(circ, &c, e, d);
            circ.zero_and_free(c);
        }
    }
    match tail_from {
        None => {

            if let Some(at) = fc(l - 1, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*at, y[l - 1]);
            }
            circ.cx(*cy[l - 2].as_ref().unwrap(), y[l - 1]);
        }
        Some(nv) => {

            add_carry_into_tail_prefix(circ, &y[nv..], cy[nv - 1].as_ref().unwrap());
        }
    }

    for i in (1..loop_end).rev() {
        if i == LAST_AND {
            let c = circ.alloc_qubit();
            circ.ccx(*e, *d, c);
            cc = Some(c);
            let dn = circ.alloc_qubit();
            toggle_dnot_e_from_intersection(circ, d, cc.as_ref().unwrap(), &dn);
            dne = Some(dn);
        }
        if i == LAST_DERIVED {
            let sx = circ.alloc_qubit();
            circ.cx(*e, sx);
            circ.cx(*d, sx);
            let so = circ.alloc_qubit();
            circ.cx(sx, so);
            circ.cx(*cc.as_ref().unwrap(), so);
            sxor = Some(sx);
            sor = Some(so);
        }
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, y[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, ci);
        }

        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(y[i], ci, bit);
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, ci);
            circ.cx(*ai, y[i]);
        }
        cy[i - 1] = Some(ci);
    }

    let cy1 = cy[0].take().unwrap();
    if let Some(a0) = fc(0, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
        circ.cx(*a0, y[0]);
        let bit = circ.alloc_bit();
        circ.hmr(cy1, bit);
        circ.zero_and_free(cy1);
        circ.cz_if_bit(y[0], *a0, bit);
        circ.cx(*a0, y[0]);
    } else {
        circ.zero_and_free(cy1);
    }

    let sx = sxor.take().unwrap();
    let so = sor.take().unwrap();
    let cc = cc.take().unwrap();
    let dne = dne.take().unwrap();
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
    circ.zero_and_free(dne);
    circ.cx(sx, so);
    circ.cx(cc, so);
    circ.zero_and_free(so);
    circ.cx(*e, sx);
    circ.cx(*d, sx);
    circ.zero_and_free(sx);
    clear_and(circ, &cc, e, d);
    circ.zero_and_free(cc);
}

fn build_fold_controls(circ: &mut B, e: &QubitId, d: &QubitId) -> (QubitId, QubitId, QubitId, QubitId) {
    let cc = circ.alloc_qubit();
    circ.ccx(*e, *d, cc);
    let sxor = circ.alloc_qubit();
    circ.cx(*e, sxor);
    circ.cx(*d, sxor);
    let sor = circ.alloc_qubit();
    circ.cx(sxor, sor);
    circ.cx(cc, sor);
    let dne = circ.alloc_qubit();
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
    (cc, sxor, sor, dne)
}

fn uncompute_fold_controls(circ: &mut B, e: &QubitId, d: &QubitId, cc: QubitId, sxor: QubitId, sor: QubitId, dne: QubitId) {
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
    circ.zero_and_free(dne);
    circ.cx(sxor, sor);
    circ.cx(cc, sor);
    circ.zero_and_free(sor);
    circ.cx(*e, sxor);
    circ.cx(*d, sxor);
    circ.zero_and_free(sxor);
    clear_and(circ, &cc, e, d);
    circ.zero_and_free(cc);
}

fn fold_ctl_map(e: QubitId, d: QubitId, cc: QubitId, sxor: QubitId, sor: QubitId, dne: QubitId, l: usize) -> Vec<Option<QubitId>> {
    (0..l).map(|p| match fold_ctl(p) { 1 => Some(e), 2 => Some(d), 3 => Some(sxor), 4 => Some(sor), 5 => Some(dne), 6 => Some(cc), _ => None }).collect()
}

fn fold_chunk_clean(circ: &mut B, ctl: &[Option<QubitId>], y: &[QubitId], cin: Option<&QubitId>, cout: &QubitId) {
    let chunk_call_index = next_fold_chunk_call_index();
    let s = y.len();
    if s == 0 {
        if let Some(c) = cin { circ.cx(*c, *cout); }
        return;
    }
    let mut cy: Vec<Option<QubitId>> = (0..s - 1).map(|_| Some(circ.alloc_qubit())).collect();
    for i in 0..s {
        let on = ctl[i].as_ref();
        if i == 0 {
            let dst: QubitId = if s == 1 { *cout } else { *cy[0].as_ref().unwrap() };
            match cin {
                Some(c) => {
                    circ.cx(*c, y[0]);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    let old_context = crate::point_add::set_op_trace_context(
                        0x0e00_0000 | (((chunk_call_index as u32) & 0xffff) << 8),
                    );
                    circ.ccx(y[0], *c, dst);
                    crate::point_add::restore_op_trace_context(old_context);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    circ.cx(*c, dst);
                }
                None => {
                    if let Some(a) = on {
                        if !skip_structural_dead_fused_carries() {
                            let old_context = crate::point_add::set_op_trace_context(
                                0x0e00_0000 | (((chunk_call_index as u32) & 0xffff) << 8),
                            );
                            circ.ccx(*a, y[0], dst);
                            crate::point_add::restore_op_trace_context(old_context);
                        }
                    }
                }
            }
        } else {
            let ci: QubitId = *cy[i - 1].as_ref().unwrap();
            let dst: QubitId = if i == s - 1 { *cout } else { *cy[i].as_ref().unwrap() };
            circ.cx(ci, y[i]);
            if let Some(a) = on { circ.cx(*a, ci); }
            if !fused_chunk_fold_has_structurally_dead_carry(chunk_call_index, i) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0e00_0000 | (((chunk_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(y[i], ci, dst);
                crate::point_add::restore_op_trace_context(old_context);
            }
            if let Some(a) = on { circ.cx(*a, ci); }
            circ.cx(ci, dst);
        }
    }
    for i in 0..s {
        if let Some(a) = ctl[i].as_ref() { circ.cx(*a, y[i]); }
    }
    for i in (0..s - 1).rev() {
        let on = ctl[i].as_ref();
        if let Some(a) = on { circ.cx(*a, y[i]); }
        let next = cy[i].take().unwrap();
        if i == 0 {
            match cin {
                Some(c) => {
                    circ.cx(*c, next);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    let bit = circ.alloc_bit();
                    circ.hmr(next, bit); circ.zero_and_free(next);
                    circ.cz_if_bit(y[0], *c, bit);
                    if let Some(a) = on { circ.cx(*a, *c); circ.cx(*a, y[0]); }
                }
                None => {
                    let bit = circ.alloc_bit();
                    circ.hmr(next, bit); circ.zero_and_free(next);
                    if let Some(a) = on { circ.cz_if_bit(y[0], *a, bit); }
                    if let Some(a) = on { circ.cx(*a, y[0]); }
                }
            }
        } else {
            let ci: QubitId = *cy[i - 1].as_ref().unwrap();
            circ.cx(ci, next);
            if let Some(a) = on { circ.cx(*a, ci); }
            let bit = circ.alloc_bit();
            circ.hmr(next, bit); circ.zero_and_free(next);
            circ.cz_if_bit(y[i], ci, bit);
            if let Some(a) = on { circ.cx(*a, ci); circ.cx(*a, y[i]); }
        }
    }
}

fn fold_boundary_erase(circ: &mut B, ctl: &[Option<QubitId>], y: &[QubitId], cin: Option<&QubitId>, carry: QubitId) {
    if std::env::var("TLM_FOLD_BOUNDARY_ZERO_DIRECT")
        .ok()
        .as_deref()
        == Some("1")
        && cin.is_some()
        && ctl.iter().all(Option::is_none)
    {
        fold_boundary_erase_zero_direct(circ, y, cin.expect("cin checked"), carry);
        return;
    }
    let s = y.len();
    let temp: Vec<QubitId> = (0..s).map(|_| circ.alloc_qubit()).collect();
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    match cin {
        Some(cin) => super::arith::erase_carry_gated_opt(circ, None, y, &temp, cin, &carry, None),
        None => {
            super::arith::erase_carry_gated_zero_cin_opt(circ, None, y, &temp, &carry, None);
            circ.zero_and_free(carry);
        }
    }
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    for q in temp { circ.zero_and_free(q); }
}

fn fold_boundary_erase_zero_direct(circ: &mut B, y: &[QubitId], cin: &QubitId, carry: QubitId) {
    let boundary_call_index = next_fold_boundary_zero_call_index();
    let n = y.len();
    assert!(n >= 1, "zero boundary erase needs >= 1 bit");
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.zero_and_free(carry);
    circ.push_condition(bit);

    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(n);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0);
    cy.push(Some(c0));
    for i in 0..n - 1 {
        let next = circ.alloc_qubit();
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, y[i]);
        circ.x(*ci);
        let old_context = crate::point_add::set_op_trace_context(
            0x1c00_0000 | (((boundary_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        if !fused_boundary_zero_has_structurally_dead_carry(boundary_call_index, i) {
            circ.ccx(y[i], *ci, next);
        }
        crate::point_add::restore_op_trace_context(old_context);
        circ.x(*ci);
        circ.cx(*ci, next);
        cy.push(Some(next));
    }
    {
        let i = n - 1;
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, y[i]);
        circ.neg();
        circ.x(*ci);
        circ.cz(y[i], *ci);
        circ.x(*ci);
        circ.z(*ci);
        circ.cx(*ci, y[i]);
    }
    for i in (0..n - 1).rev() {
        let next = cy[i + 1].take().unwrap();
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, next);
        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.x(*ci);
        circ.cz_if_bit(y[i], *ci, mbit);
        circ.x(*ci);
        circ.cx(*ci, y[i]);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0);
    circ.x(c0);
    circ.zero_and_free(c0);
    circ.pop_condition();
}

fn add_mf_fold_chunked(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], s_chunk: usize) {
    let l = y.len();
    let release_controls = std::env::var("TLM_FOLD_RELEASE_CONTROLS")
        .ok()
        .as_deref()
        == Some("1");
    let zero_cin = std::env::var("TLM_FOLD_CHUNK_ZERO_CIN")
        .ok()
        .as_deref()
        == Some("1");
    let mut controls = Some(build_fold_controls(circ, e, d));
    let (cc, sxor, sor, dne) = controls.expect("fold controls present");
    let mut ctl = fold_ctl_map(*e, *d, cc, sxor, sor, dne, l);
    let cin0 = (!zero_cin).then(|| circ.alloc_qubit());
    let nch = l.div_ceil(s_chunk);
    let last_control_chunk = 11usize.min(l - 1) / s_chunk;
    let mut boundary: Vec<QubitId> = Vec::with_capacity(nch);
    for j in 0..nch {
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let cout = circ.alloc_qubit();
        let cin = if j == 0 {
            cin0.as_ref()
        } else {
            Some(&boundary[j - 1])
        };
        fold_chunk_clean(circ, &ctl[lo..hi], &y[lo..hi], cin, &cout);
        boundary.push(cout);
        if release_controls && j == last_control_chunk && j + 1 < nch {
            let (cc, sxor, sor, dne) = controls.take().expect("fold controls present");
            uncompute_fold_controls(circ, e, d, cc, sxor, sor, dne);
        }
    }
    for j in (0..nch).rev() {
        if release_controls && j == last_control_chunk && controls.is_none() {
            let rebuilt = build_fold_controls(circ, e, d);
            ctl = fold_ctl_map(*e, *d, rebuilt.0, rebuilt.1, rebuilt.2, rebuilt.3, l);
            controls = Some(rebuilt);
        }
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let bnd = boundary.pop().expect("boundary present");
        let cin = if j == 0 {
            cin0.as_ref()
        } else {
            Some(&boundary[j - 1])
        };
        fold_boundary_erase(circ, &ctl[lo..hi], &y[lo..hi], cin, bnd);
    }
    if let Some(cin0) = cin0 {
        circ.zero_and_free(cin0);
    }
    let (cc, sxor, sor, dne) = controls.take().expect("fold controls restored");
    uncompute_fold_controls(circ, e, d, cc, sxor, sor, dne);
}

enum OnCtl {
    None,
    E,
    D,
    Owned(QubitId),
}

fn on_ctl_apply(circ: &mut B, e: &QubitId, d: &QubitId, k: u8, q: &QubitId) {
    match k {
        3 => {
            circ.cx(*e, *q);
            circ.cx(*d, *q);
        }
        4 => {
            circ.x(*e);
            circ.x(*d);
            circ.ccx(*e, *d, *q);
            circ.x(*q);
            circ.x(*e);
            circ.x(*d);
        }
        5 => {
            circ.x(*e);
            circ.ccx(*e, *d, *q);
            circ.x(*e);
        }
        6 => circ.ccx(*e, *d, *q),
        _ => {}
    }
}

fn on_ctl(circ: &mut B, e: &QubitId, d: &QubitId, p: usize) -> OnCtl {
    match fold_ctl(p) {
        1 => OnCtl::E,
        2 => OnCtl::D,
        k @ (3 | 4 | 5 | 6) => {
            let q = circ.alloc_qubit();
            on_ctl_apply(circ, e, d, k, &q);
            OnCtl::Owned(q)
        }
        _ => OnCtl::None,
    }
}

fn on_ctl_ref(c: &OnCtl, e: &QubitId, d: &QubitId) -> Option<QubitId> {
    match c {
        OnCtl::None => None,
        OnCtl::E => Some(*e),
        OnCtl::D => Some(*d),
        OnCtl::Owned(q) => Some(*q),
    }
}

fn on_ctl_clear_nonlinear_hmr(
    circ: &mut B,
    e: &QubitId,
    d: &QubitId,
    k: u8,
    q: &QubitId,
) {
    let bit = circ.alloc_bit();
    circ.hmr(*q, bit);
    match k {
        4 => {
            circ.z_if_bit(*e, bit);
            circ.z_if_bit(*d, bit);
            circ.cz_if_bit(*e, *d, bit);
        }
        5 => {
            circ.z_if_bit(*d, bit);
            circ.cz_if_bit(*e, *d, bit);
        }
        6 => circ.cz_if_bit(*e, *d, bit),
        _ => unreachable!("HMR clear requires a nonlinear fold control"),
    }
}

fn on_ctl_free(circ: &mut B, e: &QubitId, d: &QubitId, p: usize, c: OnCtl) {
    if let OnCtl::Owned(q) = c {
        let k = fold_ctl(p);
        let hmr_disabled = std::env::var("TLM_FOLD_HMR_CONTROL_CLEANUP_DISABLE")
            .ok()
            .as_deref()
            == Some("1");
        if k == 3 || hmr_disabled {

            on_ctl_apply(circ, e, d, k, &q);
        } else {
            on_ctl_clear_nonlinear_hmr(circ, e, d, k, &q);
        }
        circ.zero_and_free(q);
    }
}

fn xor_carries_perpos(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], out: &[QubitId], carry_in: Option<&QubitId>) {
    let n = y.len();
    fn ccx_cond(circ: &mut B, aq: Option<&QubitId>, c1: &QubitId, c2: &QubitId, t: &QubitId, g0: bool, g1: bool) {
        if let Some(a) = aq {
            if g0 {
                circ.cx(*a, *c1);
            }
            if g1 {
                circ.cx(*a, *c2);
            }
        }
        circ.ccx(*c1, *c2, *t);
        if let Some(a) = aq {
            if g0 {
                circ.cx(*a, *c1);
            }
            if g1 {
                circ.cx(*a, *c2);
            }
        }
    }
    for i in (1..n - 1).rev() {
        let c = on_ctl(circ, e, d, base + i);
        let aq = on_ctl_ref(&c, e, d);
        let g0 = aq.is_some();
        ccx_cond(circ, aq.as_ref(), &y[i], &out[i - 1], &out[i], g0, false);
        on_ctl_free(circ, e, d, base + i, c);
    }
    for i in 0..n - 1 {
        let c = on_ctl(circ, e, d, base + i);
        if let Some(a) = on_ctl_ref(&c, e, d) {
            circ.cx(a, out[i]);
        }
        on_ctl_free(circ, e, d, base + i, c);
    }
    {
        let c = on_ctl(circ, e, d, base);
        let aq = on_ctl_ref(&c, e, d);
        let g = aq.is_some();
        match carry_in {
            Some(cy) => ccx_cond(circ, aq.as_ref(), cy, &y[0], &out[0], g, g),
            None => {
                let cin = circ.alloc_qubit();
                ccx_cond(circ, aq.as_ref(), &cin, &y[0], &out[0], g, g);
                circ.zero_and_free(cin);
            }
        }
        on_ctl_free(circ, e, d, base, c);
    }
    for i in 1..n - 1 {
        let c = on_ctl(circ, e, d, base + i);
        let aq = on_ctl_ref(&c, e, d);
        let gi = aq.is_some();
        ccx_cond(circ, aq.as_ref(), &y[i], &out[i - 1], &out[i], gi, gi);
        on_ctl_free(circ, e, d, base + i, c);
    }
}

fn dirty_body(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], dirty: &[QubitId], carry_in: Option<&QubitId>) {
    let dirty_call_index = next_fold_dirty_call_index();
    let l = y.len();
    assert!(l >= 2);
    assert!(dirty.len() >= l - 1, "need L-1 borrowed dirty bits");
    let mut cin_owned = if carry_in.is_none() { Some(circ.alloc_qubit()) } else { None };
    let mut bits: Vec<BitId> = Vec::with_capacity(l - 1);
    let mut prev_new: Option<QubitId> = None;
    for i in 0..l - 1 {
        let new = circ.alloc_qubit();
        let anc = circ.alloc_qubit();
        let ctlh = on_ctl(circ, e, d, base + i);
        {
            let cyi: QubitId = if i == 0 {
                carry_in.copied().unwrap_or_else(|| *cin_owned.as_ref().unwrap())
            } else {
                *prev_new.as_ref().unwrap()
            };
            if let Some(ai) = on_ctl_ref(&ctlh, e, d) {
                circ.cx(ai, anc);
            }
            circ.cx(cyi, anc);
            circ.cx(cyi, y[i]);
            let old_context = crate::point_add::set_op_trace_context(
                0x0f00_0000 | (((dirty_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            if !fused_dirty_fold_has_structurally_dead_carry(dirty_call_index, i) {
                circ.ccx(y[i], anc, new);
            }
            crate::point_add::restore_op_trace_context(old_context);
            circ.cx(cyi, new);
            circ.cx(new, dirty[i]);
            circ.cx(cyi, anc);
            if let Some(ai) = on_ctl_ref(&ctlh, e, d) {
                circ.cx(ai, anc);
                circ.cx(ai, y[i]);
            }
        }
        on_ctl_free(circ, e, d, base + i, ctlh);
        circ.zero_and_free(anc);
        if i == 0 {
            if let Some(c) = cin_owned.take() {
                circ.zero_and_free(c);
            }
        } else {
            let b = circ.alloc_bit();
            circ.hmr(*prev_new.as_ref().unwrap(), b);
            circ.zero_and_free(prev_new.take().unwrap());
            bits.push(b);
        }
        prev_new = Some(new);
    }
    let cy_top = prev_new.take().unwrap();
    {
        let topc = on_ctl(circ, e, d, base + l - 1);
        if let Some(at) = on_ctl_ref(&topc, e, d) {
            circ.cx(at, y[l - 1]);
        }
        on_ctl_free(circ, e, d, base + l - 1, topc);
    }
    circ.cx(cy_top, y[l - 1]);
    let b = circ.alloc_bit();
    circ.hmr(cy_top, b);
    circ.zero_and_free(cy_top);
    bits.push(b);

    for i in 0..l - 1 {
        circ.z_if_bit(dirty[i], bits[i]);
    }
    for q in y {
        circ.x(*q);
    }
    xor_carries_perpos(circ, e, d, base, y, dirty, carry_in);
    for q in y {
        circ.x(*q);
    }
    for i in 0..l - 1 {
        circ.z_if_bit(dirty[i], bits[i]);
    }
}

fn clean_window_fwd(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], carries: &[QubitId]) {
    let clean_window_call_index = next_fold_clean_window_call_index();
    let b = y.len();
    assert_eq!(carries.len(), b);
    {
        let c0 = on_ctl(circ, e, d, base);
        if let Some(a0) = on_ctl_ref(&c0, e, d) {
            let old_context = crate::point_add::set_op_trace_context(
                0x1000_0000 | (((clean_window_call_index as u32) & 0xffff) << 8),
            );
            if !fused_clean_window_has_structurally_dead_carry(clean_window_call_index, 0) {
                circ.ccx(a0, y[0], carries[0]);
            }
            crate::point_add::restore_op_trace_context(old_context);
        }
        on_ctl_free(circ, e, d, base, c0);
    }
    for i in 1..b {
        let ci = on_ctl(circ, e, d, base + i);
        let ai = on_ctl_ref(&ci, e, d);
        circ.cx(carries[i - 1], y[i]);
        if let Some(a) = &ai {
            circ.cx(*a, carries[i - 1]);
        }
        let old_context = crate::point_add::set_op_trace_context(
            0x1000_0000 | (((clean_window_call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        if !fused_clean_window_has_structurally_dead_carry(clean_window_call_index, i) {
            circ.ccx(y[i], carries[i - 1], carries[i]);
        }
        crate::point_add::restore_op_trace_context(old_context);
        if let Some(a) = &ai {
            circ.cx(*a, carries[i - 1]);
        }
        circ.cx(carries[i - 1], carries[i]);
        on_ctl_free(circ, e, d, base + i, ci);
    }
    for i in 0..b {
        let ci = on_ctl(circ, e, d, base + i);
        if let Some(a) = on_ctl_ref(&ci, e, d) {
            circ.cx(a, y[i]);
        }
        on_ctl_free(circ, e, d, base + i, ci);
    }
}

fn clean_window_rev(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], carries: Vec<QubitId>) {
    let b = y.len();
    let mut cy: Vec<Option<QubitId>> = carries.into_iter().map(Some).collect();
    for i in (1..b).rev() {
        let ci_ctl = on_ctl(circ, e, d, base + i);
        let actl = on_ctl_ref(&ci_ctl, e, d);
        if let Some(ai) = &actl {
            circ.cx(*ai, y[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if let Some(ai) = &actl {
            circ.cx(*ai, ci);
        }
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(y[i], ci, bit);
        if let Some(ai) = &actl {
            circ.cx(*ai, ci);
            circ.cx(*ai, y[i]);
        }
        on_ctl_free(circ, e, d, base + i, ci_ctl);
        cy[i - 1] = Some(ci);
    }
    let cy0 = cy[0].take().unwrap();
    let c0 = on_ctl(circ, e, d, base);
    if let Some(a0) = on_ctl_ref(&c0, e, d) {
        circ.cx(a0, y[0]);
        let bit = circ.alloc_bit();
        circ.hmr(cy0, bit);
        circ.zero_and_free(cy0);
        circ.cz_if_bit(y[0], a0, bit);
        circ.cx(a0, y[0]);
    } else {
        circ.zero_and_free(cy0);
    }
    on_ctl_free(circ, e, d, base, c0);
}

fn build_fold_at(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], dirty: &[QubitId], nv: usize) {
    let l = y.len();
    if nv >= l - 1 {

        add_mf_fold_clean(circ, e, d, y);
        return;
    }

    const PROP_FROM: usize = 34;
    if nv >= 1 && nv >= PROP_FROM {
        add_mf_fold_clean_tail(circ, e, d, y, Some(nv));
        return;
    }
    if nv == 0 {
        dirty_body(circ, e, d, 0, y, dirty, None);
    } else {
        let carries: Vec<QubitId> = (0..nv).map(|_| circ.alloc_qubit()).collect();
        clean_window_fwd(circ, e, d, 0, &y[..nv], &carries);
        let cin = carries[nv - 1];
        dirty_body(circ, e, d, nv, &y[nv..], &dirty[nv..], Some(&cin));
        clean_window_rev(circ, e, d, 0, &y[..nv], carries);
    }
}

fn fused_fold(circ: &mut B, e: &QubitId, d: &QubitId, ylow: &[QubitId], dirty: &[QubitId]) {
    let call_index = next_fold_call_index();
    let prior_fold_call_index = enter_fold_call_index(call_index);
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let code = super::next_fold();
    let mut selected_nv = None;
    if code < 0 {
        let chunk = std::env::var("TLM_FOLD_CHUNK_FORCE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .filter(|&value| value > 0)
            .unwrap_or((-code) as usize);
        add_mf_fold_chunked(circ, e, d, ylow, chunk);
    } else {
        let default_reserve = std::env::var("TLM_TARGET_FOLD_RESERVE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(4);
        let reserve = fold_call_reserve(call_index, default_reserve);
        let nv = super::target_qubit_headroom(circ)
            .map_or(code as usize, |headroom| {
                (code as usize).min(headroom.saturating_sub(reserve))
            });
        selected_nv = Some(nv);
        build_fold_at(circ, e, d, ylow, dirty, nv);
    }
    restore_fold_call_index(prior_fold_call_index);
    if std::env::var_os("TRACE_TLM_FOLD").is_some() {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(circ.active_qubits);
        eprintln!(
            "TLM_FOLD call={} phase={} code={} nv={} entry_active={} local_peak={} ops={}",
            call_index,
            circ.phase,
            code,
            selected_nv.map_or(-1, |value| value as i32),
            entry_active,
            local_peak,
            circ.current_ops_len(),
        );
    }
}

fn fused_fold_e_only(circ: &mut B, e: &QubitId, y: &[QubitId]) {
    let call_index = next_fold_call_index();
    let prior_fold_call_index = enter_fold_call_index(call_index);
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let code = super::next_fold();
    let default_reserve = std::env::var("TLM_TARGET_FOLD_RESERVE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(4);
    let reserve = fold_call_reserve(call_index, default_reserve);
    let g = if code < 0 {
        0
    } else {
        super::target_qubit_headroom(circ)
            .map_or(code as usize, |headroom| {
                (code as usize).min(headroom.saturating_sub(reserve))
            })
            .min(LSBS - 1)
    };
    let f_bytes = F_SECP256K1.to_le_bytes();
    super::arith::add_f_window_pub(circ, e, y, LSBS, &f_bytes, Some(g));
    restore_fold_call_index(prior_fold_call_index);
    if std::env::var_os("TRACE_TLM_FOLD").is_some() {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(circ.active_qubits);
        eprintln!(
            "TLM_FOLD call={} phase={} code={} nv={} entry_active={} local_peak={} ops={}",
            call_index,
            circ.phase,
            code,
            g as i32,
            entry_active,
            local_peak,
            circ.current_ops_len(),
        );
    }
}

fn trace_fold_alloc(circ: &B, name: &str, stage: &str, i: usize) {
    if std::env::var_os("TRACE_TLM_FOLD_ALLOC").is_some() {
        eprintln!(
            "TLM_FOLD_ALLOC name={name} stage={stage} i={i} active={} ops={}",
            circ.active_qubits,
            circ.current_ops_len(),
        );
    }
}

pub fn fused_double_cdouble(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    let shift_call_index = next_fused_cdouble_fwd_shift_call_index();
    maybe_run_gradual_fold_nonlinear_control_hmr_selftest();
    let n = 256usize;
    assert_eq!(y.len(), n, "fused double expects 256-bit y (transient overflow)");
    let _ = F_SECP256K1;
    trace_fold_alloc(circ, "fwd_cdouble", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_cdouble", "after_hi", usize::MAX);
    let hi2 = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_cdouble", "after_hi2", usize::MAX);

    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    w.push(hi2);

    for i in (1..w.len()).rev() {
        circ.swap(w[i], w[i - 1]);
    }

    for i in (1..w.len()).rev() {
        let bit = i - 1;
        let old_context = crate::point_add::set_op_trace_context(
            0x1400_0000 | (((shift_call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff),
        );
        if !(bit == 0 && skip_structural_dead_fused_cdouble_shift0()) {
            circ.cswap(*s2, w[i], w[i - 1]);
        }
        crate::point_add::restore_op_trace_context(old_context);
    }

    let borrow: Vec<QubitId> = y[LSBS..2 * LSBS - 1].to_vec();
    fused_fold(circ, &w[n], &w[n + 1], &y[..LSBS], &borrow);

    circ.cx(y[0], w[n]);
    clear_and(circ, &w[n + 1], s2, &y[1]);
    circ.zero_and_free(hi);
    circ.zero_and_free(hi2);
}

pub fn fused_double_only(circ: &mut B, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused double expects 256-bit y");
    trace_fold_alloc(circ, "fwd_only", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_only", "after_hi", usize::MAX);
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    for i in (1..w.len()).rev() {
        circ.swap(w[i], w[i - 1]);
    }
    fused_fold_e_only(circ, &w[n], y);
    circ.cx(y[0], w[n]);
    circ.zero_and_free(hi);
}

pub fn fused_double_cdouble_reverse(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    let shift_call_index = next_fused_cdouble_rev_shift_call_index();
    maybe_run_gradual_fold_nonlinear_control_hmr_selftest();
    let n = 256usize;
    assert_eq!(y.len(), n, "fused halve expects 256-bit y (transient overflow)");
    trace_fold_alloc(circ, "rev_cdouble", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_cdouble", "after_hi", usize::MAX);
    let hi2 = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_cdouble", "after_hi2", usize::MAX);
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    w.push(hi2);

    circ.ccx(*s2, y[1], w[n + 1]);
    circ.cx(y[0], w[n]);

    let borrow: Vec<QubitId> = y[LSBS..2 * LSBS - 1].to_vec();
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    fused_fold(circ, &w[n], &w[n + 1], &y[..LSBS], &borrow);
    for q in &y[..LSBS] {
        circ.x(*q);
    }

    for i in 1..w.len() {
        let bit = i - 1;
        let old_context = crate::point_add::set_op_trace_context(
            0x1500_0000 | (((shift_call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff),
        );
        if !(bit == 0 && skip_structural_dead_fused_cdouble_shift0()) {
            circ.cswap(*s2, w[i], w[i - 1]);
        }
        crate::point_add::restore_op_trace_context(old_context);
    }
    for i in 1..w.len() {
        circ.swap(w[i], w[i - 1]);
    }
    circ.zero_and_free(hi);
    circ.zero_and_free(hi2);
}

pub fn fused_double_only_reverse(circ: &mut B, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused halve expects 256-bit y");
    trace_fold_alloc(circ, "rev_only", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_only", "after_hi", usize::MAX);
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    circ.cx(y[0], w[n]);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    fused_fold_e_only(circ, &w[n], y);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    for i in 1..w.len() {
        circ.swap(w[i], w[i - 1]);
    }
    circ.zero_and_free(hi);
}

fn gradual_fold_nonlinear_control_hmr_selftest() {
    use crate::circuit::OperationType;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    for &(position, kind) in &[(8usize, 4u8), (10, 5), (11, 6)] {
        assert_eq!(fold_ctl(position), kind);

        let mut circ = B::new();
        let e = circ.alloc_qubit();
        let d = circ.alloc_qubit();
        let q = circ.alloc_qubit();
        on_ctl_apply(&mut circ, &e, &d, kind, &q);
        on_ctl_free(&mut circ, &e, &d, position, OnCtl::Owned(q));

        assert_eq!(circ.active_qubits, 2, "owned control was not released");
        assert_eq!(circ.peak_qubits, 3, "cleanup increased peak width");
        assert_eq!(circ.next_bit, 1, "expected one HMR result bit");
        assert_eq!(
            circ.ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count(),
            1,
            "cleanup must add no Toffoli-class gate",
        );
        assert_eq!(
            circ.ops
                .iter()
                .filter(|op| op.kind == OperationType::Hmr)
                .count(),
            1,
        );

        let mut e_mask = 0u64;
        let mut d_mask = 0u64;
        for shot in 0..64usize {
            let state = shot & 3;
            e_mask |= ((state & 1) as u64) << shot;
            d_mask |= (((state >> 1) & 1) as u64) << shot;
        }

        let mut seed = Shake128::default();
        seed.update(b"gradual-fold-derived-control-hmr");
        seed.update(&[kind]);
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(
            circ.next_qubit as usize,
            circ.next_bit as usize,
            &mut xof,
        );
        *sim.qubit_mut(e) = e_mask;
        *sim.qubit_mut(d) = d_mask;
        sim.apply_iter(circ.ops.iter());

        assert_eq!(sim.qubit(e), e_mask, "e changed for control kind {kind}");
        assert_eq!(sim.qubit(d), d_mask, "d changed for control kind {kind}");
        assert_eq!(sim.qubit(q), 0, "owned control remained dirty for kind {kind}");
        assert_eq!(sim.phase, 0, "phase feedback failed for control kind {kind}");

        let measured = sim.bits[0];
        for state in 0..4usize {
            let mut outcomes = 0u8;
            for shot in (state..64usize).step_by(4) {
                outcomes |= 1 << ((measured >> shot) & 1);
            }
            assert_eq!(
                outcomes, 0b11,
                "HMR outcomes not exhaustive for kind {kind}, state {state}",
            );
        }
    }
}

fn maybe_run_gradual_fold_nonlinear_control_hmr_selftest() {
    if std::env::var_os("TLM_FOLD_HMR_CONTROL_SELFTEST").is_none() {
        return;
    }
    static SELFTEST: std::sync::Once = std::sync::Once::new();
    SELFTEST.call_once(|| {
        gradual_fold_nonlinear_control_hmr_selftest();
        eprintln!("TLM_FOLD_HMR_CONTROL_SELFTEST_OK");
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradual_fold_nonlinear_control_hmr_cleanup_is_exact() {
        gradual_fold_nonlinear_control_hmr_selftest();
    }
}
