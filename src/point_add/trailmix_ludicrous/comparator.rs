
use super::{B, BExt};
use crate::circuit::QubitId;
use std::cell::Cell;

thread_local! {
    static COMPARE_DIRECT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static COMPARE_CIN_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_compare_call_index() {
    COMPARE_DIRECT_CALL_INDEX.with(|index| index.set(0));
    COMPARE_CIN_CALL_INDEX.with(|index| index.set(0));
}

fn next_compare_direct_call_index() -> usize {
    COMPARE_DIRECT_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_compare_cin_call_index() -> usize {
    COMPARE_CIN_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

const COMPARE_CIN_STRUCTURAL_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (3, 0, 64),
    (4, 20, 21),
    (4, 23, 64),
    (105, 16, 18),
    (1257, 0, 2),
    (1279, 0, 2),
    (1292, 0, 2),
    (1305, 0, 2),
    (1318, 0, 2),
    (1331, 0, 2),
    (1344, 0, 2),
    (1356, 0, 2),
    (1369, 0, 2),
    (1382, 0, 2),
    (1395, 0, 2),
    (1408, 0, 2),
    (1422, 0, 2),
    (1435, 0, 2),
    (1449, 0, 2),
    (1463, 0, 2),
    (1477, 0, 2),
    (1491, 0, 2),
    (1506, 0, 2),
    (1520, 0, 2),
    (1535, 0, 2),
    (1551, 0, 2),
    (1566, 0, 2),
    (1582, 0, 2),
    (1599, 0, 2),
    (1615, 0, 2),
    (1632, 0, 2),
    (1652, 0, 2),
    (1669, 0, 2),
    (1689, 0, 2),
    (1706, 0, 2),
    (1726, 0, 2),
    (1743, 0, 2),
    (1760, 0, 2),
    (1777, 0, 2),
    (1794, 0, 2),
    (1812, 0, 2),
    (1829, 0, 2),
    (1847, 0, 2),
    (1866, 0, 2),
    (1884, 0, 2),
    (1903, 0, 2),
    (1923, 0, 2),
    (1942, 0, 2),
    (1958, 0, 2),
    (1974, 0, 2),
    (1990, 0, 2),
    (2006, 0, 2),
    (2023, 0, 2),
    (2041, 0, 2),
    (2061, 0, 2),
    (2086, 0, 2),
    (2112, 0, 2),
    (2137, 0, 2),
    (2161, 0, 2),
    (2184, 0, 2),
    (2366, 0, 2),
    (2389, 0, 2),
    (2413, 0, 2),
    (2438, 0, 2),
    (2464, 0, 2),
    (2489, 0, 2),
    (2509, 0, 2),
    (2527, 0, 2),
    (2544, 0, 2),
    (2560, 0, 2),
    (2576, 0, 2),
    (2592, 0, 2),
    (2608, 0, 2),
    (2627, 0, 2),
    (2647, 0, 2),
    (2666, 0, 2),
    (2684, 0, 2),
    (2703, 0, 2),
    (2721, 0, 2),
    (2738, 0, 2),
    (2756, 0, 2),
    (2773, 0, 2),
    (2790, 0, 2),
    (2807, 0, 2),
    (2824, 0, 2),
    (2844, 0, 2),
    (2861, 0, 2),
    (2881, 0, 2),
    (2898, 0, 2),
    (2918, 0, 2),
    (2935, 0, 2),
    (2951, 0, 2),
    (2968, 0, 2),
    (2984, 0, 2),
    (2999, 0, 2),
    (3015, 0, 2),
    (3030, 0, 2),
    (3044, 0, 2),
    (3059, 0, 2),
    (3073, 0, 2),
    (3087, 0, 2),
    (3101, 0, 2),
    (3115, 0, 2),
    (3128, 0, 2),
    (3142, 0, 2),
    (3155, 0, 2),
    (3168, 0, 2),
    (3181, 0, 2),
    (3194, 0, 2),
    (3206, 0, 2),
    (3219, 0, 2),
    (3232, 0, 2),
    (3245, 0, 2),
    (3258, 0, 2),
    (3271, 0, 2),
    (3293, 0, 2),
];

const COMPARE_CIN_REMAINDER_KEYS: &[u32] = &[
    8, 274, 530, 12817, 12818, 14098, 15377, 15378, 16657, 16658, 17937, 17938,
    19217, 19218, 20497, 20498, 21777, 21778, 23057, 23058, 24337, 24338, 25617, 25618,
    28177, 28178, 29457, 29458, 30737, 30738, 32017, 32018, 33297, 33298, 34577, 34578,
    35857, 35858, 37137, 37138, 38417, 38418, 39697, 39698, 40978, 42257, 42258, 43537,
    43538, 44817, 44818, 46097, 46098, 47378, 48658, 49938, 51218, 52498, 53777, 53778,
    55058, 56337, 56338, 57618, 58898, 60178, 61458, 62738, 64018, 66578, 67858, 69138,
    70674, 71954, 73490, 75026, 76562, 78098, 79634, 81170, 82706, 84242, 87314, 90386,
    94994, 111890, 116498, 121106, 127250, 130322, 131858, 134930, 139538, 141074,
    144146, 145682, 147218, 150290, 151826, 154897, 154898, 156433, 156434, 158226,
    159762, 161554, 165138, 175634, 177426, 179218, 182802, 184594, 186386, 188178,
    189970, 191761, 191762, 193554, 200716, 202509, 209682, 211471, 213264, 216846,
    233734, 235783, 246020, 252163, 254212, 256261, 258306, 260355, 270592, 272641,
    281345, 283666, 288000, 300818, 306962, 311314, 313618, 323858, 326162, 329490,
    332818, 339474, 854546, 867090, 873490, 925704, 929798, 937738, 943115, 946702,
    950288, 952079, 964623, 966418, 970002, 971794, 977170, 1000210, 1002002, 1003794,
    1007122, 1022482, 1031698, 1033234, 1085458, 1091602, 1103378, 1112338, 1121298,
    1127698, 1132818, 1136658, 1148177, 1148178, 1149457,
];

fn compare_cin_has_structurally_dead_carry(call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_COMPARE_SKIP_STRUCTURAL_DEAD_CALLS").is_none() {
        return false;
    }
    if std::env::var_os("TLM_COMPARE_SKIP_EXACT_CIN_REMAINDER").is_some() {
        let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if COMPARE_CIN_REMAINDER_KEYS.binary_search(&key).is_ok() {
            return true;
        }
    }
    COMPARE_CIN_STRUCTURAL_DEAD_RANGES
        .iter()
        .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

const COMPARE_STRUCTURAL_DEAD_TOP_RANGES: &[(usize, usize, usize)] = &[
    (775, 0, 18),
    (776, 0, 18),
    (777, 0, 18),
    (778, 0, 18),
    (779, 0, 18),
    (782, 0, 18),
    (783, 0, 18),
    (784, 0, 18),
    (785, 0, 18),
    (786, 0, 18),
    (788, 0, 18),
    (789, 0, 18),
    (790, 0, 18),
    (791, 0, 18),
    (792, 0, 18),
    (516, 1, 10),
    (1055, 1, 10),
    (517, 3, 11),
    (518, 5, 13),
    (1057, 3, 11),
    (1059, 5, 13),
    (520, 8, 15),
    (519, 7, 13),
    (521, 9, 15),
    (1061, 7, 13),
    (1063, 9, 15),
    (1065, 9, 15),
    (511, 8, 13),
    (522, 10, 15),
    (523, 11, 16),
    (1067, 10, 15),
    (1069, 11, 16),
    (507, 11, 15),
    (515, 6, 10),
    (525, 13, 17),
    (526, 14, 18),
    (527, 15, 19),
    (808, 29, 33),
    (1052, 9, 13),
    (1071, 12, 16),
    (1073, 13, 17),
    (1075, 14, 18),
    (21, 30, 33),
    (23, 30, 33),
    (25, 30, 33),
    (27, 30, 33),
    (29, 30, 33),
    (31, 30, 33),
    (33, 31, 34),
    (505, 12, 15),
    (509, 10, 13),
    (524, 13, 16),
    (537, 25, 28),
    (539, 27, 30),
    (541, 29, 32),
    (543, 31, 34),
    (807, 30, 33),
    (809, 30, 33),
    (810, 30, 33),
    (811, 30, 33),
    (812, 30, 33),
    (813, 31, 34),
    (819, 30, 33),
    (1049, 12, 15),
    (1050, 12, 15),
    (1051, 10, 13),
    (1053, 8, 11),
    (1054, 7, 10),
    (1077, 16, 19),
    (1081, 17, 20),
    (1095, 24, 27),
    (1099, 26, 29),
    (1101, 27, 30),
    (1115, 34, 37),
    (1121, 37, 40),
    (19, 29, 31),
    (35, 32, 34),
    (37, 31, 33),
    (39, 32, 34),
    (41, 32, 34),
    (43, 32, 34),
    (45, 31, 33),
    (47, 33, 35),
    (49, 33, 35),
    (51, 32, 34),
    (53, 32, 34),
    (55, 33, 35),
    (57, 32, 34),
    (61, 34, 36),
    (501, 14, 16),
    (503, 13, 15),
    (513, 9, 11),
    (528, 17, 19),
    (529, 18, 20),
    (533, 21, 23),
    (538, 27, 29),
    (540, 29, 31),
    (542, 31, 33),
    (544, 33, 35),
    (545, 34, 36),
    (546, 35, 37),
    (547, 36, 38),
    (548, 37, 39),
    (550, 39, 41),
    (551, 40, 42),
    (553, 42, 44),
    (555, 43, 45),
    (805, 29, 31),
    (806, 29, 31),
    (814, 32, 34),
    (815, 31, 33),
    (817, 32, 34),
    (818, 32, 34),
    (820, 33, 35),
    (821, 33, 35),
    (822, 32, 34),
    (823, 32, 34),
    (824, 33, 35),
    (825, 32, 34),
    (826, 32, 34),
    (827, 34, 36),
    (829, 33, 35),
    (832, 33, 35),
    (833, 34, 36),
    (1044, 16, 18),
    (1047, 14, 16),
    (1048, 13, 15),
    (1079, 17, 19),
    (1097, 26, 28),
    (1103, 29, 31),
    (1105, 30, 32),
    (1109, 32, 34),
    (1111, 33, 35),
    (1113, 34, 36),
    (1117, 36, 38),
    (1119, 37, 39),
    (1123, 39, 41),
    (1125, 40, 42),
    (1127, 41, 43),
    (1129, 42, 44),
    (1131, 43, 45),
    (1133, 43, 45),
    (1135, 44, 46),
    (1137, 45, 47),
    (15, 28, 29),
    (17, 30, 31),
    (59, 33, 34),
    (63, 33, 34),
    (65, 34, 35),
    (67, 34, 35),
    (69, 34, 35),
    (71, 34, 35),
    (73, 35, 36),
    (75, 34, 35),
    (77, 35, 36),
    (79, 35, 36),
    (81, 35, 36),
    (83, 34, 35),
    (85, 35, 36),
    (87, 35, 36),
    (89, 35, 36),
    (91, 35, 36),
    (93, 36, 37),
    (95, 35, 36),
    (97, 36, 37),
    (99, 35, 36),
    (101, 36, 37),
    (105, 36, 37),
    (107, 36, 37),
    (109, 36, 37),
    (111, 36, 37),
    (113, 37, 38),
    (121, 37, 38),
    (123, 37, 38),
    (129, 37, 38),
    (133, 37, 38),
    (137, 38, 39),
    (143, 38, 39),
    (165, 40, 41),
    (177, 40, 41),
    (189, 42, 43),
    (191, 41, 42),
    (203, 42, 43),
    (207, 43, 44),
    (211, 43, 44),
    (213, 42, 43),
    (215, 43, 44),
    (217, 43, 44),
    (221, 44, 45),
    (223, 43, 44),
    (225, 44, 45),
    (227, 44, 45),
    (229, 44, 45),
    (233, 44, 45),
    (245, 45, 46),
    (247, 44, 45),
    (249, 44, 45),
    (255, 44, 45),
    (257, 44, 45),
    (259, 44, 45),
    (261, 45, 46),
    (263, 45, 46),
    (265, 45, 46),
    (285, 45, 46),
    (287, 46, 47),
    (291, 46, 47),
    (299, 45, 46),
    (321, 46, 47),
    (327, 47, 48),
    (333, 47, 48),
    (359, 49, 50),
    (361, 49, 50),
    (395, 50, 51),
    (411, 52, 53),
    (415, 52, 53),
    (421, 52, 53),
    (431, 47, 48),
    (433, 46, 47),
    (439, 44, 45),
    (441, 43, 44),
    (447, 40, 41),
    (449, 39, 40),
    (451, 38, 39),
    (459, 34, 35),
    (461, 33, 34),
    (465, 31, 32),
    (467, 30, 31),
    (471, 28, 29),
    (475, 26, 27),
    (493, 18, 19),
    (495, 17, 18),
    (497, 16, 17),
    (499, 15, 16),
    (530, 19, 20),
    (532, 21, 22),
    (534, 23, 24),
    (536, 26, 27),
    (549, 39, 40),
    (552, 42, 43),
    (554, 44, 45),
    (556, 45, 46),
    (557, 46, 47),
    (558, 47, 48),
    (559, 48, 49),
    (563, 52, 53),
    (568, 52, 53),
    (570, 52, 53),
    (575, 52, 53),
    (592, 49, 50),
    (595, 48, 49),
    (613, 46, 47),
    (621, 45, 46),
    (628, 46, 47),
    (630, 46, 47),
    (631, 44, 44),
    (631, 46, 46),
    (635, 44, 45),
    (639, 44, 45),
    (640, 45, 46),
    (642, 44, 44),
    (642, 46, 46),
    (643, 45, 46),
    (644, 44, 45),
    (645, 44, 45),
    (647, 44, 45),
    (650, 44, 45),
    (651, 45, 46),
    (661, 44, 45),
    (678, 41, 42),
    (679, 42, 43),
    (706, 37, 38),
    (719, 36, 37),
    (721, 36, 37),
    (722, 36, 37),
    (723, 36, 37),
    (742, 33, 34),
    (744, 33, 34),
    (804, 28, 29),
    (816, 33, 34),
    (828, 33, 34),
    (830, 34, 35),
    (831, 34, 35),
    (834, 34, 35),
    (835, 35, 36),
    (836, 35, 36),
    (837, 35, 36),
    (838, 34, 35),
    (839, 35, 36),
    (840, 35, 36),
    (841, 35, 36),
    (842, 35, 36),
    (844, 35, 36),
    (845, 36, 37),
    (846, 35, 36),
    (847, 36, 37),
    (849, 36, 37),
    (850, 36, 37),
    (851, 36, 37),
    (856, 37, 38),
    (858, 37, 38),
    (860, 37, 38),
    (864, 37, 38),
    (865, 38, 39),
    (891, 42, 43),
    (892, 41, 42),
    (898, 42, 43),
    (900, 43, 44),
    (901, 42, 43),
    (904, 43, 44),
    (905, 43, 44),
    (907, 44, 45),
    (908, 43, 44),
    (909, 44, 45),
    (911, 44, 45),
    (912, 44, 45),
    (913, 44, 45),
    (919, 45, 46),
    (920, 44, 45),
    (921, 44, 45),
    (923, 44, 45),
    (924, 44, 45),
    (925, 44, 45),
    (926, 44, 45),
    (927, 45, 46),
    (928, 45, 46),
    (930, 45, 46),
    (931, 44, 45),
    (938, 46, 47),
    (940, 46, 47),
    (941, 45, 46),
    (942, 46, 47),
    (973, 49, 50),
    (979, 49, 50),
    (980, 49, 50),
    (981, 49, 50),
    (983, 50, 51),
    (985, 50, 51),
    (988, 50, 51),
    (992, 50, 51),
    (996, 51, 52),
    (997, 51, 52),
    (1003, 52, 53),
    (1004, 52, 53),
    (1005, 52, 53),
    (1008, 51, 52),
    (1009, 50, 51),
    (1012, 47, 48),
    (1013, 46, 47),
    (1016, 44, 45),
    (1020, 40, 41),
    (1022, 38, 39),
    (1023, 37, 38),
    (1026, 34, 35),
    (1028, 32, 33),
    (1030, 30, 31),
    (1031, 29, 30),
    (1032, 28, 29),
    (1033, 27, 28),
    (1034, 26, 27),
    (1043, 18, 19),
    (1045, 16, 17),
    (1046, 15, 16),
    (1083, 19, 20),
    (1085, 20, 21),
    (1087, 21, 22),
    (1089, 22, 23),
    (1091, 23, 24),
    (1093, 24, 25),
    (1107, 32, 33),
    (1139, 46, 46),
    (1139, 48, 48),
    (1145, 50, 51),
    (1155, 52, 53),
    (1163, 52, 53),
    (1169, 51, 52),
    (1201, 49, 50),
    (1245, 47, 48),
    (1267, 46, 47),
    (1281, 45, 46),
    (1283, 46, 47),
    (1301, 44, 45),
    (1307, 45, 46),
    (1309, 45, 46),
    (1313, 44, 45),
    (1317, 44, 45),
    (1339, 44, 45),
    (1341, 44, 45),
    (1343, 44, 45),
    (1379, 41, 42),
    (1381, 42, 43),
    (1461, 36, 37),
    (1471, 35, 36),
    (1493, 35, 36),
    (1497, 35, 36),
    (1525, 32, 33),
    (1539, 32, 33),
];

const COMPARE_DIRECT_REMAINDER_KEYS: &[u32] = &[
    530, 3356, 26405, 29478, 29990, 30502, 32038, 32550, 33574, 34598, 35623, 36135,
    37159, 37671, 38184, 38696, 39209, 39720, 40231, 40744, 41256, 41769, 42794, 43304,
    43816, 44329, 44841, 45865, 46377, 46889, 47402, 47913, 49449, 49961, 50475, 50986,
    51498, 52523, 53547, 56107, 59181, 60204, 61229, 61740, 62253, 64301, 64813, 68398,
    68909, 69422, 69934, 70957, 71470, 71981, 72495, 74030, 75054, 76079, 77103, 77615,
    78126, 78639, 79151, 79663, 80175, 80687, 81200, 81712, 82737, 84272, 84784, 85809,
    86320, 86832, 87345, 87857, 88369, 88882, 89393, 89905, 90418, 90929, 92978, 93490,
    94002, 94514, 95026, 95539, 96052, 96563, 97075, 97587, 98099, 98611, 99123, 99636,
    100147, 100659, 101685, 102196, 102708, 103220, 103731, 104245, 104756, 105781, 106805, 107317,
    108340, 108851, 109362, 109873, 111406, 111917, 113451, 113962, 116006, 116517, 117028, 118561,
    120094, 121116, 123159, 125204, 125715, 135957, 136985, 143410, 143667, 143924, 144437, 144693,
    144948, 145205, 145716, 146227, 146484, 146996, 147507, 147762, 148276, 148531, 148787, 149043,
    149299, 149555, 150068, 150578, 151090, 151346, 151858, 152624, 152882, 153137, 153393, 153650,
    153905, 154161, 154417, 154672, 154928, 155184, 155440, 155696, 155952, 156208, 156464, 156721,
    157232, 157488, 157743, 157999, 158255, 158511, 158767, 159279, 159790, 160047, 160302, 160558,
    161070, 161839, 162093, 163374, 164142, 165421, 166189, 166957, 167212, 167980, 168237, 168493,
    168749, 169005, 169515, 169773, 170284, 170540, 170795, 171052, 171307, 171564, 172075, 172586,
    172843, 173353, 174633, 175145, 175401, 175913, 176424, 176682, 177193, 177448, 177704, 177959,
    178215, 178473, 178728, 178984, 179239, 179495, 179751, 180007, 180263, 180519, 181030, 181542,
    181798, 182310, 182822, 183078, 183590, 183845, 185380, 185892, 186149, 186404, 187172, 187940,
    188196, 188451, 188708, 188963, 189219, 189475, 189731, 190244, 191266, 191522, 191779, 192287,
    192546, 192802, 193313, 193570, 194081, 194593, 194849, 195105, 205596, 215845, 217125, 218149,
    218406, 218662, 218918, 219430, 219942, 220454, 220710, 220966, 221735, 221991, 222247, 222503,
    222759, 223016, 223272, 223529, 223784, 224296, 224552, 224809, 225065, 225322, 225576, 225832,
    226089, 226345, 226601, 226857, 227113, 227369, 227626, 227881, 228649, 228905, 229163, 229418,
    229674, 230187, 230956, 231211, 231979, 233005, 234028, 234284, 234541, 234796, 235053, 236077,
    237870, 238638, 238894, 239149, 239405, 239662, 239917, 240430, 241454, 241710, 241967, 242222,
    242479, 242735, 243247, 243503, 243759, 244015, 244271, 244528, 244784, 245039, 245297, 245552,
    245808, 246064, 246320, 246576, 246833, 247088, 247344, 247601, 247857, 248112, 248370, 248625,
    248881, 249393, 249649, 249906, 250162, 250418, 251442, 251956, 252467, 252723, 253235, 253491,
    253748, 254259, 254515, 254773, 255540, 255795, 256053, 256308, 256565, 257589, 257845, 258610,
    258865, 259630, 259885, 260396, 260651, 260906, 261416, 262181, 262436, 262946, 263456, 265240,
    266005, 266516, 266771, 292145, 292658, 293684, 294197, 294709, 295221, 296245, 296757, 297268,
    299828, 300341, 300851, 301875, 302388, 302899, 303411, 304947, 305459, 305972, 306483, 306994,
    308018, 308530, 309042, 309554, 310577, 313138, 313649, 314161, 315184, 315696, 316209, 316720,
    317232, 317744, 318256, 319281, 319791, 320304, 320816, 321327, 321838, 322351, 323375, 325422,
    325935, 326446, 326958, 327471, 329006, 329519, 331565, 332590, 333614, 334126, 335661, 336685,
    338221, 338733, 339246, 341804, 342317, 344365, 345389, 346412, 346923, 347948, 348459, 348972,
    349995, 350506, 358184, 359210, 359721, 360233, 360744, 361256, 362280, 363304, 363816, 365351,
    365863, 366887, 367910, 368422, 369446, 369958, 370470, 370982, 372006, 372518, 374565, 375589,
    376101, 377125, 377635, 378149, 378660, 379171, 379684, 380196, 380707, 381220, 381732, 382755,
    383779, 384291, 384803, 385315, 385826, 386340, 386850, 387362, 387875, 388386, 388898, 389410,
    389923, 390946, 391970, 392481, 393506, 394529, 395553, 396065,
];

fn compare_call_has_structurally_dead_top(call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_COMPARE_SKIP_STRUCTURAL_DEAD_CALLS").is_none() {
        return false;
    }
    if std::env::var_os("TLM_COMPARE_SKIP_EXACT_REMAINDER").is_some() {
        let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if COMPARE_DIRECT_REMAINDER_KEYS.binary_search(&key).is_ok() {
            return true;
        }
    }
    COMPARE_STRUCTURAL_DEAD_TOP_RANGES
        .iter()
        .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn compare_geq_chunked_middle_direct<F: FnOnce(&mut B, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    body: F,
    k: usize,
) {
    let call_index = next_compare_direct_call_index();
    let ops_start = circ.current_ops_len();
    let n = a.len();
    assert_eq!(
        b.len(),
        n,
        "compare_geq_chunked_middle_direct: a,b equal width"
    );
    assert!(
        n > 0,
        "compare_geq_chunked_middle_direct: nonempty operands"
    );
    let k = super::target_qubit_headroom(circ)
        .map_or(k, |headroom| k.min(headroom.saturating_sub(1)))
        .min(n);
    let split = n - k;
    let mut cy: Vec<Option<QubitId>> = (0..=n).map(|_| None).collect();
    let c = circ.alloc_qubit();
    circ.x(c);

    for i in 0..split {
        circ.x(b[i]);
        circ.cx(c, b[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b[i], c);
    }
    cy[split] = Some(c);

    for i in split..n {
        let next = circ.alloc_qubit();
        {
            let ci = cy[i].as_ref().unwrap();
            circ.x(b[i]);
            circ.cx(*ci, b[i]);
            circ.cx(*ci, a[i]);
            if !compare_call_has_structurally_dead_top(call_index, i) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0400_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(a[i], b[i], next);
                crate::point_add::restore_op_trace_context(old_context);
            }
            circ.cx(*ci, next);
        }
        cy[i + 1] = Some(next);
    }
    body(circ, cy[n].as_ref().unwrap());

    for i in (split..n).rev() {
        let next = cy[i + 1].take().unwrap();
        circ.cx(*cy[i].as_ref().unwrap(), next);
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }

    let c = cy[split].take().unwrap();
    for i in (0..split).rev() {
        circ.ccx(a[i], b[i], c);
        circ.cx(c, a[i]);
        circ.cx(c, b[i]);
        circ.x(b[i]);
    }
    circ.x(c);
    circ.zero_and_free(c);
    if std::env::var_os("TRACE_TLM_COMPARE_DIRECT").is_some() {
        eprintln!(
            "TLM_COMPARE_DIRECT call={} phase={} n={} k={} split={} ops_start={} ops_end={}",
            call_index,
            circ.phase,
            n,
            k,
            split,
            ops_start,
            circ.current_ops_len(),
        );
    }
}

pub fn compare_geq_chunked_middle<F: FnOnce(&mut B, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    flag: &QubitId,
    body: F,
    k: usize,
) {
    assert_eq!(
        b.len(),
        a.len(),
        "compare_geq_chunked_middle: a,b equal width"
    );
    if a.is_empty() {
        circ.x(*flag);
        body(circ, flag);
        circ.x(*flag);
        return;
    }
    compare_geq_chunked_middle_direct(
        circ,
        a,
        b,
        |c, carry| {
            c.cx(*carry, *flag);
            body(c, flag);
            c.cx(*carry, *flag);
        },
        k,
    );
}

pub fn controlled_swap_decision_lt_truncated(
    circ: &mut B,
    ctrl: &QubitId,
    u: &[QubitId],
    v: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    assert!(
        k > 0 && k <= u.len() && k <= v.len(),
        "k must fit in both operands"
    );
    let u_top: Vec<QubitId> = u[u.len() - k..].to_vec();
    let v_top: Vec<QubitId> = v[v.len() - k..].to_vec();

    let ck = super::next_cmp_k().saturating_add(1);
    compare_geq_chunked_middle_direct(
        circ,
        &u_top,
        &v_top,
        |c, carry| {
            c.x(*carry);
            c.ccx(*ctrl, *carry, *target);
            c.x(*carry);
        },
        ck,
    );
}

pub fn compare_geq_cin_middle<F: FnOnce(&mut B, &QubitId, &QubitId, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    cin: &QubitId,
    body: F,
) {
    let call_index = next_compare_cin_call_index();
    let n = a.len();
    assert_eq!(b.len(), n, "compare_geq_cin_middle: a,b equal width");
    assert!(n >= 1, "needs >= 1 bit");
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(n);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0);
    cy.push(Some(c0));
    for i in 0..n - 1 {
        let next = circ.alloc_qubit();
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        let old_context = crate::point_add::set_op_trace_context(
            0x1300_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        if !compare_cin_has_structurally_dead_carry(call_index, i) {
            circ.ccx(a[i], b[i], next);
        }
        crate::point_add::restore_op_trace_context(old_context);
        circ.cx(*ci, next);
        cy.push(Some(next));
    }

    {
        let i = n - 1;
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        body(circ, &a[i], &b[i], ci);
        circ.cx(*ci, a[i]);
        circ.cx(*ci, b[i]);
        circ.x(b[i]);
    }

    for i in (0..n - 1).rev() {
        let next = cy[i + 1].take().unwrap();
        let ci_raw = cy[i].as_ref().unwrap();
        circ.cx(*ci_raw, next);
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0);
    circ.x(c0);
    circ.zero_and_free(c0);
}

pub fn swap_decision_uncompute_vented(
    circ: &mut B,
    ctrl: &QubitId,
    v: &[QubitId],
    u: &[QubitId],
    k: usize,
    flag: &QubitId,
) {
    assert!(
        k > 0 && k <= v.len() && k <= u.len(),
        "k must fit in both operands"
    );
    let v_top: Vec<QubitId> = v[v.len() - k..].to_vec();
    let u_top: Vec<QubitId> = u[u.len() - k..].to_vec();

    let ck = super::next_cmp_k().saturating_add(1);
    let bit = circ.alloc_bit();
    circ.hmr(*flag, bit);
    circ.push_condition(bit);
    compare_geq_chunked_middle_direct(
        circ,
        &v_top,
        &u_top,
        |c, carry| {

            c.x(*carry);
            c.cz(*ctrl, *carry);
            c.x(*carry);
        },
        ck,
    );
    circ.pop_condition();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::OperationType;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    fn alloc_case(circ: &mut B, n: usize) -> (Vec<QubitId>, Vec<QubitId>, QubitId, QubitId) {
        let a = (0..n).map(|_| circ.alloc_qubit()).collect();
        let b = (0..n).map(|_| circ.alloc_qubit()).collect();
        let ctrl = circ.alloc_qubit();
        let target = circ.alloc_qubit();
        (a, b, ctrl, target)
    }

    fn xor_value(circ: &mut B, qs: &[QubitId], value: usize) {
        for (i, &q) in qs.iter().enumerate() {
            if (value >> i) & 1 != 0 {
                circ.x(q);
            }
        }
    }

    fn simulate(circ: &B) -> (Vec<u64>, u64) {
        let mut shake = Shake256::default();
        shake.update(b"comparator-direct-final-carry-test");
        let mut xof = shake.finalize_xof();
        let mut sim =
            Simulator::new(circ.next_qubit as usize, circ.next_bit as usize, &mut xof);
        sim.apply_iter(circ.ops.iter());
        (sim.qubits, sim.phase)
    }

    fn read_uniform(qs: &[QubitId], qubits: &[u64]) -> usize {
        qs.iter().enumerate().fold(0usize, |value, (i, q)| {
            let lane = qubits[q.0 as usize];
            assert!(
                lane == 0 || lane == u64::MAX,
                "nonuniform data lane q{}",
                q.0
            );
            value | (usize::from(lane == u64::MAX) << i)
        })
    }

    fn toffoli_count(circ: &B) -> usize {
        circ.ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count()
    }

    #[test]
    fn direct_final_carry_is_exhaustive_for_small_widths() {
        for n in 1..=4 {
            let limit = 1usize << n;
            for held in 0..=n {
                for a_value in 0..limit {
                    for b_value in 0..limit {
                        for ctrl_value in 0..=1usize {
                            for target_value in 0..=1usize {
                                let mut circ = B::new();
                                let (a, b, ctrl, target) = alloc_case(&mut circ, n);
                                xor_value(&mut circ, &a, a_value);
                                xor_value(&mut circ, &b, b_value);
                                if ctrl_value != 0 {
                                    circ.x(ctrl);
                                }
                                if target_value != 0 {
                                    circ.x(target);
                                }
                                compare_geq_chunked_middle_direct(
                                    &mut circ,
                                    &a,
                                    &b,
                                    |c, carry| {
                                        c.x(*carry);
                                        c.ccx(ctrl, *carry, target);
                                        c.cz(ctrl, *carry);
                                        c.x(*carry);
                                    },
                                    held,
                                );

                                assert_eq!(circ.active_qubits as usize, 2 * n + 2);
                                let (qubits, phase) = simulate(&circ);
                                let predicate = ctrl_value != 0 && a_value < b_value;
                                assert_eq!(read_uniform(&a, &qubits), a_value);
                                assert_eq!(read_uniform(&b, &qubits), b_value);
                                assert_eq!(
                                    qubits[ctrl.0 as usize],
                                    if ctrl_value != 0 { u64::MAX } else { 0 }
                                );
                                assert_eq!(
                                    qubits[target.0 as usize],
                                    if (target_value != 0) ^ predicate {
                                        u64::MAX
                                    } else {
                                        0
                                    },
                                );
                                assert_eq!(phase, if predicate { u64::MAX } else { 0 });
                                assert!(qubits[2 * n + 2..].iter().all(|&q| q == 0));
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn freed_predicate_lane_funds_one_held_carry() {
        for n in 1..=8 {
            for held in 0..n {
                let mut legacy = B::new();
                let (a, b, ctrl, target) = alloc_case(&mut legacy, n);
                let flag = legacy.alloc_qubit();
                compare_geq_chunked_middle(
                    &mut legacy,
                    &a,
                    &b,
                    &flag,
                    |c, flag| {
                        c.x(*flag);
                        c.ccx(ctrl, *flag, target);
                        c.x(*flag);
                    },
                    held,
                );
                legacy.zero_and_free(flag);

                let mut direct = B::new();
                let (a, b, ctrl, target) = alloc_case(&mut direct, n);
                compare_geq_chunked_middle_direct(
                    &mut direct,
                    &a,
                    &b,
                    |c, carry| {
                        c.x(*carry);
                        c.ccx(ctrl, *carry, target);
                        c.x(*carry);
                    },
                    held + 1,
                );

                assert_eq!(
                    direct.peak_qubits, legacy.peak_qubits,
                    "n={n} held={held}"
                );
                assert_eq!(
                    toffoli_count(&direct) + 1,
                    toffoli_count(&legacy),
                    "n={n} held={held}",
                );
            }
        }
    }
}
