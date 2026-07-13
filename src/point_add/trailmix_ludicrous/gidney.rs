use super::comparator::compare_geq_cin_middle;
use super::{BExt, B};
use crate::circuit::QubitId;
use std::cell::{Cell, RefCell};

thread_local! {
    static DIRTY_VENT_POOL: RefCell<Vec<QubitId>> = const { RefCell::new(Vec::new()) };
    static THREADED_ADD_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static HYBRID_ADD_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static ERASE_GATED_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static ERASE_GATED_CAPPED_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_gidney_call_index() {
    THREADED_ADD_CALL_INDEX.with(|index| index.set(0));
    HYBRID_ADD_CALL_INDEX.with(|index| index.set(0));
    ERASE_GATED_CALL_INDEX.with(|index| index.set(0));
    ERASE_GATED_CAPPED_CALL_INDEX.with(|index| index.set(0));
}

fn next_threaded_add_call_index() -> usize {
    THREADED_ADD_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_hybrid_add_call_index() -> usize {
    HYBRID_ADD_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_erase_gated_call_index() -> usize {
    ERASE_GATED_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_erase_gated_capped_call_index() -> usize {
    ERASE_GATED_CAPPED_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

const GIDNEY_THREAD_FWD_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (2591, 30, 30),
    (2591, 52, 52),
    (2591, 54, 253),
    (5, 0, 123),
    (1, 18, 18),
    (1, 20, 121),
    (0, 53, 122),
    (4, 0, 65),
    (3, 21, 21),
    (3, 23, 65),
    (2850, 0, 14),
    (2851, 0, 13),
    (2852, 0, 12),
    (2853, 0, 11),
    (2854, 0, 10),
    (2334, 0, 9),
    (2855, 0, 9),
    (2865, 0, 9),
    (2335, 3, 11),
    (2856, 0, 8),
    (2878, 3, 11),
    (2337, 6, 13),
    (2857, 0, 7),
    (2909, 6, 13),
    (2336, 5, 11),
    (2858, 0, 6),
    (2894, 5, 11),
    (2338, 8, 13),
    (2859, 0, 5),
    (2925, 8, 13),
    (2303, 7, 11),
    (2339, 9, 13),
    (2340, 10, 14),
    (2860, 0, 4),
    (2944, 9, 13),
    (2964, 10, 14),
    (2272, 10, 13),
    (2341, 11, 14),
    (2344, 14, 17),
    (2602, 249, 252),
    (2844, 9, 9),
    (2844, 11, 13),
    (2846, 8, 11),
    (2861, 0, 3),
    (2985, 11, 14),
    (3007, 12, 15),
    (62, 52, 55),
    (68, 53, 56),
    (74, 54, 57),
    (80, 51, 54),
    (86, 52, 55),
    (92, 53, 56),
    (98, 50, 53),
    (104, 52, 54),
    (110, 53, 55),
    (116, 50, 52),
    (122, 51, 53),
    (128, 52, 54),
    (134, 49, 51),
    (140, 50, 52),
    (146, 51, 53),
    (152, 48, 50),
    (158, 49, 51),
    (164, 50, 52),
    (170, 47, 49),
    (182, 49, 51),
    (2256, 11, 13),
    (2287, 9, 11),
    (2316, 7, 9),
    (2332, 6, 8),
    (2342, 13, 15),
    (2343, 14, 16),
    (2358, 29, 31),
    (2601, 251, 253),
    (2603, 249, 251),
    (2604, 248, 250),
    (2605, 247, 249),
    (2606, 246, 248),
    (2607, 245, 247),
    (2613, 239, 241),
    (2843, 11, 13),
    (2845, 9, 11),
    (2847, 7, 9),
    (2848, 6, 8),
    (2862, 0, 2),
    (3030, 14, 16),
    (3052, 15, 17),
    (56, 51, 53),
];

const GIDNEY_THREAD_SUM_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (2334, 0, 10),
    (2865, 0, 10),
    (2335, 3, 12),
    (2878, 3, 12),
    (2337, 6, 14),
    (2909, 6, 14),
    (2336, 5, 12),
    (2894, 5, 12),
    (2338, 8, 14),
    (2925, 8, 14),
    (2339, 9, 14),
    (2340, 10, 15),
    (2844, 9, 14),
    (2944, 9, 14),
    (2964, 10, 15),
    (2256, 10, 14),
    (2272, 10, 14),
    (2332, 5, 9),
    (2341, 11, 15),
    (2344, 14, 18),
    (2602, 249, 253),
    (2846, 8, 12),
    (2985, 11, 15),
    (3007, 12, 16),
    (128, 51, 54),
    (2287, 9, 12),
    (2303, 9, 12),
    (2316, 7, 10),
    (2342, 13, 16),
    (2343, 14, 17),
    (2358, 29, 32),
    (2601, 251, 254),
    (2604, 248, 251),
    (2605, 247, 250),
    (2606, 246, 249),
    (2607, 245, 248),
    (2613, 239, 242),
    (2843, 11, 14),
    (2845, 9, 12),
    (2847, 7, 10),
    (2848, 6, 9),
    (2849, 6, 9),
    (3030, 14, 17),
    (3052, 15, 18),
    (3185, 22, 25),
    (62, 52, 55),
    (68, 53, 56),
    (74, 54, 57),
    (80, 51, 54),
    (86, 52, 55),
    (92, 53, 56),
    (98, 50, 53),
    (104, 52, 54),
    (110, 53, 55),
    (116, 50, 52),
    (122, 51, 53),
    (140, 50, 52),
    (146, 51, 53),
    (152, 48, 50),
    (158, 49, 51),
    (164, 50, 52),
    (170, 47, 49),
    (182, 49, 51),
    (2217, 13, 15),
    (2237, 12, 14),
    (2345, 16, 18),
    (2346, 17, 19),
    (2354, 26, 28),
    (2355, 27, 29),
    (2356, 28, 30),
    (2357, 29, 31),
    (2359, 31, 33),
    (2360, 32, 34),
    (2361, 33, 35),
    (2362, 34, 36),
    (2363, 35, 37),
    (2364, 36, 38),
    (2365, 37, 39),
    (2368, 40, 42),
    (2370, 42, 44),
    (2599, 252, 254),
    (2600, 252, 254),
    (2608, 245, 247),
    (2609, 244, 246),
    (2611, 242, 244),
    (2612, 241, 243),
    (2614, 239, 241),
    (2615, 238, 240),
    (2617, 236, 238),
    (2618, 235, 237),
    (2620, 233, 235),
    (2621, 232, 234),
    (2626, 227, 229),
    (2627, 226, 228),
    (2838, 15, 17),
    (2839, 14, 16),
    (2841, 13, 15),
    (2842, 12, 14),
    (3069, 16, 18),
    (3084, 17, 19),
    (3098, 17, 19),
    (3201, 24, 26),
    (3216, 25, 27),
    (3232, 24, 26),
    (3247, 25, 27),
    (3261, 26, 28),
    (3290, 26, 28),
    (3304, 27, 29),
    (3318, 26, 28),
    (3332, 27, 29),
    (3350, 28, 30),
    (3364, 27, 29),
    (3397, 29, 31),
    (3415, 28, 30),
    (3430, 29, 31),
    (3444, 30, 32),
    (3459, 29, 31),
    (56, 51, 53),
];

const GIDNEY_THREAD_BOUNDARY_DEAD_CALLS: &[usize] = &[
    0, 1, 1002, 1010, 1018, 1026, 1034, 104, 1043, 1060, 1069, 1078, 1087, 1096, 110, 1105, 1114,
    1123, 1132, 1141, 1150, 116, 1168, 1177, 1186, 1195, 1204, 1213, 122, 1222, 1231, 1240, 1249,
    1258, 1259, 1268, 1278, 128, 1287, 1296, 1297, 1306, 1307, 1325, 1334, 1335, 1344, 1353, 1362,
    1371, 1380, 1389, 1399, 14, 140, 1408, 1418, 1428, 1438, 1448, 1458, 146, 1468, 1478, 1488,
    1498, 1508, 152, 158, 164, 170, 176, 182, 188, 194, 2, 20, 200, 206, 212, 218, 224, 230, 236,
    242, 248, 254, 26, 260, 266, 272, 278, 284, 2850, 2851, 2852, 2853, 2854, 2855, 2856, 2857,
    2858, 2859, 2860, 2861, 2862, 2863, 2864, 290, 296, 302, 308, 314, 32, 320, 326, 333, 339, 346,
    353, 360, 367, 3674, 3684, 3694, 3704, 3714, 3724, 3734, 374, 3744, 3754, 3764, 3774, 3783,
    3793, 38, 3802, 381, 3811, 3820, 3829, 3838, 3847, 3848, 3857, 3866, 3876, 388, 3886, 3895,
    3905, 3915, 3924, 3934, 3943, 395, 3952, 3961, 3970, 3979, 3988, 3997, 4006, 4015, 402, 4024,
    4033, 4042, 4051, 4060, 4069, 4078, 4087, 409, 4096, 4105, 4114, 4123, 4132, 4140, 4149, 4157,
    416, 4165, 4173, 4181, 4189, 4197, 4205, 4213, 4221, 4229, 423, 4237, 4245, 4253, 4261, 4269,
    4277, 4285, 4293, 430, 4301, 4309, 4317, 4325, 4333, 4341, 4349, 4357, 4365, 437, 4373, 4389,
    4397, 44, 4404, 4412, 4420, 4428, 4436, 444, 4444, 4451, 4459, 4466, 4480, 4487, 4494, 4501,
    4508, 451, 4515, 4522, 4529, 4543, 4550, 4557, 4564, 4571, 4578, 458, 4585, 4592, 4599, 4606,
    4613, 4620, 4627, 4634, 4641, 4648, 465, 4655, 4662, 4669, 4676, 4683, 4690, 4697, 4704, 4711,
    4718, 472, 4725, 4732, 4739, 4746, 4753, 4760, 4767, 4774, 4781, 4788, 479, 4795, 4802, 4809,
    4816, 4823, 4830, 4837, 4844, 4850, 4857, 486, 4863, 4869, 4875, 4881, 4887, 4893, 4899, 4905,
    4911, 4917, 4923, 4929, 493, 4935, 4941, 4947, 4953, 4959, 4965, 4971, 4977, 4983, 4989, 4995,
    50, 500, 5001, 5007, 5013, 5019, 5025, 5031, 5037, 5043, 5049, 5055, 5061, 5067, 507, 5073,
    5079, 5085, 5091, 5097, 5103, 5109, 5115, 5121, 514, 521, 528, 535, 542, 549, 556, 56, 563,
    570, 577, 584, 591, 598, 605, 612, 619, 62, 626, 633, 640, 647, 654, 661, 668, 675, 68, 682,
    689, 696, 703, 710, 717, 724, 732, 739, 74, 747, 755, 763, 771, 779, 786, 794, 8, 80, 802, 810,
    818, 826, 834, 842, 850, 86, 866, 874, 882, 890, 898, 906, 914, 92, 922, 930, 938, 946, 954,
    962, 970, 978, 98, 986, 994,
];

const GIDNEY_THREAD_FWD_REMAINDER_KEYS: &[u32] = &[
    521, 9767, 11304, 11305, 12846, 12847, 45105, 45106, 48175, 48176, 49712, 49713, 51249, 51250,
    52782, 52783, 54319, 54320, 55856, 55857, 57389, 57390, 58926, 58927, 60463, 60464, 61996,
    61997, 63533, 63534, 65070, 65071, 66603, 66604, 68140, 68141, 69677, 69678, 71210, 71211,
    72747, 72748, 74284, 74285, 75817, 75818, 77354, 77355, 78892, 80424, 80425, 81961, 81962,
    83498, 83499, 85287, 85288, 86824, 86825, 88618, 90407, 92200, 93992, 93993, 95781, 95782,
    97575, 99368, 101156, 101157, 102950, 104742, 104743, 106532, 108324, 108325, 110118, 111907,
    113699, 113700, 115493, 117282, 119075, 120868, 122657, 124450, 126243, 128032, 129825, 131618,
    133406, 133407, 135200, 136993, 138782, 140575, 142368, 144156, 144157, 145950, 147743, 149532,
    151325, 153118, 154911, 154912, 156704, 156705, 158493, 160282, 162075, 163868, 165657, 167454,
    167455, 169243, 171036, 171037, 172830, 174622, 174623, 176406, 176407, 178204, 178205, 179997,
    179998, 181782, 183579, 183580, 185372, 185373, 187417, 187418, 189210, 189211, 191259, 191260,
    193305, 195353, 195354, 197398, 201236, 203285, 205330, 207383, 207384, 209432, 209433, 211477,
    211478, 213522, 215576, 217620, 217621, 219669, 219670, 221718, 221719, 223763, 223764, 225812,
    225813, 227861, 227862, 229907, 231956, 234000, 236045, 240143, 242188, 244237, 246291, 248335,
    248336, 250384, 250385, 252429, 254478, 254479, 256523, 260617, 262665, 262666, 264715, 267016,
    269065, 271370, 273671, 275976, 278281, 280582, 282887, 285192, 287493, 289798, 294404, 296709,
    299014, 301314, 301315, 303620, 305925, 308226, 310531, 312836, 319747, 324609, 327170, 331817,
    334593, 336936, 339240, 344102, 346406, 348710, 351012, 353316, 355620, 358178, 360482, 363042,
    365600, 368160, 370720, 373278, 375837, 375838, 378398, 380956, 383516, 386076, 397092, 402722,
    411681, 427296, 430623, 441118, 444959, 456990, 460829, 465436, 480796, 484379, 491803, 495386,
    503322, 511256, 545041, 550672, 556559, 562190, 567565, 567566, 572684, 572685, 600336, 600337,
    600593, 600594, 600850, 601364, 601621, 601878, 602394, 602650, 602651, 602907, 602908, 603164,
    603165, 603421, 603422, 603935, 603936, 604192, 604193, 604449, 604450, 604706, 604707, 604963,
    604964, 605220, 605221, 605478, 605735, 605992, 606249, 606506, 608048, 626553, 638889, 665341,
    665596, 665597, 665852, 665853, 667893, 667894, 668148, 668149, 668404, 668658, 668659, 668913,
    668914, 669423, 669424, 669678, 669679, 669933, 669934, 670188, 670189, 670443, 670444, 670698,
    670699, 670953, 670954, 671208, 671209, 671464, 671718, 671719, 671974, 672229, 672483, 672484,
    672738, 672739, 672994, 673249, 673504, 673759, 674014, 674269, 674524, 674779, 675034, 675544,
    675799, 676054, 676309, 676819, 677074, 677329, 678604, 679114, 679624, 680644, 680899, 687530,
    687785, 689315, 689825, 690080, 690845, 691100, 691610, 691865, 692120, 692630, 692885, 693140,
    694670, 694925, 695180, 695690, 695945, 696200, 696455, 696710, 696965, 697475, 697730, 699515,
    700025, 700279, 700535, 707419, 708439, 709969, 710224, 710479, 710989, 711499, 712264, 713284,
    714304, 714559, 716089, 716344, 716599, 717364, 717619, 718384, 718639, 719404, 720424, 720934,
    721189, 721954, 722464, 722974, 723229, 723484, 723739, 723994, 726289, 726543, 726544, 726799,
    727054, 727309, 727310, 727564, 727565, 785680, 785681, 789521, 789522, 793105, 793106, 796435,
    799764, 803093, 806422, 810775, 815383, 815384, 819480, 819481, 823321, 823322, 827416, 827417,
    831257, 831258, 834842, 834843, 838682, 842266, 842267, 845851, 845852, 849434, 849435, 853019,
    853020, 857628, 857629, 861211, 861212, 865821, 869661, 869662, 874268, 874269, 878109, 878110,
    881694, 885533, 885534, 889119, 892447, 892448, 896030, 896031, 955936, 958496, 961056, 963618,
    973348, 975652, 985088, 992257, 1020933, 1023236, 1032452, 1043975, 1046278, 1055498, 1059848,
    1070347, 1072399, 1074445, 1076496, 1076497, 1078544, 1080595, 1092880, 1094932, 1096979,
    1099030, 1101077, 1103124, 1105175, 1107222, 1109269, 1115414, 1117465, 1131546, 1133593,
    1135643, 1135644, 1137691, 1139482, 1143324, 1146910, 1152287, 1154078, 1155869, 1159455,
    1170209, 1172000, 1195298, 1198880, 1222183, 1223974, 1227560, 1229351, 1248044, 1251114,
    1255723, 1257262, 1260332, 1261871, 1263406, 1272624, 1274159, 1275698, 1278768, 1280307,
    1284916, 1294134, 1298743, 1301813, 1303352,
];

const GIDNEY_THREAD_SUM_REMAINDER_KEYS: &[u32] = &[
    9767, 11304, 11305, 12846, 12847, 45105, 45106, 48175, 48176, 49712, 49713, 51249, 51250,
    52782, 52783, 54319, 54320, 55856, 55857, 57389, 57390, 58926, 58927, 60463, 60464, 61996,
    61997, 63533, 63534, 65070, 65071, 66604, 68140, 68141, 69677, 69678, 71210, 71211, 72747,
    72748, 74284, 74285, 75817, 75818, 77354, 77355, 78892, 80424, 80425, 81961, 81962, 83498,
    83499, 85287, 85288, 86824, 86825, 88617, 88618, 90407, 92200, 93993, 95781, 95782, 97575,
    99368, 101156, 101157, 102950, 104742, 104743, 106532, 108324, 108325, 110118, 111907, 113699,
    113700, 117282, 120868, 122657, 124450, 126243, 128032, 129825, 131618, 133406, 133407, 135200,
    136993, 138782, 140575, 142368, 144156, 144157, 145950, 147743, 149532, 151325, 153118, 154911,
    154912, 156704, 156705, 158493, 160282, 162075, 163868, 165657, 167454, 167455, 169243, 171036,
    171037, 172830, 174622, 174623, 176406, 176407, 178204, 178205, 179997, 179998, 181782, 183579,
    183580, 185372, 185373, 187417, 187418, 189210, 189211, 191259, 191260, 193305, 195353, 195354,
    197398, 201236, 203285, 205330, 207383, 207384, 209433, 211477, 211478, 213522, 215576, 217620,
    217621, 221718, 221719, 223763, 223764, 225812, 225813, 227861, 227862, 229907, 231956, 234000,
    236045, 240143, 242188, 246291, 248336, 250384, 250385, 252429, 254478, 254479, 256523, 260617,
    262666, 264715, 267016, 271370, 273671, 275975, 275976, 278281, 280582, 282887, 285192, 287492,
    287493, 289798, 294403, 294404, 299014, 301314, 301315, 305925, 308226, 310531, 312836, 315137,
    317442, 319747, 322304, 324609, 327170, 331817, 332032, 334592, 334593, 339240, 341760, 344102,
    346406, 348710, 351012, 353316, 355620, 358178, 360482, 363042, 365600, 368160, 370720, 373278,
    375837, 375838, 378398, 380956, 383516, 386076, 388902, 391461, 394276, 397093, 399907, 399908,
    402722, 405796, 408610, 408611, 411682, 414755, 441118, 441119, 444959, 444960, 448543, 452382,
    456990, 456991, 460829, 460830, 465436, 465437, 469021, 469022, 473629, 477212, 480796, 480797,
    484379, 484380, 487963, 491803, 491804, 495386, 495387, 499226, 503322, 503323, 507162, 511256,
    511257, 523286, 536851, 540690, 545041, 545042, 550672, 550673, 556559, 556560, 562190, 562191,
    600850, 600851, 601108, 601364, 601365, 601621, 601622, 601878, 601879, 602136, 602394, 602395,
    605735, 605736, 605992, 605993, 606506, 606507, 607020, 607021, 608820, 609077, 609333, 609334,
    610362, 610619, 610876, 611133, 612161, 612418, 612675, 613703, 613960, 614217, 614474, 614731,
    615245, 615759, 616273, 616530, 616787, 617044, 618072, 618329, 618586, 619100, 619614, 620128,
    620642, 620899, 621413, 621670, 621927, 622184, 622441, 622955, 623212, 623469, 623983, 625268,
    625525, 625782, 626039, 626040, 626553, 626554, 627068, 628867, 629123, 629124, 629381, 629638,
    629895, 630152, 630409, 630666, 630923, 631437, 631693, 631694, 633493, 633750, 634007, 634264,
    634521, 635549, 635806, 636320, 636834, 637348, 638118, 638889, 638890, 639147, 639917, 640431,
    640688, 641202, 641716, 642487, 642744, 643001, 643258, 643772, 644029, 644286, 644543, 645057,
    645314, 645828, 646856, 647113, 647627, 648141, 648398, 648912, 649169, 649426, 649939, 649940,
    650197, 650454, 650711, 651225, 651482, 653281, 653538, 653795, 654052, 654823, 655594, 656622,
    656879, 657907, 660477, 665086, 665341, 665342, 668404, 668405, 670699, 670700, 671464, 671465,
    671719, 671720, 671974, 671975, 672229, 672230, 672994, 672995, 673249, 673250, 673504, 673505,
    673759, 673760, 674014, 674015, 674269, 674270, 674524, 674525, 674779, 674780, 675034, 675035,
    675290, 675544, 675545, 675799, 675800, 676054, 676055, 676309, 676310, 676565, 676819, 676820,
    677074, 677075, 677329, 677330, 677585, 677840, 678095, 678350, 678604, 678605, 678860, 679114,
    679115, 679370, 679624, 679625, 679880, 680135, 680390, 680900, 681155, 681410, 681665, 681920,
    682175, 682430, 682685, 682940, 683195, 683705, 683960, 684214, 684215, 684470, 684980, 685235,
    685490, 685745, 686000, 686255, 686510, 686765, 687020, 687275, 687530, 687531, 687785, 687786,
    688040, 688295, 688550, 688805, 689060, 689315, 689316, 689570, 689825, 689826, 690080, 690081,
    690336, 690590, 690845, 690846, 691100, 691101, 691355, 691610, 691611, 691865, 691866, 692120,
    692121, 692376, 692630, 692631, 692885, 692886, 693140, 693141, 693395, 693650, 693905, 694160,
    694415, 694671, 694925, 694926, 695180, 695181, 695435, 695690, 695691, 695945, 695946, 696200,
    696201, 696455, 696456, 696710, 696711, 696965, 696966, 697475, 697476, 697730, 697984, 697985,
    698495, 698750, 699005, 699260, 699515, 699516, 699770, 699771, 700025, 700026, 700279, 700280,
    700535, 700790, 701045, 701300, 701555, 701810, 702065, 702575, 702830, 703085, 703340, 703595,
    703850, 704105, 704360, 704615, 704869, 704870, 705125, 705380, 705635, 706145, 706400, 706655,
    706910, 707165, 707675, 707930, 708185, 708439, 708440, 708950, 709205, 709460, 709715, 710224,
    710225, 710479, 710480, 710735, 710989, 710990, 711245, 711499, 711500, 711755, 712010, 712264,
    712265, 712520, 712775, 713030, 713285, 713794, 713795, 714050, 714305, 714560, 714815, 715070,
    715325, 715580, 716089, 716090, 716344, 716345, 716599, 716600, 716855, 717109, 717110, 717364,
    717365, 717619, 717620, 717875, 719404, 719405, 719660, 719915, 720170, 720424, 720425, 720680,
    720934, 720935, 721189, 721190, 721445, 721700, 721954, 721955, 722210, 722464, 722465, 722720,
    722974, 722975, 723229, 723230, 723484, 723485, 723739, 723740, 723994, 723995, 724503, 724758,
    725013, 725268, 725779, 726034, 726289, 726290, 727054, 727055, 796435, 796436, 799764, 799765,
    803093, 803094, 806422, 806423, 810775, 810776, 838682, 838683, 865821, 865822, 905761, 908834,
    911907, 914978, 918051, 923939, 926756, 929573, 932388, 935205, 950814, 955936, 958496, 961056,
    963618, 973348, 975652, 985088, 987432, 992257, 997162, 1014020, 1020933, 1023236, 1032452,
    1041672, 1043975, 1046278, 1055498, 1059848, 1064202, 1068300, 1070347, 1072399, 1074445,
    1076496, 1076497, 1078544, 1080595, 1092880, 1094932, 1096979, 1099030, 1101076, 1101077,
    1103124, 1105175, 1107221, 1107222, 1109269, 1115414, 1117465, 1131546, 1133592, 1133593,
    1135643, 1135644, 1137690, 1137691, 1139482, 1143324, 1146910, 1152287, 1154078, 1155869,
    1159455, 1170209, 1171999, 1172000, 1193503, 1195298, 1197089, 1198880, 1222183, 1223974,
    1227560, 1229351, 1244970, 1248044, 1251114, 1252653, 1255723, 1257262, 1260332, 1261871,
    1263406, 1272624, 1274159, 1275698, 1278768, 1280307, 1284916, 1292595, 1294134, 1298743,
    1301813, 1303352, 1304887,
];

fn gidney_structural_dead_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_STRUCTURAL_DEAD_CALLS").is_some()
}

fn gidney_skip_top2_thread_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_TOP2_THREAD").is_some()
}

fn gidney_skip_fullvent_top2_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_FULLVENT_TOP2").is_some()
}

fn gidney_skip_exact_remainder_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_EXACT_REMAINDER").is_some()
}

fn gidney_skip_exact_fwd_remainder_enabled() -> bool {
    gidney_skip_exact_remainder_enabled()
        || std::env::var_os("TLM_GIDNEY_SKIP_EXACT_FWD_REMAINDER").is_some()
}

fn gidney_skip_exact_sum_remainder_enabled() -> bool {
    gidney_skip_exact_remainder_enabled()
        || std::env::var_os("TLM_GIDNEY_SKIP_EXACT_SUM_REMAINDER").is_some()
}

fn gidney_skip_exact_erase_ccz_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_EXACT_ERASE_CCZ").is_some()
        || std::env::var_os("TLM_GIDNEY_SKIP_EXACT_ERASE_ALL_CCZ").is_some()
}

fn gidney_skip_exact_erase_capped_ccz_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_EXACT_ERASE_ALL_CCZ").is_some()
        || std::env::var_os("TLM_GIDNEY_SKIP_EXACT_ERASE_CAPPED_CCZ").is_some()
}

fn gidney_skip_small_residual_enabled() -> bool {
    std::env::var_os("TLM_GIDNEY_SKIP_SMALL_RESIDUAL_DEAD").is_some()
}

fn gidney_key(call_index: usize, bit: usize) -> u32 {
    (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff)
}

const GIDNEY_THREAD_BOUNDARY_RESIDUAL_CALLS: &[usize] = &[
    134, 858, 1051, 1159, 1316, 3885, 3923, 4381, 4473, 4536, 5127,
];

const GIDNEY_ERASE_CCZ_RESIDUAL_CALLS: &[usize] =
    &[2418, 2425, 2431, 2444, 2522, 2552, 2591, 2641, 2646, 2727];

const GIDNEY_ERASE_CAPPED_CCZ_RESIDUAL_CALLS: &[usize] = &[
    15, 595, 605, 607, 626, 801, 807, 834, 849, 867, 870, 876, 888, 897, 903, 909, 933, 939, 942,
    957, 960, 975, 996, 999, 1002, 1035, 1047, 1056, 1071, 1089, 1095, 1107, 1116, 1128, 1137,
];

const GIDNEY_ERASE_CCZ_REMAINDER_CALLS: &[usize] = &[
    0, 1, 2, 308, 320, 337, 370, 379, 384, 389, 394, 399, 404, 409, 414, 419, 424, 429, 434, 439,
    444, 449, 454, 460, 465, 471, 477, 483, 489, 495, 501, 507, 513, 519, 531, 537, 543, 549, 555,
    561, 567, 573, 579, 585, 591, 597, 604, 610, 617, 630, 637, 656, 1518, 1519, 1520, 1521, 1522,
    1523, 1524, 1525, 1526, 1527, 1528, 1529, 1530, 1531, 1532, 2379, 2398, 2405, 2450, 2456, 2462,
    2468, 2474, 2480, 2498, 2504, 2510, 2516, 2528, 2534, 2540, 2558, 2564, 2575, 2581, 2586, 2596,
    2601, 2606, 2611, 2616, 2621, 2626, 2651, 2656, 2665,
];

const GIDNEY_ERASE_CAPPED_CCZ_REMAINDER_CALLS: &[usize] = &[
    1, 2, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81,
    84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141,
    144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198,
    201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255,
    258, 261, 264, 267, 270, 273, 276, 279, 282, 285, 288, 291, 294, 297, 300, 303, 306, 309, 312,
    315, 318, 321, 324, 327, 330, 333, 336, 339, 342, 345, 348, 351, 357, 360, 366, 369, 372, 378,
    381, 384, 387, 390, 393, 396, 402, 537, 542, 544, 547, 549, 551, 553, 555, 557, 559, 561, 563,
    565, 567, 569, 571, 573, 575, 577, 579, 581, 583, 585, 587, 589, 591, 593, 597, 599, 601, 603,
    609, 611, 613, 615, 617, 619, 624, 634, 765, 771, 774, 777, 780, 783, 786, 789, 795, 798, 816,
    819, 822, 825, 828, 831, 837, 843, 846, 852, 855, 858, 864, 879, 882, 885, 891, 900, 906, 912,
    915, 921, 924, 927, 930, 936, 945, 951, 954, 963, 969, 972, 978, 981, 984, 987, 990, 993, 1005,
    1008, 1011, 1014, 1017, 1020, 1023, 1026, 1029, 1032, 1038, 1041, 1044, 1050, 1053, 1059, 1062,
    1065, 1068, 1074, 1077, 1080, 1083, 1086, 1092, 1098, 1101, 1104, 1110, 1113, 1119, 1122, 1125,
    1131,
];

fn gidney_erase_ccz_has_exact_dead_call(call_index: usize) -> bool {
    if gidney_skip_small_residual_enabled()
        && GIDNEY_ERASE_CCZ_RESIDUAL_CALLS
            .binary_search(&call_index)
            .is_ok()
    {
        return true;
    }
    gidney_skip_exact_erase_ccz_enabled()
        && GIDNEY_ERASE_CCZ_REMAINDER_CALLS
            .binary_search(&call_index)
            .is_ok()
}

fn gidney_erase_capped_ccz_has_exact_dead_call(call_index: usize) -> bool {
    if gidney_skip_small_residual_enabled()
        && GIDNEY_ERASE_CAPPED_CCZ_RESIDUAL_CALLS
            .binary_search(&call_index)
            .is_ok()
    {
        return true;
    }
    gidney_skip_exact_erase_capped_ccz_enabled()
        && GIDNEY_ERASE_CAPPED_CCZ_REMAINDER_CALLS
            .binary_search(&call_index)
            .is_ok()
}

fn threaded_add_call_has_structurally_dead_forward(
    call_index: usize,
    bit: usize,
    total: usize,
    width: usize,
    vents: usize,
) -> bool {
    if gidney_skip_exact_fwd_remainder_enabled()
        && GIDNEY_THREAD_FWD_REMAINDER_KEYS
            .binary_search(&gidney_key(call_index, bit))
            .is_ok()
    {
        return true;
    }
    if gidney_skip_fullvent_top2_enabled() && vents >= width && bit + 2 >= width {
        return true;
    }
    if gidney_skip_top2_thread_enabled() && bit + 2 >= total {
        return true;
    }
    gidney_structural_dead_enabled()
        && GIDNEY_THREAD_FWD_DEAD_RANGES
            .iter()
            .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

fn threaded_add_call_has_structurally_dead_boundary(call_index: usize) -> bool {
    if gidney_skip_small_residual_enabled()
        && GIDNEY_THREAD_BOUNDARY_RESIDUAL_CALLS
            .binary_search(&call_index)
            .is_ok()
    {
        return true;
    }
    gidney_structural_dead_enabled() && GIDNEY_THREAD_BOUNDARY_DEAD_CALLS.contains(&call_index)
}

fn threaded_add_call_has_structurally_dead_sum(
    call_index: usize,
    bit: usize,
    total: usize,
    width: usize,
    vents: usize,
) -> bool {
    if gidney_skip_exact_sum_remainder_enabled()
        && GIDNEY_THREAD_SUM_REMAINDER_KEYS
            .binary_search(&gidney_key(call_index, bit))
            .is_ok()
    {
        return true;
    }
    if gidney_skip_fullvent_top2_enabled() && vents >= width && bit + 2 >= width {
        return true;
    }
    if gidney_skip_top2_thread_enabled() && bit + 2 >= total {
        return true;
    }
    gidney_structural_dead_enabled()
        && GIDNEY_THREAD_SUM_DEAD_RANGES
            .iter()
            .any(|&(call, lo, hi)| call == call_index && (lo..=hi).contains(&bit))
}

pub fn with_dirty_vent_pool<R>(dirty: &[QubitId], body: impl FnOnce() -> R) -> R {
    let count = std::env::var("TLM_DIRTY_VENTS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
        .min(dirty.len());
    let prior = DIRTY_VENT_POOL
        .with(|pool| std::mem::replace(&mut *pool.borrow_mut(), dirty[..count].to_vec()));
    let result = body();
    DIRTY_VENT_POOL.with(|pool| {
        *pool.borrow_mut() = prior;
    });
    result
}

fn dirty_vent_pool() -> Vec<QubitId> {
    DIRTY_VENT_POOL.with(|pool| pool.borrow().clone())
}

fn trace_schedule_fit(
    trace_env: &str,
    family: &str,
    mode: &str,
    fit: super::ScheduleFit,
    effective: usize,
    width: usize,
    entry_active: u32,
    timeline_start: usize,
    ops_start: usize,
    circ: &B,
) {
    if std::env::var_os(trace_env).is_none() {
        return;
    }
    let local_peak = circ.active_timeline[timeline_start..]
        .iter()
        .map(|(_, active)| *active)
        .max()
        .unwrap_or(entry_active);
    eprintln!(
        "TLM_{family} call={} phase={} mode={} width={} base={} selected={} effective={} entry_active={} local_peak={} ops_added={} ops={}",
        fit.call_index,
        circ.phase,
        mode,
        width,
        fit.base,
        fit.selected,
        effective,
        entry_active,
        local_peak,
        circ.current_ops_len().saturating_sub(ops_start),
        circ.current_ops_len(),
    );
}

pub fn controlled_hybrid_add_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId]) {
    controlled_hybrid_add_refs_impl(circ, ctrl, a, b, false);
}

fn controlled_hybrid_add_refs_skiplow(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
) {
    controlled_hybrid_add_refs_impl(circ, ctrl, a, b, true);
}

fn controlled_hybrid_add_refs_impl(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    skip_low_ctrl_sum: bool,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_hybrid_add: a, b must match width");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.ccx(*ctrl, *b[0], *a[0]);
        return;
    }
    let call_index = next_hybrid_add_call_index();

    let fit = super::next_hyb_v_fit();
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let ops_start = circ.current_ops_len();
    let vents = super::target_qubit_headroom(circ)
        .map_or(fit.selected, |headroom| fit.selected.min(headroom));

    for i in 1..n {
        circ.cx(*b[i], *a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(*b[i], *b[i + 1]);
    }

    #[derive(Clone, Copy)]
    enum VentLane {
        Clean(QubitId),
        Dirty(QubitId),
    }
    let dirty_pool = dirty_vent_pool();
    let mut vent_ancs: Vec<Option<VentLane>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            if let Some(&dirty) = dirty_pool.get(i) {
                debug_assert_ne!(dirty, *a[i]);
                debug_assert_ne!(dirty, *b[i]);
                debug_assert_ne!(dirty, *b[i + 1]);
                circ.cx(dirty, *b[i + 1]);
                let old_context = crate::point_add::set_op_trace_context(
                    0x1600_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(*a[i], *b[i], dirty);
                crate::point_add::restore_op_trace_context(old_context);
                circ.cx(dirty, *b[i + 1]);
                vent_ancs[i] = Some(VentLane::Dirty(dirty));
            } else {
                let anc = circ.alloc_qubit();
                let old_context = crate::point_add::set_op_trace_context(
                    0x1600_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(*a[i], *b[i], anc);
                crate::point_add::restore_op_trace_context(old_context);
                circ.cx(anc, *b[i + 1]);
                vent_ancs[i] = Some(VentLane::Clean(anc));
            }
        } else {
            let old_context = crate::point_add::set_op_trace_context(
                0x1600_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            circ.ccx(*a[i], *b[i], *b[i + 1]);
            crate::point_add::restore_op_trace_context(old_context);
        }
    }

    for i in (0..n - 1).rev() {
        let old_context = crate::point_add::set_op_trace_context(
            0x1700_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        circ.ccx(*ctrl, *b[i + 1], *a[i + 1]);
        crate::point_add::restore_op_trace_context(old_context);
        if let Some(lane) = vent_ancs[i].take() {
            match lane {
                VentLane::Clean(anc) => {
                    circ.cx(anc, *b[i + 1]);
                    let bit = circ.alloc_bit();
                    circ.hmr(anc, bit);
                    circ.zero_and_free(anc);
                    circ.cz_if_bit(*a[i], *b[i], bit);
                }
                VentLane::Dirty(dirty) => {
                    circ.cx(dirty, *b[i + 1]);
                    let old_context = crate::point_add::set_op_trace_context(
                        0x1800_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                    );
                    circ.ccx(*a[i], *b[i], dirty);
                    crate::point_add::restore_op_trace_context(old_context);
                    circ.cx(dirty, *b[i + 1]);
                }
            }
        } else {
            let old_context = crate::point_add::set_op_trace_context(
                0x1900_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            circ.ccx(*a[i], *b[i], *b[i + 1]);
            crate::point_add::restore_op_trace_context(old_context);
        }
    }

    for i in 1..n - 1 {
        circ.cx(*b[i], *b[i + 1]);
    }
    if !skip_low_ctrl_sum {
        let old_context = crate::point_add::set_op_trace_context(
            0x1a00_0000 | (((call_index as u32) & 0xffff) << 8),
        );
        circ.ccx(*ctrl, *b[0], *a[0]);
        crate::point_add::restore_op_trace_context(old_context);
    }
    for i in 1..n {
        circ.cx(*b[i], *a[i]);
    }
    trace_schedule_fit(
        "TRACE_TLM_HYB",
        "HYB",
        if skip_low_ctrl_sum {
            "skiplow"
        } else {
            "plain"
        },
        fit,
        vents,
        n,
        entry_active,
        timeline_start,
        ops_start,
        circ,
    );
}

fn controlled_clean_add_threaded(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
    vents: usize,
) {
    let call_index = next_threaded_add_call_index();
    let ops_start = circ.current_ops_len();
    let s = a.len();
    if s == 0 {
        if let (Some(ci), Some(co)) = (cin, cout) {
            circ.ccx(*ctrl, *ci, *co);
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
            circ.cx(*ci, *a[i]);
            circ.cx(*ci, *b[i]);
            if !threaded_add_call_has_structurally_dead_forward(call_index, i, n_inner, s, vents) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0500_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(*a[i], *b[i], *co);
                crate::point_add::restore_op_trace_context(old_context);
            }
            circ.cx(*ci, *co);
        } else {
            if !threaded_add_call_has_structurally_dead_forward(call_index, i, n_inner, s, vents) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0500_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(*a[i], *b[i], *co);
                crate::point_add::restore_op_trace_context(old_context);
            }
        }
    }

    if let Some(cout) = cout {
        let old_context = crate::point_add::set_op_trace_context(
            0x0600_0000
                | (((call_index as u32) & 0xffff) << 8)
                | ((s.saturating_sub(1)) as u32 & 0xff),
        );
        if !threaded_add_call_has_structurally_dead_boundary(call_index) {
            circ.ccx(*ctrl, *inner[s - 1].as_ref().unwrap(), *cout);
        }
        crate::point_add::restore_op_trace_context(old_context);
    }

    for i in (0..s).rev() {
        if !produces(i) {
            let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
            if let Some(ci) = ci {
                circ.cx(*ci, *b[i]);
            }
            if !threaded_add_call_has_structurally_dead_sum(call_index, i, s, s, vents) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x0700_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
                );
                circ.ccx(*ctrl, *b[i], *a[i]);
                crate::point_add::restore_op_trace_context(old_context);
            }
            if let Some(ci) = ci {
                circ.cx(*ci, *b[i]);
            }
            continue;
        }
        let co = inner[i].take().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, co);
        }
        if i < vents {
            let bit = circ.alloc_bit();
            circ.hmr(co, bit);
            circ.zero_and_free(co);
            circ.cz_if_bit(*a[i], *b[i], bit);
        } else {
            circ.ccx(*a[i], *b[i], co);
            circ.zero_and_free(co);
        }
        if let Some(ci) = ci {
            circ.cx(*ci, *a[i]);
        }
        if !threaded_add_call_has_structurally_dead_sum(call_index, i, s, s, vents) {
            let old_context = crate::point_add::set_op_trace_context(
                0x0700_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
            );
            circ.ccx(*ctrl, *b[i], *a[i]);
            crate::point_add::restore_op_trace_context(old_context);
        }
        if let Some(ci) = ci {
            circ.cx(*ci, *b[i]);
        }
    }
    if std::env::var_os("TRACE_TLM_GIDNEY_THREAD").is_some() {
        eprintln!(
            "TLM_GIDNEY_THREAD call={} phase={} width={} cin={} cout={} vents={} ops_start={} ops_end={}",
            call_index,
            circ.phase,
            s,
            usize::from(cin.is_some()),
            usize::from(cout.is_some()),
            vents,
            ops_start,
            circ.current_ops_len(),
        );
    }
}

fn deref(s: &[&QubitId]) -> Vec<QubitId> {
    s.iter().map(|q| **q).collect()
}

fn controlled_erase_carry_gated_impl(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: Option<&QubitId>,
    carry: QubitId,
) {
    let call_index = next_erase_gated_call_index();
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.push_condition(bit);
    let (av, bv) = (deref(a), deref(b));
    let ctrl = *ctrl;
    let cin = match cin {
        Some(cin) => {
            circ.loan_zero_qubit(carry);
            *cin
        }
        None => carry,
    };
    compare_geq_cin_middle(circ, &av, &bv, &cin, |c, ta, tb, c_prev| {
        c.z(ctrl);
        let old_context = crate::point_add::set_op_trace_context(
            0x0900_0000 | (((call_index as u32) & 0xffff) << 8),
        );
        if !gidney_erase_ccz_has_exact_dead_call(call_index) {
            c.ccz(ctrl, *ta, *tb);
        }
        crate::point_add::restore_op_trace_context(old_context);
        c.cz(ctrl, *c_prev);
    });
    circ.pop_condition();
    if cin == carry {
        circ.zero_and_free(carry);
    }
}

fn controlled_erase_carry_gated(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: &QubitId,
    carry: QubitId,
) {
    controlled_erase_carry_gated_impl(circ, ctrl, a, b, Some(cin), carry);
}

fn controlled_erase_carry_gated_zero_cin(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    carry: QubitId,
) {
    controlled_erase_carry_gated_impl(circ, ctrl, a, b, None, carry);
}

fn controlled_erase_carry_gated_capped(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: &QubitId,
    carry: QubitId,
    cap: usize,
) {
    let call_index = next_erase_gated_capped_call_index();
    let s = a.len();
    if s <= cap {
        controlled_erase_carry_gated(circ, ctrl, a, b, cin, carry);
        return;
    }
    let lo = s - cap;
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.push_condition(bit);
    let (av, bv) = (deref(&a[lo..]), deref(&b[lo..]));
    let ctrl = *ctrl;
    compare_geq_cin_middle(circ, &av, &bv, &carry, |c, ta, tb, c_prev| {
        c.z(ctrl);
        let old_context = crate::point_add::set_op_trace_context(
            0x0a00_0000 | (((call_index as u32) & 0xffff) << 8),
        );
        if !gidney_erase_capped_ccz_has_exact_dead_call(call_index) {
            c.ccz(ctrl, *ta, *tb);
        }
        crate::point_add::restore_op_trace_context(old_context);
        c.cz(ctrl, *c_prev);
    });
    circ.pop_condition();
    circ.zero_and_free(carry);
}

fn controlled_erase_carry_gated_capped_zero_cin(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    carry: QubitId,
    cap: usize,
) {
    if a.len() <= cap {
        controlled_erase_carry_gated_zero_cin(circ, ctrl, a, b, carry);
    } else {
        controlled_erase_carry_gated_capped(circ, ctrl, a, b, &carry, carry, cap);
    }
}

fn controlled_vented_chunk_add(
    circ: &mut B,
    ctrl: &QubitId,
    a_chunk: &[&QubitId],
    b_chunk: &[&QubitId],
    cin: &QubitId,
    cout: &QubitId,
) {
    let one = circ.alloc_qubit();
    circ.x(one);
    let zero = circ.alloc_qubit();
    let mut aext: Vec<&QubitId> = Vec::with_capacity(a_chunk.len() + 2);
    aext.push(&one);
    aext.extend_from_slice(a_chunk);
    aext.push(cout);
    let mut bext: Vec<&QubitId> = Vec::with_capacity(b_chunk.len() + 2);
    bext.push(cin);
    bext.extend_from_slice(b_chunk);
    bext.push(&zero);

    controlled_hybrid_add_refs_skiplow(circ, ctrl, &aext, &bext);
    circ.x(one);
    circ.zero_and_free(one);
    circ.zero_and_free(zero);
}

fn varchunk_schedule(n: usize, k: usize) -> Vec<usize> {
    const RESERVE: usize = 4;
    let mut sizes = Vec::new();
    let (mut covered, mut held) = (0usize, 0usize);
    while covered < n {
        let room = k.saturating_sub(held + RESERVE);
        if room == 0 {
            return Vec::new();
        }
        let s = room.min(n - covered);
        sizes.push(s);
        covered += s;
        held += 1;
    }
    sizes
}

fn varchunk_cost(n: usize, k: usize, cap: usize) -> usize {
    let sizes = varchunk_schedule(n, k);
    if sizes.is_empty() {
        return usize::MAX;
    }
    let erase: usize = sizes.iter().map(|&s| s.min(cap) / 2).sum();
    n + erase
}

pub(crate) struct AdaptiveLayout {
    pub(crate) c: usize,
    pub(crate) chunked_len: usize,
    pub(crate) plain_len: usize,
}
pub(crate) const ADAPTIVE_RES: usize = 5;

fn adaptive_chunk_size(n: usize) -> usize {
    std::env::var("TLM_ADAPTIVE_CHUNK")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v > 0)
        .unwrap_or_else(|| (n as f64).sqrt() as usize)
        .clamp(1, n)
}

pub(crate) fn adaptive_layout(n: usize, k: usize) -> AdaptiveLayout {
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    adaptive_layout_for_chunk(n, k, c)
}

fn adaptive_layout_for_chunk(n: usize, k: usize, c: usize) -> AdaptiveLayout {
    let mut plain = 0usize;
    while plain < n {
        let l = n - (plain + 1);
        let nch = l.div_ceil(c);
        if nch + (plain + 1) <= k {
            plain += 1;
        } else {
            break;
        }
    }
    AdaptiveLayout {
        c,
        chunked_len: n - plain,
        plain_len: plain,
    }
}

fn searched_cout_layout(n: usize, k: usize) -> Option<AdaptiveLayout> {
    if std::env::var_os("TLM_COUT_LAYOUT_SEARCH").is_none() {
        return None;
    }
    let mut margin = std::env::var("TLM_COUT_LAYOUT_MARGIN")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1);
    if margin == 0
        && std::env::var("TLM_COUT_LAYOUT_FORCE_M1_KS")
            .ok()
            .map(|s| {
                s.split(',')
                    .filter_map(|part| part.trim().parse::<usize>().ok())
                    .any(|force_k| force_k == k)
            })
            .unwrap_or(false)
    {
        margin = 1;
    }
    let mut best: Option<(usize, AdaptiveLayout)> = None;
    for c in 1..=n {
        for plain_len in 0..=n {
            let chunked_len = n - plain_len;
            let nchunks = chunked_len.div_ceil(c);
            if nchunks + plain_len + margin > k {
                continue;
            }
            if nchunks + c.min(chunked_len.max(1)) + margin > k {
                continue;
            }
            let cost = 2 * n + chunked_len + nchunks + 1;
            let layout = AdaptiveLayout {
                c,
                chunked_len,
                plain_len,
            };
            match best {
                Some((best_cost, _)) if best_cost <= cost => {}
                _ => best = Some((cost, layout)),
            }
        }
    }
    best.map(|(_, layout)| layout)
}

fn emit_cout_layout(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cout: &QubitId,
    layout: AdaptiveLayout,
) {
    let n = a.len();
    let l = layout.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + layout.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: Option<&QubitId> = if j == 0 { None } else { Some(&carries[j - 1]) };
        controlled_clean_add_threaded(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, Some(&cy), hi - lo);
        carries.push(cy);
    }
    controlled_clean_add_threaded(
        circ,
        ctrl,
        &a[l..n],
        &b[l..n],
        carries.last(),
        Some(cout),
        layout.plain_len,
    );
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_zero_cin(circ, ctrl, &a[lo..hi], &b[lo..hi], carry);
        } else {
            controlled_erase_carry_gated(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                &carries[j - 1],
                carry,
            );
        }
    }
}

fn searched_gcd_adaptive_layout(n: usize, k: usize) -> Option<AdaptiveLayout> {
    if std::env::var_os("TLM_GCD_ADAPTIVE_LAYOUT_SEARCH").is_none() {
        return None;
    }
    let margin = std::env::var("TLM_GCD_ADAPTIVE_LAYOUT_MARGIN")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1);
    let mut best: Option<(usize, AdaptiveLayout)> = None;
    for c in 1..=n {
        for plain_len in 0..=n {
            let chunked_len = n - plain_len;
            let nchunks = chunked_len.div_ceil(c);
            if nchunks + plain_len + margin > k {
                continue;
            }
            if nchunks + c.min(chunked_len.max(1)) + margin > k {
                continue;
            }
            let cost = 2 * n + chunked_len + nchunks - 1;
            let layout = AdaptiveLayout {
                c,
                chunked_len,
                plain_len,
            };
            match best {
                Some((best_cost, _)) if best_cost <= cost => {}
                _ => best = Some((cost, layout)),
            }
        }
    }
    best.map(|(_, layout)| layout)
}

fn emit_adaptive_layout_no_cout(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    layout: AdaptiveLayout,
) {
    let n = a.len();
    let l = layout.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + layout.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    let cin0 = circ.alloc_qubit();
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cout = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_clean_add_threaded(
            circ,
            ctrl,
            &a[lo..hi],
            &b[lo..hi],
            Some(cin),
            Some(&cout),
            hi - lo,
        );
        carries.push(cout);
    }
    if layout.plain_len > 0 {
        let cin: &QubitId = carries.last().unwrap_or(&cin0);
        controlled_clean_add_threaded(
            circ,
            ctrl,
            &a[l..n],
            &b[l..n],
            Some(cin),
            None,
            layout.plain_len,
        );
    }
    circ.zero_and_free(cin0);
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_zero_cin(circ, ctrl, &a[lo..hi], &b[lo..hi], carry);
        } else {
            controlled_erase_carry_gated(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                &carries[j - 1],
                carry,
            );
        }
    }
}

fn adaptive_add_cost_tof(n: usize, k: usize, controlled: bool) -> u64 {
    if n == 0 {
        return 0;
    }
    let base = if controlled { 3 * n } else { 2 * n };
    let s2 = 2 * (n as f64).sqrt() as usize;
    let saved = if k >= n {
        n
    } else if k < s2 {
        (k * k) / 8
    } else {
        n / 2 + (k - s2) / 2
    };
    (base.saturating_sub(saved)) as u64
}

fn controlled_chunked_then_cuccaro(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cout: Option<&QubitId>,
    k: usize,
) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let cin0 = circ.alloc_qubit();
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (mut lo, mut i) = (0usize, 0usize);
    while lo < n && k > i + 2 {
        let cc = (k - 2 - i).min(n - lo);
        bounds.push((lo, lo + cc));
        lo += cc;
        i += 1;
    }
    let chunked_len = lo;
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(clo, chi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_clean_add_threaded(
            circ,
            ctrl,
            &a[clo..chi],
            &b[clo..chi],
            Some(cin),
            Some(&cy),
            chi - clo,
        );
        carries.push(cy);
    }
    if chunked_len < n {
        let cin: &QubitId = carries.last().unwrap_or(&cin0);

        let at = deref(&a[chunked_len..n]);
        let bt = deref(&b[chunked_len..n]);
        super::arith::cuccaro_carry(circ, Some(ctrl), &bt, &at, Some(cin), cout);
    } else if let Some(co) = cout {
        circ.cx(*carries.last().unwrap_or(&cin0), *co);
    }
    circ.zero_and_free(cin0);
    for j in (0..bounds.len()).rev() {
        let (clo, chi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_zero_cin(circ, ctrl, &a[clo..chi], &b[clo..chi], carry);
        } else {
            controlled_erase_carry_gated(
                circ,
                ctrl,
                &a[clo..chi],
                &b[clo..chi],
                &carries[j - 1],
                carry,
            );
        }
    }
}

fn controlled_hybrid_add_adaptive_refs(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    k: usize,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled adaptive add: a,b width mismatch");
    if n == 0 {
        return;
    }
    if let Some(layout) = searched_gcd_adaptive_layout(n, k) {
        emit_adaptive_layout_no_cout(circ, ctrl, a, b, layout);
        return;
    }
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    if n <= 4 || k.saturating_add(2 * c) >= n {
        controlled_hybrid_add_refs(circ, ctrl, a, b);
        return;
    }
    let tight = k < n.div_ceil(c) + c + ADAPTIVE_RES;
    let cov = (k.saturating_sub(2).saturating_mul(k.saturating_sub(1)) / 2).min(n);
    if tight && cov < n {
        if cov > 2 * k {
            controlled_chunked_then_cuccaro(circ, ctrl, a, b, None, k);
        } else {
            controlled_hybrid_add_refs(circ, ctrl, a, b);
        }
        return;
    }
    let cin0 = circ.alloc_qubit();
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (l, plain_len) = if tight {
        let (mut lo, mut i) = (0usize, 0usize);
        while lo < n && k > i + 2 {
            let cc = (k - 2 - i).min(n - lo);
            bounds.push((lo, lo + cc));
            lo += cc;
            i += 1;
        }
        (n, 0)
    } else {
        let lay = adaptive_layout(n, k);
        let mut lo = 0;
        while lo < lay.chunked_len {
            let hi = (lo + lay.c).min(lay.chunked_len);
            bounds.push((lo, hi));
            lo = hi;
        }
        (lay.chunked_len, lay.plain_len)
    };
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cout = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_clean_add_threaded(
            circ,
            ctrl,
            &a[lo..hi],
            &b[lo..hi],
            Some(cin),
            Some(&cout),
            hi - lo,
        );
        carries.push(cout);
    }
    if plain_len > 0 {
        let cin: &QubitId = carries.last().unwrap_or(&cin0);
        controlled_clean_add_threaded(circ, ctrl, &a[l..n], &b[l..n], Some(cin), None, plain_len);
    }
    circ.zero_and_free(cin0);
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_zero_cin(circ, ctrl, &a[lo..hi], &b[lo..hi], carry);
        } else {
            controlled_erase_carry_gated(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                &carries[j - 1],
                carry,
            );
        }
    }
}

fn controlled_hybrid_add_varchunk_gated_refs(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    k: usize,
    cap: usize,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "varchunk add: a,b width mismatch");
    if n == 0 {
        return;
    }
    let sizes = varchunk_schedule(n, k);
    assert!(!sizes.is_empty(), "varchunk infeasible at k={k} for n={n}");
    let direct = std::env::var("TLM_DIRECT_VARCHUNK").ok().as_deref() == Some("1");
    let cin0 = (!direct).then(|| circ.alloc_qubit());
    let mut carries: Vec<QubitId> = Vec::with_capacity(sizes.len());
    let mut bounds: Vec<(usize, usize)> = Vec::with_capacity(sizes.len());
    let mut lo = 0usize;
    for (j, &s) in sizes.iter().enumerate() {
        let hi = lo + s;
        let cout = circ.alloc_qubit();
        if direct {
            let fit = super::next_hyb_v_fit();
            let timeline_start = circ.active_timeline.len();
            let entry_active = circ.active_qubits;
            let ops_start = circ.current_ops_len();
            let cin = (j != 0).then(|| &carries[j - 1]);
            controlled_clean_add_threaded(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                cin,
                Some(&cout),
                hi - lo,
            );
            trace_schedule_fit(
                "TRACE_TLM_HYB",
                "HYB",
                "direct-varchunk",
                fit,
                hi - lo,
                hi - lo,
                entry_active,
                timeline_start,
                ops_start,
                circ,
            );
        } else {
            let cin: &QubitId = if j == 0 {
                cin0.as_ref().expect("legacy cin0")
            } else {
                &carries[j - 1]
            };
            controlled_vented_chunk_add(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, &cout);
        }
        carries.push(cout);
        bounds.push((lo, hi));
        lo = hi;
    }
    if let Some(cin0) = cin0 {
        circ.zero_and_free(cin0);
    }
    for j in (0..sizes.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_capped_zero_cin(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                carry,
                cap,
            );
        } else {
            controlled_erase_carry_gated_capped(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                &carries[j - 1],
                carry,
                cap,
            );
        }
    }
}

fn controlled_hybrid_add_knob_capped_refs(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    k: usize,
    cap: usize,
) {
    let n = a.len();
    if cap < n
        && !varchunk_schedule(n, k).is_empty()
        && (varchunk_cost(n, k, cap) as u64 + n as u64) < adaptive_add_cost_tof(n, k, true)
    {
        controlled_hybrid_add_varchunk_gated_refs(circ, ctrl, a, b, k, cap);
    } else {
        controlled_hybrid_add_adaptive_refs(circ, ctrl, a, b, k);
    }
}

pub fn controlled_hybrid_add_capped_branch(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    k: usize,
    cap: usize,
    branch: u8,
) {
    let n = a.len();
    let k = super::target_qubit_headroom(circ).map_or(k, |headroom| k.min(headroom));
    if std::env::var("TLM_GCD_RESELECT_LAYOUT").ok().as_deref() == Some("1") {
        controlled_hybrid_add_knob_capped_refs(circ, ctrl, a, b, k, cap);
        return;
    }
    if branch == 1 && n > 0 && !varchunk_schedule(n, k).is_empty() {
        controlled_hybrid_add_varchunk_gated_refs(circ, ctrl, a, b, k, cap);
    } else if branch == 0 {
        controlled_hybrid_add_refs(circ, ctrl, a, b);
    } else if branch == 255 {
        controlled_hybrid_add_knob_capped_refs(circ, ctrl, a, b, k, cap);
    } else {
        controlled_hybrid_add_adaptive_refs(circ, ctrl, a, b, k);
    }
}

pub fn controlled_hybrid_add_cout_refs(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cout: &QubitId,
    k: usize,
) {
    let fit = super::take_cout_fit(k);
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let ops_start = circ.current_ops_len();
    let forced = std::env::var("TLM_FORCE_COUT_K")
        .ok()
        .and_then(|v| v.parse::<usize>().ok());
    let effective = match forced {
        Some(v) => v,
        None => super::target_qubit_headroom(circ)
            .map_or(fit.selected, |headroom| fit.selected.min(headroom)),
    };
    controlled_hybrid_add_cout_refs_impl(circ, ctrl, a, b, cout, effective);
    trace_schedule_fit(
        "TRACE_TLM_COUT",
        "COUT",
        "dispatch",
        fit,
        effective,
        a.len(),
        entry_active,
        timeline_start,
        ops_start,
        circ,
    );
}

fn controlled_hybrid_add_cout_refs_impl(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cout: &QubitId,
    k: usize,
) {
    let n = a.len();
    if std::env::var("TRACE_TLM_COUT_N").is_ok() {
        eprintln!("TLM_COUT_N n={} k={}", n, k);
    }
    assert_eq!(b.len(), n, "controlled cout add: a,b width mismatch");
    assert!(n >= 1, "controlled cout add: empty operands");
    if let Some(layout) = searched_cout_layout(n, k) {
        if std::env::var("TRACE_TLM_COUT_PATH").is_ok() {
            eprintln!("TLM_COUT_PATH n={} k={} path=searched_env c={} plain={}", n, k, layout.c, layout.plain_len);
        }
        emit_cout_layout(circ, ctrl, a, b, cout, layout);
        return;
    }
    let c = adaptive_chunk_size(n);
    let lay = adaptive_layout_for_chunk(n, k, c);
    let tight = k < n.div_ceil(c) + c + ADAPTIVE_RES;
    let cov = (k.saturating_sub(2).saturating_mul(k.saturating_sub(1)) / 2).min(n);
    if n > 4 && k.saturating_add(2 * c) < n && tight && cov > 2 * k {
        if std::env::var("TRACE_TLM_COUT_PATH").is_ok() {
            eprintln!("TLM_COUT_PATH n={} k={} c={} tight={} path=chunked_cuccaro", n, k, c, tight);
        }
        controlled_chunked_then_cuccaro(circ, ctrl, a, b, Some(cout), k);
        return;
    }
    if n <= 4
        || k < n.div_ceil(c) + c + ADAPTIVE_RES
        || k.saturating_add(2 * c) >= n
        || lay.plain_len == 0
    {
        if std::env::var("TRACE_TLM_COUT_PATH").is_ok() {
            eprintln!("TLM_COUT_PATH n={} k={} path=pure_ks", n, k);
        }
        let zpad = circ.alloc_qubit();
        let mut aref: Vec<&QubitId> = a.to_vec();
        aref.push(cout);
        let mut bref: Vec<&QubitId> = b.to_vec();
        bref.push(&zpad);
        controlled_hybrid_add_refs(circ, ctrl, &aref, &bref);
        circ.zero_and_free(zpad);
        return;
    }
    let l = lay.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + lay.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: Option<&QubitId> = if j == 0 { None } else { Some(&carries[j - 1]) };
        controlled_clean_add_threaded(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, Some(&cy), hi - lo);
        carries.push(cy);
    }
    controlled_clean_add_threaded(
        circ,
        ctrl,
        &a[l..n],
        &b[l..n],
        carries.last(),
        Some(cout),
        lay.plain_len,
    );
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            controlled_erase_carry_gated_zero_cin(circ, ctrl, &a[lo..hi], &b[lo..hi], carry);
        } else {
            controlled_erase_carry_gated(
                circ,
                ctrl,
                &a[lo..hi],
                &b[lo..hi],
                &carries[j - 1],
                carry,
            );
        }
    }
}
