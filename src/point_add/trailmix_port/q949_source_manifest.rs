// Every packaged source that can change Q949 allocation, operation
// serialization, Fiat-Shamir sampling, or embedded support data belongs here.
// Auxiliary proof scripts and WMI job files are deliberately excluded because
// the challenge packages only src/point_add.
pub(super) const SOURCES: &[(&str, &[u8])] = &[
    ("Cargo.toml", include_bytes!("../../../Cargo.toml")),
    ("Cargo.lock", include_bytes!("../../../Cargo.lock")),
    ("rust-toolchain", include_bytes!("../../../rust-toolchain")),
    ("src/lib.rs", include_bytes!("../../lib.rs")),
    ("src/circuit.rs", include_bytes!("../../circuit.rs")),
    ("src/sim.rs", include_bytes!("../../sim.rs")),
    (
        "src/weierstrass_elliptic_curve.rs",
        include_bytes!("../../weierstrass_elliptic_curve.rs"),
    ),
    ("src/point_add/mod.rs", include_bytes!("../mod.rs")),
    ("src/point_add/emit.rs", include_bytes!("../emit.rs")),
    ("src/point_add/venting.rs", include_bytes!("../venting.rs")),
    ("src/point_add/trailmix_port/mod.rs", include_bytes!("mod.rs")),
    (
        "src/point_add/trailmix_port/q949_source_manifest.rs",
        include_bytes!("q949_source_manifest.rs"),
    ),
    (
        "src/point_add/trailmix_port/circuit.rs",
        include_bytes!("circuit.rs"),
    ),
    (
        "src/point_add/trailmix_port/mod_arith.rs",
        include_bytes!("mod_arith.rs"),
    ),
    (
        "src/point_add/trailmix_port/rfold_mbu.rs",
        include_bytes!("rfold_mbu.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/compare.rs",
        include_bytes!("arith/compare.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/const_add.rs",
        include_bytes!("arith/const_add.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/cuccaro.rs",
        include_bytes!("arith/cuccaro.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/gidney_const_adder.rs",
        include_bytes!("arith/gidney_const_adder.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/khattar_gidney.rs",
        include_bytes!("arith/khattar_gidney.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/mcx.rs",
        include_bytes!("arith/mcx.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/qshift_sub.rs",
        include_bytes!("arith/qshift_sub.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/ripple_add.rs",
        include_bytes!("arith/ripple_add.rs"),
    ),
    (
        "src/point_add/trailmix_port/arith/shift.rs",
        include_bytes!("arith/shift.rs"),
    ),
    (
        "src/point_add/trailmix_port/ec/point_add.rs",
        include_bytes!("ec/point_add.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/shrunken_pz_primitives.rs",
        include_bytes!("inversion/shrunken_pz_primitives.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/shrunken_pz_schedule.rs",
        include_bytes!("inversion/shrunken_pz_schedule.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/shrunken_pz_state_machine.rs",
        include_bytes!("inversion/shrunken_pz_state_machine.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/q949_robust_envelope.rs",
        include_bytes!("inversion/q949_robust_envelope.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/q949_robust_envelope_data.rs",
        include_bytes!("inversion/q949_robust_envelope_data.rs"),
    ),
    (
        "src/point_add/trailmix_port/inversion/q949_robust_projection_metadata.json",
        include_bytes!("inversion/q949_robust_projection_metadata.json"),
    ),
];
