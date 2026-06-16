//! Independent kmx parser + bitsliced (64-shot) simulator, using
//! trailmix's OWN gate semantics. Mirrors the zenodo fuzz harness:
//!   kmx_resim <circuit.kmx>   with "a b -> c d" cases on stdin.
//!
//! Purpose: if trailmix RE-SIMULATES its own emitted kmx and gets a
//! DIFFERENT answer than its construction-time sim (the cases' expected
//! column is the construction output), the serialization is lossy — some
//! op's effect is applied during construction but not captured in the
//! emitted op stream. If it matches, the text round-trips faithfully and
//! any zenodo divergence is a simulator-semantics gap instead.

use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
use std::io::{stdin, BufRead};

#[derive(Clone, Copy)]
enum Op {
    X(u32),
    Z(u32),
    Cx(u32, u32),
    Cz(u32, u32),
    Ccx(u32, u32, u32),
    Ccz(u32, u32, u32),
    Swap(u32, u32),
    Hmr(u32, u32),
    R(u32),
    Neg,
    Push(u32),
    Pop,
    BitInvert(u32),
    BitStore0(u32),
    BitStore1(u32),
}

struct Parsed {
    ops: Vec<(Op, Option<u32>)>, // (op, inline-condition-bit)
    regs: Vec<Vec<(bool, u32)>>, // per register: (is_qubit, id)
    max_q: u32,
    max_b: u32,
}

fn pq(s: &str) -> u32 {
    s.trim_start_matches('q').parse().unwrap()
}
fn pb(s: &str) -> u32 {
    s.trim_start_matches('b').parse().unwrap()
}

fn parse(kmx: &str) -> Parsed {
    let mut ops = Vec::new();
    let mut regs: Vec<Vec<(bool, u32)>> = Vec::new();
    let mut max_q = 0u32;
    let mut max_b = 0u32;
    for line in kmx.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        // inline "... if b<bit>" suffix (but NOT for PUSH_CONDITION,
        // whose native form is literally "PUSH_CONDITION if b<bit>").
        let (toks, cond): (&[&str], Option<u32>) =
            if t[0] != "PUSH_CONDITION" && t.len() >= 2 && t[t.len() - 2] == "if" {
                (&t[..t.len() - 2], Some(pb(t[t.len() - 1])))
            } else {
                (&t[..], None)
            };
        let track_q = |id: u32, mq: &mut u32| {
            if id > *mq {
                *mq = id;
            }
        };
        match toks[0] {
            "REGISTER" => regs.push(Vec::new()),
            "APPEND_TO_REGISTER" => {
                let item = toks[1];
                let reg: usize = toks[2].trim_start_matches('r').parse().unwrap();
                while regs.len() <= reg {
                    regs.push(Vec::new());
                }
                if let Some(ids) = item.strip_prefix('q') {
                    let id: u32 = ids.parse().unwrap();
                    track_q(id, &mut max_q);
                    regs[reg].push((true, id));
                } else {
                    let id: u32 = item.trim_start_matches('b').parse().unwrap();
                    if id > max_b {
                        max_b = id;
                    }
                    regs[reg].push((false, id));
                }
            }
            "X" => {
                let q = pq(toks[1]);
                track_q(q, &mut max_q);
                ops.push((Op::X(q), cond));
            }
            "Z" => {
                let q = pq(toks[1]);
                track_q(q, &mut max_q);
                ops.push((Op::Z(q), cond));
            }
            "CX" => {
                let (a, b) = (pq(toks[1]), pq(toks[2]));
                track_q(a, &mut max_q);
                track_q(b, &mut max_q);
                ops.push((Op::Cx(a, b), cond));
            }
            "CZ" => {
                let (a, b) = (pq(toks[1]), pq(toks[2]));
                track_q(a, &mut max_q);
                track_q(b, &mut max_q);
                ops.push((Op::Cz(a, b), cond));
            }
            "CCX" => {
                let (a, b, c) = (pq(toks[1]), pq(toks[2]), pq(toks[3]));
                track_q(a, &mut max_q);
                track_q(b, &mut max_q);
                track_q(c, &mut max_q);
                ops.push((Op::Ccx(a, b, c), cond));
            }
            "CCZ" => {
                let (a, b, c) = (pq(toks[1]), pq(toks[2]), pq(toks[3]));
                track_q(a, &mut max_q);
                track_q(b, &mut max_q);
                track_q(c, &mut max_q);
                ops.push((Op::Ccz(a, b, c), cond));
            }
            "SWAP" => {
                let (a, b) = (pq(toks[1]), pq(toks[2]));
                track_q(a, &mut max_q);
                track_q(b, &mut max_q);
                ops.push((Op::Swap(a, b), cond));
            }
            "HMR" => {
                let q = pq(toks[1]);
                let b = pb(toks[2]);
                track_q(q, &mut max_q);
                if b > max_b {
                    max_b = b;
                }
                ops.push((Op::Hmr(q, b), cond));
            }
            "R" => {
                let q = pq(toks[1]);
                track_q(q, &mut max_q);
                ops.push((Op::R(q), cond));
            }
            "NEG" => ops.push((Op::Neg, cond)),
            "PUSH_CONDITION" => {
                let b = pb(toks[2]);
                if b > max_b {
                    max_b = b;
                }
                ops.push((Op::Push(b), cond));
            }
            "POP_CONDITION" => ops.push((Op::Pop, cond)),
            "BIT_INVERT" => {
                let b = pb(toks[1]);
                if b > max_b {
                    max_b = b;
                }
                ops.push((Op::BitInvert(b), cond));
            }
            "BIT_STORE0" => {
                let b = pb(toks[1]);
                if b > max_b {
                    max_b = b;
                }
                ops.push((Op::BitStore0(b), cond));
            }
            "BIT_STORE1" => {
                let b = pb(toks[1]);
                if b > max_b {
                    max_b = b;
                }
                ops.push((Op::BitStore1(b), cond));
            }
            other => panic!("unknown kmx op: {other} (line: {line})"),
        }
    }
    Parsed {
        ops,
        regs,
        max_q,
        max_b,
    }
}

struct Sim {
    q: Vec<u64>,
    b: Vec<u64>,
    phase: u64,
    rng: u64,
}

impl Sim {
    fn next_rng(&mut self) -> u64 {
        // splitmix64
        self.rng = self.rng.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.rng;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    fn run(&mut self, ops: &[(Op, Option<u32>)]) {
        let mut stack: Vec<u64> = Vec::new();
        let mut base: u64 = u64::MAX;
        for (op, inline) in ops {
            let mut cond = base;
            if let Some(cb) = inline {
                cond &= self.b[*cb as usize];
            }
            match *op {
                Op::X(q) => self.q[q as usize] ^= cond,
                Op::Z(q) => self.phase ^= cond & self.q[q as usize],
                Op::Cx(c, t) => self.q[t as usize] ^= cond & self.q[c as usize],
                Op::Cz(a, b) => self.phase ^= cond & self.q[a as usize] & self.q[b as usize],
                Op::Ccx(a, b, t) => {
                    self.q[t as usize] ^= cond & self.q[a as usize] & self.q[b as usize]
                }
                Op::Ccz(a, b, c) => {
                    self.phase ^=
                        cond & self.q[a as usize] & self.q[b as usize] & self.q[c as usize]
                }
                Op::Swap(a, b) => {
                    let (ai, bi) = (a as usize, b as usize);
                    let mut qa = self.q[ai];
                    let mut qb = self.q[bi];
                    qa ^= qb;
                    qb ^= cond & qa;
                    qa ^= qb;
                    self.q[ai] = qa;
                    self.q[bi] = qb;
                }
                Op::Hmr(q, bbit) => {
                    let r = self.next_rng() & cond;
                    self.b[bbit as usize] = r;
                    self.phase ^= self.q[q as usize] & r;
                    self.q[q as usize] = 0;
                }
                Op::R(q) => {
                    let r = self.next_rng();
                    self.phase ^= self.q[q as usize] & r & cond;
                    self.q[q as usize] = 0;
                }
                Op::Neg => self.phase ^= cond,
                Op::Push(b) => {
                    stack.push(base);
                    base &= self.b[b as usize];
                }
                Op::Pop => {
                    if let Some(v) = stack.pop() {
                        base = v;
                    }
                }
                Op::BitInvert(b) => self.b[b as usize] ^= cond,
                Op::BitStore0(b) => self.b[b as usize] &= !cond,
                Op::BitStore1(b) => self.b[b as usize] |= cond,
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let kmx = std::fs::read_to_string(&args[1]).expect("read kmx");
    let parsed = parse(&kmx);
    let nq = (parsed.max_q + 1) as usize;
    let nb = (parsed.max_b + 1) as usize;
    eprintln!(
        "[resim] {} ops, {} regs, {} qubits, {} bits",
        parsed.ops.len(),
        parsed.regs.len(),
        nq,
        nb
    );

    let mut block_in: Vec<Vec<BigUint>> = Vec::new();
    let mut block_exp: Vec<Vec<BigUint>> = Vec::new();
    let mut shot = 0usize;
    let mut shots_total = 0usize;
    let mut failures = 0usize;

    let set_reg = |sim: &mut Sim, reg: &[(bool, u32)], val: &BigUint, s: usize| {
        for (i, (is_q, id)) in reg.iter().enumerate() {
            let bit = val.bit(i as u64) as u64;
            if *is_q {
                sim.q[*id as usize] = (sim.q[*id as usize] & !(1u64 << s)) | (bit << s);
            } else {
                sim.b[*id as usize] = (sim.b[*id as usize] & !(1u64 << s)) | (bit << s);
            }
        }
    };
    let get_reg = |sim: &Sim, reg: &[(bool, u32)], s: usize| -> BigUint {
        let mut v = BigUint::from(0u32);
        for (i, (is_q, id)) in reg.iter().enumerate() {
            let bit = if *is_q {
                (sim.q[*id as usize] >> s) & 1
            } else {
                (sim.b[*id as usize] >> s) & 1
            };
            if bit == 1 {
                v |= BigUint::from(1u32) << i;
            }
        }
        v
    };

    let run_block = |block_in: &Vec<Vec<BigUint>>,
                     block_exp: &Vec<Vec<BigUint>>,
                     failures: &mut usize,
                     shots_total: usize| {
        let seed: u64 = std::env::var("SEED")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0x1234_5678);
        let mut sim = Sim {
            q: vec![0u64; nq],
            b: vec![0u64; nb],
            phase: 0,
            rng: seed,
        };
        for (s, ins) in block_in.iter().enumerate() {
            for (k, reg) in parsed.regs.iter().enumerate() {
                set_reg(&mut sim, reg, &ins[k], s);
            }
        }
        sim.run(&parsed.ops);
        for s in 0..block_in.len() {
            let mut bad = (sim.phase >> s) & 1 != 0;
            let mut got = Vec::new();
            for (k, reg) in parsed.regs.iter().enumerate() {
                let g = get_reg(&sim, reg, s);
                if g != block_exp[s][k] {
                    bad = true;
                }
                got.push(g);
            }
            if bad {
                eprintln!(
                    "RESIM FAIL shot {}: got {:?} exp {:?} phase_bit {}",
                    shots_total + s,
                    got,
                    block_exp[s],
                    (sim.phase >> s) & 1
                );
                *failures += 1;
                if *failures >= 3 {
                    return;
                }
            }
        }
    };

    for line in stdin().lock().lines() {
        let line = line.expect("stdin");
        let (inp, out) = line.split_once(" -> ").unwrap();
        let ins: Vec<BigUint> = inp.split_whitespace().map(|s| s.parse().unwrap()).collect();
        let outs: Vec<BigUint> = out.split_whitespace().map(|s| s.parse().unwrap()).collect();
        block_in.push(ins);
        block_exp.push(outs);
        shot += 1;
        if shot == 64 {
            run_block(&block_in, &block_exp, &mut failures, shots_total);
            shots_total += 64;
            shot = 0;
            block_in.clear();
            block_exp.clear();
            if failures > 0 {
                break;
            }
        }
    }
    if shot > 0 && failures == 0 {
        run_block(&block_in, &block_exp, &mut failures, shots_total);
        shots_total += shot;
    }
    if failures == 0 {
        println!(
            "RESIM PASS ({} shots) — kmx round-trips faithfully under trailmix semantics",
            shots_total
        );
    } else {
        println!(
            "RESIM FAILED ({} failures) — kmx is lossy under trailmix's own semantics",
            failures
        );
    }
}
