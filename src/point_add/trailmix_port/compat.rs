pub mod num_bigint {
    use core::fmt;
    use core::ops::{
        Add, AddAssign, BitAnd, BitOr, BitOrAssign, Div, Mul, Rem, Shl, Shr, Sub, SubAssign,
    };
    use ruint::aliases::U512;

    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct BigUint(pub U512);

    impl BigUint {
        #[must_use]
        pub fn from_bytes_le(bytes: &[u8]) -> Self {
            Self(U512::from_le_slice(bytes))
        }

        #[must_use]
        pub fn parse_bytes(bytes: &[u8], radix: u32) -> Option<Self> {
            let s = core::str::from_utf8(bytes).ok()?;
            U512::from_str_radix(s, radix.into()).ok().map(Self)
        }

        #[must_use]
        pub fn to_bytes_le(&self) -> Vec<u8> {
            self.0.to_le_bytes_trimmed_vec()
        }

        #[must_use]
        pub fn modpow(&self, exp: &Self, modulus: &Self) -> Self {
            if modulus.0 == U512::ZERO {
                return Self(U512::ZERO);
            }
            let mut base = self.0 % modulus.0;
            let mut acc = U512::from(1u64) % modulus.0;
            for i in 0..512 {
                if exp.0.bit(i) {
                    acc = (acc * base) % modulus.0;
                }
                base = base.wrapping_mul(base) % modulus.0;
            }
            Self(acc)
        }

        #[must_use]
        pub fn zero() -> Self {
            Self(U512::ZERO)
        }

        #[must_use]
        pub fn is_zero(&self) -> bool {
            self.0 == U512::ZERO
        }
    }

    macro_rules! from_int {
        ($($t:ty),* $(,)?) => {
            $(impl From<$t> for BigUint {
                fn from(value: $t) -> Self {
                    Self(U512::from(value as u64))
                }
            })*
        };
    }
    from_int!(u8, u16, u32, u64, usize, i32);

    macro_rules! binop {
        ($trait:ident, $method:ident, $op:tt) => {
            impl $trait for BigUint {
                type Output = BigUint;
                fn $method(self, rhs: BigUint) -> BigUint {
                    BigUint(self.0 $op rhs.0)
                }
            }
            impl $trait<&BigUint> for BigUint {
                type Output = BigUint;
                fn $method(self, rhs: &BigUint) -> BigUint {
                    BigUint(self.0 $op rhs.0)
                }
            }
            impl $trait<BigUint> for &BigUint {
                type Output = BigUint;
                fn $method(self, rhs: BigUint) -> BigUint {
                    BigUint(self.0 $op rhs.0)
                }
            }
            impl $trait<&BigUint> for &BigUint {
                type Output = BigUint;
                fn $method(self, rhs: &BigUint) -> BigUint {
                    BigUint(self.0 $op rhs.0)
                }
            }
        };
    }

    binop!(Add, add, +);
    binop!(Sub, sub, -);
    binop!(Mul, mul, *);
    binop!(Div, div, /);
    binop!(Rem, rem, %);
    binop!(BitAnd, bitand, &);
    binop!(BitOr, bitor, |);

    macro_rules! binop_int {
        ($trait:ident, $method:ident, $op:tt, $($t:ty),* $(,)?) => {
            $(
                impl $trait<$t> for BigUint {
                    type Output = BigUint;
                    fn $method(self, rhs: $t) -> BigUint {
                        BigUint(self.0 $op U512::from(rhs as u64))
                    }
                }
                impl $trait<$t> for &BigUint {
                    type Output = BigUint;
                    fn $method(self, rhs: $t) -> BigUint {
                        BigUint(self.0 $op U512::from(rhs as u64))
                    }
                }
            )*
        };
    }

    binop_int!(Add, add, +, u8, u16, u32, u64, usize, i32);
    binop_int!(Sub, sub, -, u8, u16, u32, u64, usize, i32);
    binop_int!(Mul, mul, *, u8, u16, u32, u64, usize, i32);
    binop_int!(Div, div, /, u8, u16, u32, u64, usize, i32);
    binop_int!(Rem, rem, %, u8, u16, u32, u64, usize, i32);

    impl Shl<usize> for BigUint {
        type Output = BigUint;
        fn shl(self, rhs: usize) -> BigUint {
            BigUint(if rhs >= 512 { U512::ZERO } else { self.0 << rhs })
        }
    }
    impl Shl<u32> for BigUint {
        type Output = BigUint;
        fn shl(self, rhs: u32) -> BigUint {
            self << rhs as usize
        }
    }
    impl Shl<usize> for &BigUint {
        type Output = BigUint;
        fn shl(self, rhs: usize) -> BigUint {
            *self << rhs
        }
    }
    impl Shl<u32> for &BigUint {
        type Output = BigUint;
        fn shl(self, rhs: u32) -> BigUint {
            *self << rhs
        }
    }
    impl Shr<usize> for BigUint {
        type Output = BigUint;
        fn shr(self, rhs: usize) -> BigUint {
            BigUint(if rhs >= 512 { U512::ZERO } else { self.0 >> rhs })
        }
    }
    impl Shr<u32> for BigUint {
        type Output = BigUint;
        fn shr(self, rhs: u32) -> BigUint {
            self >> rhs as usize
        }
    }
    impl Shr<usize> for &BigUint {
        type Output = BigUint;
        fn shr(self, rhs: usize) -> BigUint {
            *self >> rhs
        }
    }
    impl Shr<u32> for &BigUint {
        type Output = BigUint;
        fn shr(self, rhs: u32) -> BigUint {
            *self >> rhs
        }
    }

    impl AddAssign for BigUint {
        fn add_assign(&mut self, rhs: BigUint) {
            self.0 = self.0 + rhs.0;
        }
    }
    impl SubAssign for BigUint {
        fn sub_assign(&mut self, rhs: BigUint) {
            self.0 = self.0 - rhs.0;
        }
    }
    impl BitOrAssign for BigUint {
        fn bitor_assign(&mut self, rhs: BigUint) {
            self.0 = self.0 | rhs.0;
        }
    }

    impl fmt::Display for BigUint {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            fmt::Display::fmt(&self.0, f)
        }
    }
    impl fmt::LowerHex for BigUint {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            fmt::LowerHex::fmt(&self.0, f)
        }
    }
}

pub mod num_traits {
    pub trait Zero {
        fn zero() -> Self;
        fn is_zero(&self) -> bool;
    }

    pub trait One {
        fn one() -> Self;
    }

    impl Zero for super::num_bigint::BigUint {
        fn zero() -> Self {
            super::num_bigint::BigUint::zero()
        }

        fn is_zero(&self) -> bool {
            self.is_zero()
        }
    }

    impl One for super::num_bigint::BigUint {
        fn one() -> Self {
            super::num_bigint::BigUint::from(1u32)
        }
    }
}

pub mod rand {
    use core::ops::Range;

    pub trait RngCore {
        fn next_u64(&mut self) -> u64;

        fn fill_bytes(&mut self, dest: &mut [u8]) {
            for chunk in dest.chunks_mut(8) {
                let word = self.next_u64().to_le_bytes();
                chunk.copy_from_slice(&word[..chunk.len()]);
            }
        }
    }

    pub trait SeedableRng {
        fn seed_from_u64(seed: u64) -> Self;
    }

    pub trait RandGen {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self;
    }

    pub trait Rng: RngCore {
        fn gen<T: RandGen>(&mut self) -> T {
            T::gen_from(self)
        }

        fn gen_range(&mut self, range: Range<usize>) -> usize {
            let span = range.end.saturating_sub(range.start);
            if span == 0 {
                return range.start;
            }
            range.start + (self.next_u64() as usize % span)
        }
    }

    impl<T: RngCore + ?Sized> Rng for T {}

    pub mod rngs {
        #[derive(Clone, Copy)]
        pub struct StdRng {
            pub(crate) state: u64,
        }

        pub type ThreadRng = StdRng;
    }

    impl SeedableRng for rngs::StdRng {
        fn seed_from_u64(seed: u64) -> Self {
            let state = if seed == 0 { 0xA5A5_5A5A_D1CE_BA5E } else { seed };
            Self { state }
        }
    }

    impl RngCore for rngs::StdRng {
        fn next_u64(&mut self) -> u64 {
            self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }

    pub fn thread_rng() -> rngs::ThreadRng {
        rngs::StdRng::seed_from_u64(0xC0DE_5EED_5AFE_F00D)
    }

    impl RandGen for u64 {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self {
            rng.next_u64()
        }
    }

    impl RandGen for usize {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self {
            rng.next_u64() as usize
        }
    }

    impl RandGen for u8 {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self {
            rng.next_u64() as u8
        }
    }

    impl RandGen for bool {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self {
            rng.next_u64() & 1 == 1
        }
    }

    impl RandGen for [u8; 32] {
        fn gen_from<R: RngCore + ?Sized>(rng: &mut R) -> Self {
            let mut out = [0u8; 32];
            rng.fill_bytes(&mut out);
            out
        }
    }
}

pub mod rng {
    #[derive(Clone, Copy)]
    pub struct SplitMix64 {
        state: u64,
    }

    impl SplitMix64 {
        #[must_use]
        pub fn seed_from_u64(seed: u64) -> Self {
            Self { state: seed }
        }

        pub fn next_u64(&mut self) -> u64 {
            self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }
}
