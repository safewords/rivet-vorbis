//! Runtime selection of the vector code paths.
//!
//! The kernels are written once, as plain Rust whose loops vectorise, and
//! [`multiversion!`] compiles each twice on x86-64: for the baseline
//! (SSE2) and with AVX2 enabled, picked at run time by CPUID. On aarch64
//! NEON is the baseline, so the one build already uses it. The kernels are
//! integer code, or floating-point code with each sum in a fixed order and
//! no fused multiply-add, so every path gives the same result to the bit;
//! the tests hold each to the plain definition.
//!
//! The `force-scalar` feature compiles the dispatch out, leaving only the
//! baseline build of every kernel (what CI tests the fallback with).

/// Whether the AVX2 builds of the kernels may run.
#[cfg(all(target_arch = "x86_64", not(feature = "force-scalar")))]
#[inline]
pub(crate) fn avx2() -> bool {
    std::is_x86_feature_detected!("avx2")
}

/// `fn name(args) -> ret { body }` becomes that function, which runs the
/// body compiled with AVX2 when the CPU has it and as written otherwise.
macro_rules! multiversion {
    ($(#[$attr:meta])* $vis:vis fn $name:ident($($arg:ident: $ty:ty),* $(,)?) $(-> $ret:ty)? $body:block) => {
        $(#[$attr])*
        $vis fn $name($($arg: $ty),*) $(-> $ret)? {
            #[inline(always)]
            fn generic($($arg: $ty),*) $(-> $ret)? $body
            #[cfg(all(target_arch = "x86_64", not(feature = "force-scalar")))]
            {
                #[target_feature(enable = "avx2")]
                fn avx2($($arg: $ty),*) $(-> $ret)? {
                    generic($($arg),*)
                }
                if $crate::simd::avx2() {
                    // SAFETY: the CPU has AVX2 (checked just above), the
                    // only requirement of a `target_feature` function; the
                    // body itself is safe code.
                    return unsafe { avx2($($arg),*) };
                }
            }
            generic($($arg),*)
        }
    };
}
pub(crate) use multiversion;
