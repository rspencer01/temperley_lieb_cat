//! Algebraic structures and definitions of numerical operations
//!
//! Really, this should be handled by the `Num` crate, but it assumes
//! silly things like remainders for all "numbers".
//!
//! In general, you should not need to interface with the traits defined in
//! this module, except to implement them on any new rings you define.
//!
//! ## Example
//! Here we define a basic implementation of the field with two elements
//! ```
//! # use temperley_lieb_cat::structures::*;
//! # use std::ops::*;
//! #[derive(Clone, PartialEq, Eq)]
//! struct GF2(bool);
//! impl Add for GF2 {
//!     // ...
//! #    type Output = GF2;
//! #    fn add(self, other : GF2) -> GF2 {
//! #        GF2(self.0 ^ other.0)
//! #    }
//! }
//! impl Mul for GF2 {
//!     // ...
//! #    type Output = GF2;
//! #    fn mul(self, other : GF2) -> GF2 {
//! #        GF2(self.0 & other.0)
//! #    }
//! }
//! impl Sub for GF2 {
//!     // ...
//! #    type Output = GF2;
//! #    fn sub(self, other : GF2) -> GF2 {
//! #        GF2(self.0 ^ other.0)
//! #    }
//! }
//! impl Div for GF2 {
//!     // ...
//! #    type Output = GF2;
//! #    fn div(self, other : GF2) -> GF2 {
//! #        GF2(self.0 & other.0)
//! #    }
//! }
//! impl Neg for GF2 {
//!     // ...
//! #    type Output = GF2;
//! #    fn neg(self) -> GF2 {
//! #        self
//! #    }
//! }
//!
//! // Since we have defined all necessary operations, we can mark this type
//! // with the `RingOps` makrer.
//! impl RingOps<GF2,GF2> for GF2 {};
//!
//! // We need to define the units to give this a ring structure
//! impl Ring for GF2 {
//!     fn zero() -> Self {
//!         GF2(false)
//!     }
//!     fn is_zero(&self) -> bool {
//!         !self.0
//!     }
//!     fn one() -> Self {
//!         GF2(true)
//!     }
//!     fn is_one(&self) -> bool {
//!         self.0
//!     }
//! };
//!
//! // This happens to be a domain...
//! impl Domain for GF2 {}
//! // ...and a field
//! impl Field for GF2 {}
//! ```

use crate::fraction::Fraction;
use crate::gcd::{GCD, PartialGCD};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Rational numbers to 128 bits of accuracy
pub type Q = Fraction<i128>;

/// A trait describing numerical operations `+ - * /` on a type
/// with another type
///
/// You can use this to bind your types, for example `A : RingOps<Rhs=A, Out=A>` for an algebra
/// or `A : RingOps<Rhs=B, Out=B>` for a ring homomorphism $A \to B$.
///
/// Note that many rings don't define a division operator.  It is acceptable to leave this
/// as unimplemented.
pub trait RingOps<Rhs, Out>:
    Add<Rhs, Output = Out>
    + Sub<Rhs, Output = Out>
    + Mul<Rhs, Output = Out>
    + Div<Rhs, Output = Out>
    + Neg<Output = Out>
{
}

/// A trait for types with well defined concepts of "negative" and "positive"
pub trait Signed: Neg<Output = Self> {
    fn is_positive(&self) -> bool;
    fn is_negative(&self) -> bool;
    fn abs(&self) -> Self;
}

/// A trait for objects that behave like elements of a ring
///
/// These items can be added, subtracted, divided and multiplied.
/// There is a well defined concept of negatives and zero and one.
pub trait Ring: RingOps<Self, Self> + Eq + Clone {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn is_one(&self) -> bool;
}

/// A marker trait to assert that this ring is a domain
///
/// Domains have no zero divisors
pub trait Domain: Ring {}

/// A marker trait to assert that this domain is a GCD domain
///
/// A GCD domain is one for which each pair of elements has
/// a greatest common divisor.  Such domains have a well defined
/// concept of division with remainder.
///
/// Implimentations of [GCDDomain] do not need to be able to find
/// the greatest common divisor, and hence only need implement [PartialGCD].
pub trait GCDDomain: Domain + PartialGCD {}

/// A marker trait to assert that this domain is a field
///
/// Fields have well defined inverses.
pub trait Field: Domain {}

macro_rules! ring_ops_for_integers_impl {
    ($($t:ty)*) => {$(
        impl RingOps<$t, $t> for $t {}
        impl RingOps<&$t, $t> for &$t {}

        impl Signed for $t {
            fn is_positive(&self) -> bool {
                *self > 0
            }
            fn is_negative(&self) -> bool {
                *self < 0
            }
            fn abs(&self) -> Self {
                Self::abs(*self)
            }
        }

        impl Ring for $t {
            fn zero() -> Self {
                0
            }
            fn is_zero(&self) -> bool {
                *self == 0
            }
            fn one() -> Self {
                1
            }
            fn is_one(&self) -> bool {
                *self == 1
            }
        }

        impl GCD for $t {
            fn gcd(&self, other: &Self) -> Self {
                let mut a = self.abs();
                let mut b = other.abs();
                if a < b {
                    let c = b;
                    b = a;
                    a = c;
                }
                while b != 0 {
                    let c = a % b;
                    a = b;
                    b = c;
                }
                a
            }
        }

        impl Domain for $t {}
        impl GCDDomain for $t {}
    )*}
}
ring_ops_for_integers_impl![i128 i64 i32 isize];
