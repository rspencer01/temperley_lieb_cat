//! Algebraic structures and definitions of numerical operations
//!
//! Really, this should be handled by the `Num` crate, but it assumes
//! silly things like remainders for all "numbers".

use crate::fraction::Fraction;
use crate::gcd::PartialGCD;
use std::ops::{Add, Sub, Mul, Div, Neg};

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
pub trait Ring : RingOps<Self, Self> + Eq + Clone {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn is_one(&self) -> bool;
}

/// A marker trait to assert that this ring is a domain
///
/// Domains have no zero divisors
pub trait Domain : Ring {}

/// A marker trait to assert that this domain is a GCD domain
///
/// A GCD domain is one for which each pair of elements has
/// a greatest common divisor.  Such domains have a well defined
/// concept of division with remainder.
///
/// Implimentations of [GCDDomain] do not need to be able to find
/// the greatest common divisor, and hence only need implement [gcd::PartialGCD].
pub trait GCDDomain : Domain + PartialGCD {}

/// A marker trait to assert that this domain is a field
///
/// Fields have well defined inverses.
pub trait Field : Domain {}

impl RingOps<i128, i128> for i128 {}
impl Signed for i128 {
    fn is_positive(&self) -> bool {
        *self > 0
    }
    fn is_negative(&self) -> bool {
        *self < 0
    }
    fn abs(&self) -> Self {
        i128::abs(*self)
    }
}
impl Ring for i128 {
    fn zero() -> Self { 0 }
    fn is_zero(&self) -> bool { *self == 0 }
    fn one() -> Self { 1 }
    fn is_one(&self) -> bool { *self == 1 }
}
impl Domain for i128 {}
impl GCDDomain for i128 {}

impl RingOps<&i128, i128> for &i128 {}
