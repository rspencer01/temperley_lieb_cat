//! Algebraic structures and definitions of numerical operations
//!
//! Really, this should be handled by the `Num` crate, but it assumes
//! silly things like remainders for all "numbers".

use num::{Zero, One};
use crate::fraction::Fraction;
use std::ops::{Add, Sub, Mul, Div, Neg};

/// Rational numbers to 128 bits of accuracy
pub type Q = Fraction<i128>;

/// A trait describing numerical operations `+ - * /` on a type
/// with another type
///
/// You can use this to bind your types, for example `A : NumOps<Rhs=A, Out=A>` for an algebra
/// or `A : NumOps<Rhs=B, Out=B>` for a ring homomorphism $A \to B$.
pub trait NumOps<Rhs, Out>:
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
pub trait Ring : NumOps<Self, Self> + Zero + One + Eq + Clone {}

/// A marker trait to assert that this ring is a domain
///
/// Domains have no zero divisors
pub trait Domain : Ring {}

impl NumOps<i128, i128> for i128 {}
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
impl Ring for i128 {}
impl Domain for i128 {}

impl NumOps<&i128, i128> for &i128 {}
