//! Algebraic structures and definitions of numerical operations
//!
//! Really, this should be handled by the `Num` crate, but it assumes
//! silly things like remainders for all "numbers".

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
impl NumOps<&i128, i128> for &i128 {}
