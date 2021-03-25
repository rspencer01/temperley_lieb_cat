//! Algebraic structures and definitions of numerical operations

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

impl NumOps<i128, i128> for i128 {}
impl NumOps<&i128, i128> for &i128 {}
