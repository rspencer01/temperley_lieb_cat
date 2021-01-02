use crate::fraction::Fraction;
//use crate::poly::Polynomial;
use num::Num;
use std::ops::{Add, Sub, Mul, Div, Neg};

pub trait Field
    : Num +
      std::ops::Neg<Output=Self> +
      std::fmt::Display
{}

pub type Q = Fraction<i128>;

impl Field for Q {}
//impl Field for Fraction<Polynomial<Q>> {}

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
