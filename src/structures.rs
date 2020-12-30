use crate::fraction::Fraction;
use crate::poly::Polynomial;
use num::Num;

pub trait Field
    : Num +
      std::ops::Neg<Output=Self> +
      std::fmt::Display
{}

pub type Q = Fraction<i128>;

impl Field for Q {}
impl Field for Fraction<Polynomial<Q>> {}
