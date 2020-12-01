use num::{Zero, Num};

use crate::gcd::GCD;

#[derive(Copy, Clone, Debug)]
pub struct Fraction<T : Copy + Num + GCD> {
    num : T,
    den : T
}

impl<T> Fraction<T>
where T : Num + GCD + Copy {
    pub fn new(n : T, d : T) -> Fraction<T> {
        let g = n.gcd(&d);
        Fraction {
            num : n / g,
            den : d / g,
        }
    }
}

impl<T> PartialEq for Fraction<T>
where T : Num + GCD + Copy {
    fn eq(&self, other: &Self) -> bool {
        self.num * other.den == self.den * other.num
    }
}

impl<T> std::ops::Add for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den + other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Sub for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den - other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Mul for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.num,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Div for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den,
            self.den * other.num
        )
    }
}

impl<T> std::convert::From<T> for Fraction<T>
where T : Num + GCD + Copy {
    fn from(x : T) -> Fraction<T> {
        Fraction::new(x, T::one())
    }
}

impl<T> std::ops::Add<T> for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn add(self, other: T) -> Self {
        self + Fraction::from(other)
    }
}

impl<T> std::ops::Sub<T> for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn sub(self, other: T) -> Self {
        self - Fraction::from(other)
    }
}

impl<T> std::ops::Mul<T> for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        self * Fraction::from(other)
    }
}

impl<T> std::ops::Div<T>for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn div(self, other: T) -> Self {
        self / Fraction::from(other)
    }
}

impl<T> num::One for Fraction<T>
where T : Num + GCD + Copy {
    fn one() -> Self {
        Fraction::from(T::one())
    }
}

impl<T> num::Zero for Fraction<T>
where T : Num + GCD + Copy {
    fn zero() -> Self {
        Fraction::from(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.num == T::zero()
    }
}

impl<T> std::ops::Rem for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn rem(self, _other:Self) -> Self {
        Fraction::new(T::zero(), T::one())
    }
}

impl<T> std::ops::Neg for Fraction<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn neg(self) -> Self {
        Fraction::zero() - self
    }
}


impl<T> Num for Fraction<T>
where T : Num + GCD + Copy {
    type FromStrRadixErr = ();

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl<T> std::fmt::Display for Fraction<T>
where T : std::fmt::Display + Num + GCD + Copy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}
