use num::{Zero, Num, Signed};

use crate::gcd::GCD;
use crate::tex::Tex;

#[derive(Copy, Clone, Debug)]
pub struct Fraction<T : Copy + Num + GCD + Signed> {
    num : T,
    den : T
}

impl<T> Fraction<T>
where T : Num + GCD + Copy + Signed {
    pub fn new(n : T, d : T) -> Fraction<T> {
        if n.is_zero() {
            Fraction {
                num: T::zero(),
                den: T::one(),
            }
        } else if d.is_negative() {
            Fraction::new(-n, -d)
        } else {
            let g = n.gcd(&d);
            Fraction {
                num : n / g,
                den : d / g,
            }
        }
    }

    pub fn num(&self) -> T {
        self.num
    }

    pub fn den(&self) -> T {
        self.den
    }
}

impl<T> PartialEq for Fraction<T>
where T : Num + GCD + Copy + Signed {
    fn eq(&self, other: &Self) -> bool {
        self.num * other.den == self.den * other.num
    }
}

impl<T> std::ops::Add for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den + other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Sub for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den - other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Mul for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.num,
            self.den * other.den
        )
    }
}

impl<T> std::ops::Div for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den,
            self.den * other.num
        )
    }
}

impl<T> std::convert::From<T> for Fraction<T>
where T : Num + GCD + Copy + Signed {
    fn from(x : T) -> Fraction<T> {
        Fraction::new(x, T::one())
    }
}

impl<T> std::ops::Add<T> for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn add(self, other: T) -> Self {
        self + Fraction::from(other)
    }
}

impl<T> std::ops::Sub<T> for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn sub(self, other: T) -> Self {
        self - Fraction::from(other)
    }
}

impl<T> std::ops::Mul<T> for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        self * Fraction::from(other)
    }
}

impl<T> std::ops::Div<T>for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn div(self, other: T) -> Self {
        self / Fraction::from(other)
    }
}

impl<T> num::One for Fraction<T>
where T : Num + GCD + Copy + Signed {
    fn one() -> Self {
        Fraction::from(T::one())
    }
}

impl<T> num::Zero for Fraction<T>
where T : Num + GCD + Copy + Signed {
    fn zero() -> Self {
        Fraction::from(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.num == T::zero()
    }
}

impl<T> std::ops::Rem for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn rem(self, _other:Self) -> Self {
        Fraction::new(T::zero(), T::one())
    }
}

impl<T> std::ops::Neg for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type Output = Self;

    fn neg(self) -> Self {
        Fraction::zero() - self
    }
}


impl<T> Num for Fraction<T>
where T : Num + GCD + Copy + Signed {
    type FromStrRadixErr = ();

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl<T> std::fmt::Display for Fraction<T>
where T : std::fmt::Display + Num + GCD + Copy + Signed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}

impl<T> Tex for Fraction<T>
where T : Copy + Num + GCD + Tex + Signed {
    fn into_tex(&self) -> String {
        if self.den.is_one() {
            self.num.into_tex()
        } else {
            format!("\\frac{{ {} }}{{ {} }}", self.num.into_tex(), self.den.into_tex())
        }
    }

    fn is_multiterm(&self) -> bool {
        if self.den.is_one() {
            self.num.is_multiterm()
        } else {
            false
        }
    }
}

impl<T> Signed for Fraction<T>
where T : Copy + Num + GCD + Tex + Signed {
    fn is_positive(&self) -> bool {
        self.num.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.num.is_negative()
    }

    fn abs(&self) -> Fraction<T> {
        if self.is_positive() {
            *self
        } else {
            -*self
        }
    }

    fn abs_sub(&self, _ : &Fraction<T>) -> Fraction<T> {
        unimplemented!()
    }

    fn signum(&self) -> Fraction<T> {
        unimplemented!()
    }
}
