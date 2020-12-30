use num::{Num, Signed, One, Zero};
use std::ops::{Add, Sub, Mul, Div, Rem, Neg};
use std::fmt::{Debug, Display};

use crate::gcd::PartialGCD;
use crate::tex::Tex;
use crate::serial::Serialisable;

#[derive(Copy, Clone, Debug)]
pub struct Fraction<T : Copy + Num + PartialGCD> {
    num : T,
    den : T
}

impl<T> Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    pub fn new(n : T, d : T) -> Fraction<T> {
        if n.is_zero() {
            Fraction {
                num: T::zero(),
                den: T::one(),
            }
        } else if d.is_negative() {
            Fraction::new(-n, -d)
        } else {
            let g = n.partial_gcd(&d);
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
where T : Num + PartialGCD + Copy + Signed {
    fn eq(&self, other: &Self) -> bool {
        self.num * other.den == self.den * other.num
    }
}

impl<T> Eq for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
}

impl<T> Add for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den + other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> Sub for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Fraction::new(
            self.num * other.den - other.num * self.den,
            self.den * other.den
        )
    }
}

impl<T> Mul for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let g1 = self.num.partial_gcd(&other.den);
        let g2 = other.num.partial_gcd(&self.den);
        Fraction::new(
            (self.num / g1) * (other.num / g2),
            (self.den / g2) * (other.den / g1)
        )
    }
}

impl<T> Div for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let g1 = self.num.partial_gcd(&other.num);
        let g2 = other.den.partial_gcd(&self.den);
        Fraction::new(
            (self.num / g1) * (other.den / g2),
            (self.den / g2) * (other.num / g1)
        )
    }
}

impl<T> std::convert::From<T> for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    fn from(x : T) -> Fraction<T> {
        Fraction {
            num : x,
            den : T::one()
        }
    }
}

impl<T> Add<T> for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn add(self, other: T) -> Self {
        self + Fraction::from(other)
    }
}

impl<T> Sub<T> for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn sub(self, other: T) -> Self {
        self - Fraction::from(other)
    }
}

impl<T> Mul<T> for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Fraction::new(self.num * other, self.den)
    }
}

impl<T> Div<T> for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Fraction::new(self.num, self.den * other)
    }
}

impl<T> One for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    fn one() -> Self {
        Fraction::from(T::one())
    }
}

impl<T> Zero for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    fn zero() -> Self {
        Fraction::from(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.num == T::zero()
    }
}

impl<T> Rem for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn rem(self, _other:Self) -> Self {
        unimplemented!()
    }
}

impl<T> Neg for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type Output = Self;

    fn neg(self) -> Self {
        Fraction::new(-self.num, self.den)
    }
}


impl<T> Num for Fraction<T>
where T : Num + PartialGCD + Copy + Signed {
    type FromStrRadixErr = ();

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl<T> Display for Fraction<T>
where T : Display + Num + PartialGCD + Copy + Signed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}

impl<T> Tex for Fraction<T>
where T : Copy + Num + PartialGCD + Tex + Signed {
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
where T : Copy + Num + PartialGCD + Tex + Signed {
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

impl<T> Serialisable for Fraction<T>
where T: Serialisable + PartialGCD + Signed + Copy {

    fn serialise(&self) -> String {
        format!("({})/({})", self.num.serialise(), self.den.serialise())
    }

    fn deserialise(input : &str) -> Self {
        let mut parts : Vec<&str> = Vec::new();
        assert!(input.bytes().nth(0) == Some('(' as u8));
        let mut i = 0;
        let mut c = 1;
        while c > 0 {
            if input.bytes().nth(i) == None { panic!("Cannot read fraction"); }
            i += 1;
            if input.bytes().nth(i) == Some('(' as u8) { c += 1; }
            if input.bytes().nth(i) == Some(')' as u8) { c -= 1; }
        }
        i += 1;
        assert!(input.bytes().nth(i) == Some('/' as u8));
        parts.push(&input[1..i-1]);
        i += 1;
        assert!(input.bytes().nth(i) == Some('(' as u8));
        let j = i+1;
        let mut c = 1;
        while c > 0 {
            if input.bytes().nth(i) == None { panic!("Cannot read fraction"); }
            i += 1;
            if input.bytes().nth(i) == Some('(' as u8) { c += 1; }
            if input.bytes().nth(i) == Some(')' as u8) { c -= 1; }
        }
        parts.push(&input[j..i]);
        Fraction::new(
            T::deserialise(parts[0]),
            T::deserialise(parts[1]),
        )
    }
}
