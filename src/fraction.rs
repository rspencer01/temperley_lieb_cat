//! Fractions of arbitrary types

use crate::structures::{Domain, Field, GCDDomain, Ring, RingOps, Signed};
use std::fmt::{Debug, Display};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::gcd::PartialGCD;
use crate::serial::Serialisable;
use crate::tex::Tex;

/// A fraction of `T`s
///
/// Fractions, as constructed by `Fraction::new` are not always stored in lowest terms,
/// but their denominators will always be positive.
#[derive(Copy, Clone, Debug)]
pub struct Fraction<T: Domain> {
    num: T,
    den: T,
}

impl<T: GCDDomain + Signed> Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    pub fn new(n: T, d: T) -> Fraction<T> {
        if n.is_zero() {
            Fraction {
                num: T::zero(),
                den: T::one(),
            }
        } else if d.is_negative() {
            Fraction::new(-n, -d)
        } else {
            if n.is_small() && d.is_small() {
                Fraction { num: n, den: d }
            } else {
                let g = n.partial_gcd(&d);
                let num = &n / &g;
                let den = &d / &g;
                Fraction { num, den }
            }
        }
    }

    #[inline(always)]
    pub fn num(&self) -> &T {
        &self.num
    }

    #[inline(always)]
    pub fn den(&self) -> &T {
        &self.den
    }
}

impl<T: GCDDomain + Signed> Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T> + Rem<&'r T, Output = T>,
{
    pub fn is_integral(&self) -> bool {
        (self.num() % self.den()).is_zero()
    }
}

/*
impl<T> Rem<&Fraction<T>> for &Fraction<T>
where T : Clone + PartialGCD + Signed + RingOps<T,T> + std::fmt::Display,
      for <'r> &'r T : RingOps<&'r T, T> + Rem<&'r T, Output=T> {
    type Output = Fraction<T>;
    /// Reduces the numerator and denominator modulo some element
    ///
    /// WARNING: It is not true that if fractions `a/b = c/d` then
    /// `(a/b).mod(e) = (c/d).mod(e)`.  For this check use `equal_mod`.
    fn rem(self, other: &Fraction<T>) -> Fraction<T> {
        assert!(other.is_integral());
        let reduced_num = self.num() % other.num();
        let reduced_den = self.den() % other.num();
//        if reduced_den.is_zero() {
//            assert!(reduced_num.is_zero(), "Cannot take remainder of fraction");
//            Fraction::new(
//                self.num() / &other,
//                self.den() / &other,
//            ).rem(other)
//        } else {
            Fraction::new(
                reduced_num,
                reduced_den,
            )
//        }
    }
}
*/

impl<T: GCDDomain + Signed> PartialEq for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn eq(&self, other: &Fraction<T>) -> bool {
        self.num() * other.den() == self.den() * other.num()
    }
}

impl<T: GCDDomain + Signed> Eq for Fraction<T> where for<'r> &'r T: RingOps<&'r T, T> {}

impl<T: GCDDomain + Signed> Add for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn add(self, other: Fraction<T>) -> Fraction<T> {
        &self + &other
    }
}

impl<T: GCDDomain + Signed> Add for &Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn add(self, other: &Fraction<T>) -> Fraction<T> {
        Fraction::new(
            self.num() * other.den() + other.num() * self.den(),
            self.den() * other.den(),
        )
    }
}

impl<T: GCDDomain + Signed> Sub for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn sub(self, other: Fraction<T>) -> Fraction<T> {
        &self - &other
    }
}

impl<T: GCDDomain + Signed> Sub for &Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn sub(self, other: &Fraction<T>) -> Fraction<T> {
        Fraction::new(
            self.num() * other.den() - other.num() * self.den(),
            self.den() * other.den(),
        )
    }
}

impl<T: GCDDomain + Signed,> Mul for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn mul(self, other: Fraction<T>) -> Fraction<T> {
        &self * &other
    }
}

impl<T: GCDDomain + Signed> Mul for &Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn mul(self, other: &Fraction<T>) -> Fraction<T> {
        Fraction::new(self.num() * other.num(), self.den() * other.den())
    }
}

impl<T: GCDDomain + Signed> Mul<&Fraction<T>> for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn mul(self, other: &Fraction<T>) -> Fraction<T> {
        &self * other
    }
}

impl<T: GCDDomain + Signed> Div for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn div(self, other: Fraction<T>) -> Fraction<T> {
        &self / &other
    }
}

impl<T: GCDDomain + Signed> Div for &Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn div(self, other: &Fraction<T>) -> Fraction<T> {
        Fraction::new(self.num() * other.den(), self.den() * other.num())
    }
}

impl<T: Domain> From<T> for Fraction<T> {
    fn from(x: T) -> Fraction<T> {
        Fraction {
            num: x,
            den: T::one(),
        }
    }
}

//impl<'a, T:'static> Add<T> for Fraction<T>
//where T : Clone + Num + PartialGCD + Signed + Clone,
//      &'a T : MyNum<T> {
//    type Output = Self;
//
//    fn add(self, other: T) -> Self {
//        self + Fraction::from(other)
//    }
//}
//
//impl<'a, T:'static> Sub<T> for Fraction<T>
//where T : Clone + Num + PartialGCD + Signed + Clone,
//      &'a T : MyNum<T> {
//    type Output = Self;
//
//    fn sub(self, other: T) -> Self {
//        self - Fraction::from(other)
//    }
//}
//
impl<T: GCDDomain + Signed> Mul<T> for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Fraction::new(self.num * other, self.den)
    }
}

//impl<'a, T:'static> Div<T> for Fraction<T>
//where T : Clone + Num + PartialGCD + Signed + Clone,
//      &'a T : MyNum<T> {
//    type Output = Self;
//
//    fn div(self, other: T) -> Self {
//        Fraction::new(self.num, self.den * other)
//    }
//}

impl<T: Domain> Rem for Fraction<T> {
    type Output = Self;

    fn rem(self, _other: Self) -> Self {
        unimplemented!("Use remainder by reference");
    }
}

impl<T: GCDDomain + Signed> Neg for Fraction<T> {
    type Output = Fraction<T>;

    fn neg(self) -> Fraction<T> {
        Fraction {
            num: -self.num,
            den: self.den,
        }
    }
}

impl<T: GCDDomain + Signed> Neg for &Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Fraction<T>;

    fn neg(self) -> Fraction<T> {
        Fraction {
            num: -self.num().clone(),
            den: self.den().clone(),
        }
    }
}

impl<T: Domain + Display> Display for Fraction<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}

impl<T: Domain + Tex> Tex for Fraction<T> {
    fn into_tex(&self) -> String {
        if self.den.is_one() {
            self.num.into_tex()
        } else {
            format!(
                "\\frac{{ {} }}{{ {} }}",
                self.num.into_tex(),
                self.den.into_tex()
            )
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

impl<T: GCDDomain + Signed> Signed for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn is_positive(&self) -> bool {
        self.num.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.num.is_negative()
    }

    fn abs(&self) -> Fraction<T> {
        if self.is_positive() {
            self.clone()
        } else {
            -self.clone()
        }
    }
}

impl<T: GCDDomain> PartialGCD for Fraction<T> {
    fn partial_gcd(&self, other: &Self) -> Self {
        (*other).clone()
    }

    fn is_small(&self) -> bool {
        false
    }
}

impl<T: GCDDomain + Signed + Serialisable> Serialisable for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn serialise(&self) -> String {
        if self.den.is_one() {
            format!("({})", self.num.serialise())
        } else {
            format!("({})/({})", self.num.serialise(), self.den.serialise())
        }
    }

    fn deserialise(input: &str) -> Self {
        let mut parts = Vec::new();
        assert!(input.bytes().nth(0) == Some('(' as u8));
        let mut i = 0;
        let mut c = 1;
        while c > 0 {
            if input.bytes().nth(i) == None {
                panic!("Cannot read fraction");
            }
            i += 1;
            if input.bytes().nth(i) == Some('(' as u8) {
                c += 1;
            }
            if input.bytes().nth(i) == Some(')' as u8) {
                c -= 1;
            }
        }
        parts.push(T::deserialise(&input[1..i]));
        i += 1;
        if input.len() <= i {
            parts.push(T::one());
        } else {
            assert!(input.bytes().nth(i) == Some('/' as u8));
            i += 1;
            assert!(input.bytes().nth(i) == Some('(' as u8));
            let j = i + 1;
            let mut c = 1;
            while c > 0 {
                if input.bytes().nth(i) == None {
                    panic!("Cannot read fraction");
                }
                i += 1;
                if input.bytes().nth(i) == Some('(' as u8) {
                    c += 1;
                }
                if input.bytes().nth(i) == Some(')' as u8) {
                    c -= 1;
                }
            }
            parts.push(T::deserialise(&input[j..i]));
        }
        let den = parts.pop().unwrap();
        let num = parts.pop().unwrap();
        Fraction::new(num, den)
    }
}

impl<T: GCDDomain + Signed> RingOps<Fraction<T>, Fraction<T>> for Fraction<T> where
    for<'r> &'r T: RingOps<&'r T, T>
{
}

impl<T: GCDDomain + Signed> Ring for Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn zero() -> Self {
        Fraction::from(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.num == T::zero()
    }
    fn one() -> Self {
        Fraction::from(T::one())
    }

    fn is_one(&self) -> bool {
        self.num == self.den
    }
}

impl<T: GCDDomain + Signed> Domain for Fraction<T> where for<'r> &'r T: RingOps<&'r T, T> {}

impl<T: GCDDomain + Signed> Field for Fraction<T> where for<'r> &'r T: RingOps<&'r T, T> {}

impl<'a, T: 'a+ GCDDomain + Signed,> RingOps<&'a Fraction<T>, Fraction<T>> for &'a Fraction<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic_arithmetic() {
        assert_eq!(
            Fraction::new(1, 2) + Fraction::new(1, 2),
            Fraction::new(1, 1)
        );
        assert_eq!(
            Fraction::new(1, 3) - Fraction::new(1, 4),
            Fraction::new(1, 12)
        );
        assert_eq!(
            Fraction::new(1, 3) * Fraction::new(3, 1),
            Fraction::new(1, 1)
        );
    }

    #[test]
    fn is_integral() {
        assert!(Fraction::new(1, 1).is_integral());
        assert!(Fraction::new(2, 1).is_integral());
        assert!(Fraction::new(4, 2).is_integral());
        assert!(Fraction::new(100, 5).is_integral());
        assert!(!Fraction::new(101, 5).is_integral());
        assert!(Fraction::new(-10, 5).is_integral());
    }
}
