use num::Num;
use std::ops::{Add, Sub, Mul, Div, Neg, Rem};
use std::fmt::{Debug, Display};

use crate::gcd::PartialGCD;
use crate::num::{One, Zero, Signed};
use crate::tex::Tex;
use crate::serial::Serialisable;
use crate::structures::Field;
use crate::structures;
use crate::fraction::Fraction;

const MAX_DEGREE : usize = 200;

/// A polynomial in a single variable over a _field_.
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Polynomial<T : Copy> {
    coeffs : [T; MAX_DEGREE],
    degree : usize,
}

impl<T> Polynomial<T>
where T : Copy {
    fn degree(&self) -> usize {
        self.degree
    }

    fn hightest_term(&self) -> T {
        self.coeffs[self.degree()]
    }
}

impl<T> Polynomial<T>
where T : Copy + Zero {
    pub fn new<V>(in_coeffs: &[V]) -> Polynomial<T>
    where V : Into<T> + Copy {
        let mut coeffs = [T::zero() ; MAX_DEGREE];
        let mut degree = 0;
        let mut it = in_coeffs.iter();
        for i in 0..MAX_DEGREE {
            if let Some(x) = it.next() {
                coeffs[i] = (*x).into();
                if !coeffs[i].is_zero() {
                    degree = i;
                }
            } else {
                break;
            }
       }
       assert!(it.next().is_none(), "Overflow in polynomial construction");
       Polynomial {
           coeffs,
           degree,
       }
    }

    fn shift(&mut self, n : usize) {
        for i in (n..MAX_DEGREE).rev() {
            self.coeffs[i] = self.coeffs[i-n];
        }
        for i in 0..n {
            self.coeffs[i] = T::zero();
        }
        self.degree += n;
    }
}

impl<T> Polynomial<T>
where T : Copy + Zero + One {
    pub fn gen() -> Self {
        Polynomial::new(&[T::zero(), T::one()])
    }
}


impl<T> Polynomial<T>
where T : Copy + Field + Eq {
    fn div_mod(&self, other: &Self) -> (Self, Self) {
        assert!(!other.is_zero(), "Dividing polynomial by zero");
        if self.is_zero() {
            return (Polynomial::zero(), Polynomial::zero());
        }
        if *self == *other {
            return (Polynomial::one(), Polynomial::zero());
        }
        if *self == -*other {
            return (-Polynomial::one(), Polynomial::zero());
        }
        if self.degree() < other.degree() {
            return (Polynomial::zero(), *self);
        }
        let mut rem = *self;
        let mut ans = [T::zero() ; MAX_DEGREE];
        for i in (other.degree()..self.degree()+1).rev() {
            if rem.coeffs[i].is_zero() {
                continue;
            }
            // Because T is guaranteed to be a field, this gives no remainder
            let m = rem.hightest_term() / other.hightest_term();
            let mut diff = *other * m;
            diff.shift(i-other.degree());
            rem = rem - diff;
            ans[i - other.degree] = m;
            assert!(rem.degree() < i || i == 0);
        }
        (Polynomial::new(&ans), rem)
    }
}

impl<T> Add for Polynomial<T>
where T : Copy + Add + Zero {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut ans = self.coeffs;
        for i in 0..other.degree()+1 {
            ans[i] = ans[i] + other.coeffs[i];
        }
        let mut degree = 0;
        for i in 0..self.degree().max(other.degree())+1 {
            if !ans[i].is_zero() {
                degree = i;
            }
        }
        Polynomial{
            coeffs: ans,
            degree,
        }
    }
}

impl<T> Neg for Polynomial<T>
where T : Copy + Zero + Neg<Output=T> {
    type Output = Self;

    fn neg(self) -> Self {
        let mut ans = [T::zero(); MAX_DEGREE];
        for i in 0..self.degree()+1 { ans[i] = -self.coeffs[i] }
        Polynomial{
            coeffs:ans,
            degree:self.degree(),
        }
    }
}

impl<T> Sub for Polynomial<T>
where T : Copy + Zero + Sub<Output=T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut ans = self.coeffs;
        for i in 0..other.degree()+1 {
            ans[i] = ans[i] - other.coeffs[i];
        }
        Polynomial::new(&ans)
    }
}

impl<T> Mul for Polynomial<T>
where T : Copy + Zero + Add + Mul<Output=T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        assert!(self.degree() + other.degree() < MAX_DEGREE, "Attempt to multiply polynomials with overflow");
        let mut ans = [T::zero() ; MAX_DEGREE];
        let max_resulting_degree = self.degree() + other.degree();
        for i in 0..max_resulting_degree + 1 {
            for j in i.saturating_sub(other.degree())..i.min(self.degree())+1 {
                    ans[i] = ans[i] + self.coeffs[j] * other.coeffs[i-j];
            }
        }
        Polynomial::new(&ans)
    }
}

impl<T> Mul<T> for Polynomial<T>
where T : Copy + One + Mul + Eq {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        if other.is_one() {
            return self
        }
        let mut ans = self.coeffs;
        for i in 0..self.degree()+1 {
            ans[i] = ans[i] * other;
        }
        Polynomial{
            coeffs:ans,
            degree:self.degree()
        }
    }
}

impl<T> Div<T> for Polynomial<T>
where T : Copy + Zero + Div<Output=T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        assert!(!other.is_zero(), "Cannot divide polynomial by zero");
        let mut ans = self.coeffs;
        for i in 0..self.degree()+1 {
            ans[i] = ans[i] / other;
        }
        Polynomial{
            coeffs:ans,
            degree:self.degree()
        }
    }
}

impl<T> Div for Polynomial<T>
where T : Field + Copy + Eq {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.div_mod(&other).0
    }
}

impl<T> One for Polynomial<T>
where T : Copy + One + Zero {
    fn one() -> Self {
        Polynomial::new(&[T::one()])
    }
}

impl<T> Zero for Polynomial<T>
where T : Copy + Zero {
    fn zero() -> Self {
        let c : [T;0] = [];
        Polynomial::new(&c)
    }

    fn is_zero(&self) -> bool {
        self.degree() == 0 && self.hightest_term().is_zero()
    }
}

impl<T> Rem for Polynomial<T>
where T : Copy + Field + Eq {
    type Output = Self;

    fn rem(self, other:Self) -> Self {
        self.div_mod(&other).1
    }
}


impl<T> Num for Polynomial<T>
where T : Copy + Field + Eq {
    type FromStrRadixErr = ();

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl<T> Signed for Polynomial<T>
where T : Copy + Field + Eq + Signed {
    fn is_positive(&self) -> bool {
        self.hightest_term().is_positive()
    }

    fn is_negative(&self) -> bool {
        self.hightest_term().is_negative()
    }

    fn abs(&self) -> Polynomial<T> {
        if self.is_positive() {
            self.clone()
        } else {
            -self.clone()
        }
    }

    fn abs_sub(&self, _ : &Polynomial<T>) -> Polynomial<T> {
        unimplemented!()
    }

    fn signum(&self) -> Polynomial<T> {
        unimplemented!()
    }
}

impl PartialGCD for Polynomial<structures::Q> {
    fn partial_gcd(&self, other: &Polynomial<structures::Q>) -> Polynomial<structures::Q> {
        if (self.coeffs.iter().fold(0, |a, x| x.den().abs().max(a)) > 0xffffffff) ||
           (self.coeffs.iter().fold(0, |a, x| x.num().abs().max(a)) > 0xffffffff) ||
           (other.coeffs.iter().fold(0, |a, x| x.den().abs().max(a)) > 0xfffffffff) ||
           (other.coeffs.iter().fold(0, |a, x| x.num().abs().max(a)) > 0xfffffffff) {
            return Polynomial::one();
        }
        if self.degree() < other.degree() {
            return other.partial_gcd(self);
        }
        if other.is_zero() {
            return *self;
        }
        let div_rem = self.div_mod(other);
        if div_rem.0.is_zero() {
            Polynomial::one()
        } else {
            let mut re = div_rem.1;
            //let c = re.coeffs.iter().fold(re.hightest_term().den(), |x,y | x * y.den() / x.gcd(&y.den()));
            let c = re.hightest_term();
            if !c.is_zero() {
                re = re / c;
            }
            let sub = other.partial_gcd(&re);
            sub / sub.hightest_term()
        }
    }
}

impl<T> Display for Polynomial<T>
where T: Copy + One + Zero + Display + Eq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut s = false;
        for i in 0..self.degree().max(0)+1 {
            if !self.coeffs[i].is_zero() || (i == 0 && self.degree() == 0) {
                if s {
                    write!(f, " + ")?;
                }
                s = true;
                if !self.coeffs[i].is_one() || i == 0 {
                    write!(f, "{}", self.coeffs[i])?;
                }
                if i > 1 {
                    write!(f, "X^{}", i)?;
                } else if i == 1 {
                    write!(f, "X")?;
                }
            }
        }
        Ok(())
    }
}

impl<T> Debug for Polynomial<T>
where T: Copy + Display + One + Zero + Eq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl<T> Tex for Polynomial<T>
where T: Copy + Zero + One + Eq + Signed + Tex {
    fn into_tex(&self) -> String {
        let mut ans = String::new();
        let mut s = false;
        for i in 0..self.degree().max(0)+1 {
            if !self.coeffs[i].is_zero() || (i == 0 && self.degree() == 0) {
                if s {
                    if self.coeffs[i].is_positive() {
                        ans += " + ";
                    }
                }
                s = true;
                if !self.coeffs[i].abs().is_one() || i == 0 {
                    if self.coeffs[i].is_multiterm() {
                        ans += "\\left(";
                    }
                    ans += &self.coeffs[i].into_tex();
                    if self.coeffs[i].is_multiterm() {
                        ans += "\\right)";
                    }
                }
                if self.coeffs[i].is_negative() && (-self.coeffs[i]).is_one() && i > 0 {
                    ans += " - ";
                }
                if i > 1 {
                    ans += &format!("X^{{ {} }}", i);
                } else if i == 1 {
                    ans += "X";
                }
            }
        }
        ans
    }

    fn is_multiterm(&self) -> bool {
        self.coeffs.iter().filter(|x| !x.is_zero()).count() > 1
    }
}

pub fn quantum(n : i128) -> Polynomial<structures::Q> {
    if n == 0 {
        return Polynomial::zero();
    }
    fn choose(n : i128, r: i128) -> i128 {
        let mut ans = 1;
        for i in n-r+1..n+1 {
            ans = ans * i;
        }
        for i in 1..r+1 {
            ans = ans / i;
        }
        ans
    }
    let mut coeffs = Vec::new();
    for _ in 0..n {
        coeffs.push(Fraction::zero())
    }
    for i in 0..(n-1)/2+1 {
        coeffs[(n-1-2*i) as usize] = Fraction::from(
            (1-2*(i%2)) * choose(n-1-i,i)
        );
    }
    Polynomial::new(coeffs.as_slice())
}

impl Serialisable for Polynomial<structures::Q> {
    fn serialise(&self) -> String {
        self.coeffs[..self.degree+1].to_vec().serialise()
    }

    fn deserialise(inpt : &str) -> Self {
        let res : Vec<structures::Q> = Vec::deserialise(inpt);
        Polynomial::new(
            res.as_slice()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::structures::Q;

    #[test]
    fn new () {
        assert_eq!(Polynomial::<Q>::new(&[1,2,3,4,0,0,0]).degree(), 3);
        assert_eq!(Polynomial::<Q>::new(&[1,2,3,4,0,0,0]),
                         Polynomial::new(&[1,2,3,4])
            );
        assert_eq!(Polynomial::<Q>::new(&[0,0,0]).degree(), 0);
    }

    #[test]
    fn gen() {
        assert_eq!(Polynomial::<Q>::gen(), Polynomial::new(&[Q::zero(),Q::one()]));
        assert_eq!(Polynomial::<Q>::gen().degree(), 1);
        assert_eq!("X", format!("{}", Polynomial::<Q>::gen()));
        assert_ne!(Polynomial::<Q>::gen(), Polynomial::zero());
        assert_ne!(Polynomial::<Q>::gen(), Polynomial::one());
    }

    #[test]
    fn additive() {
        let a = Polynomial::<Q>::new(&[1,2,3]);
        let b = Polynomial::<Q>::new(&[1,2,3]);
        assert_eq!(Polynomial::zero(), a - b);
    }

    #[test]
    fn neg() {
        assert_eq!(-Polynomial::<Q>::new(&[1,2,3]), Polynomial::new(&[-1,-2,-3]));
    }

    #[test]
    fn sub() {
        assert_eq!(Polynomial::<Q>::new(&[1,2,3]) - Polynomial::new(&[2,0,1,5]),
                   Polynomial::<Q>::new(&[-1,2,2,-5]));
    }

    #[test]
    fn mul() {
        let a = Polynomial::<Q>::new(&[1,2,3]);
        let b = Polynomial::<Q>::new(&[0,4,6]);
        assert_eq!(Polynomial::new(&[0,4,14,24,18]), a * b);
    }

    #[test]
    fn div() {
        let a = Polynomial::<Q>::new(&[0, 0, -55, 0, 1320, -176, -10560, 3520, 38720, -19712, -70400, 45056, 61440, -45056, -20480, 16384]);
        let b = Polynomial::new(&[0,5,0,-20,16]);
        assert_eq!(Polynomial::<Q>::new(&[0, -11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024])
                , a / b);
    }

    #[test]
    fn quantum_num() {
        assert_eq!(quantum(0), Polynomial::zero());
        assert_eq!(quantum(1), Polynomial::one());
        assert_eq!(quantum(2), Polynomial::gen());
        assert_eq!(quantum(3), Polynomial::new(&[-1,0,1]));
        assert_eq!(quantum(8), Polynomial::new(&[0,-4,0,10,0,-6,0,1]));
    }

    #[test]
    fn gcd() {
        assert_eq!(Polynomial::new(&[1,2,1]).partial_gcd(&Polynomial::new(&[1,1])) , Polynomial::new(&[1,1]));
        assert_eq!(Polynomial::new(&[2,4,2]).partial_gcd(&Polynomial::new(&[4,4])) , Polynomial::new(&[1,1]));
        assert_eq!(quantum(17).partial_gcd(&quantum(15)) , Polynomial::one());
        assert_eq!(quantum(17).partial_gcd(&quantum(15)) , Polynomial::one());
        assert_eq!(quantum(27).partial_gcd(&quantum(3)) , quantum(3));
        assert_eq!(quantum(27).partial_gcd(&quantum(6)) , quantum(3));
    }

    #[test]
    fn big_gcd() {
        let a = Polynomial::<Q>::new(&[
            24,
            0,
            -524,
            0,
            5056,
            0,
            -27152,
            0,
            84534,
            0,
            -138895,
            0,
            23826,
            0,
            445831,
            0,
            -1177984,
            0,
            1742628,
            0,
            -1768127,
            0,
            1310596,
            0,
            -728166,
            0,
            305827,
            0,
            -96767,
            0,
            22722,
            0,
            -3841,
            0,
            442,
            0,
            -31,
            0,
            1
        ]);

        let b = Polynomial::<Q>::new(&[
            Fraction::new(304,5),
            Fraction::new(0,1),
            Fraction::new(21856,-15),
            Fraction::new(0,1),
            Fraction::new(238652,15),
            Fraction::new(0,1),
            Fraction::new(1527248,-15),
            Fraction::new(0,1),
            Fraction::new(1258064,3),
            Fraction::new(0,1),
            Fraction::new(17683936,-15),
            Fraction::new(0,1),
            Fraction::new(35414599,15),
            Fraction::new(0,1),
            Fraction::new(-52134644,15),
            Fraction::new(0,1),
            Fraction::new(57703364,15),
            Fraction::new(0,1),
            Fraction::new(48757748,-15),
            Fraction::new(0,1),
            Fraction::new(10579261,5),
            Fraction::new(0,1),
            Fraction::new(-1064413,1),
            Fraction::new(0,1),
            Fraction::new(6190517,15),
            Fraction::new(0,1),
            Fraction::new(-366422,3),
            Fraction::new(0,1),
            Fraction::new(406151,15),
            Fraction::new(0,1),
            Fraction::new(65266,-15),
            Fraction::new(0,1),
            Fraction::new(7177,15),
            Fraction::new(0,1),
            Fraction::new(161,-5),
            Fraction::new(0,1),
            Fraction::new(1,1)
        ]);
        print!("{}", a.partial_gcd(&b));
    }
}
