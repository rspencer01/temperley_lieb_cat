use num::Num;

use crate::gcd::GCD;
use crate::num::{One, Zero};
use crate::tex::Tex;

const MAX_DEGREE : usize = 256;

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Polynomial<T : Num + GCD + Copy> {
    coeffs : [T; MAX_DEGREE],
    degree : usize,
}

impl<T> Polynomial<T>
where T : Num + GCD + Copy {
    fn new<I,S>(mut in_coeffs: I) -> Polynomial<T>
    where I : Iterator<Item=S>,
          S : std::ops::Deref<Target=T> {
        let mut coeffs = [T::zero() ; MAX_DEGREE];
        let mut degree = 0;
        for i in 0..MAX_DEGREE {
            if let Some(x) = in_coeffs.next() {
                coeffs[i] = *x;
                if !coeffs[i].is_zero() {
                    degree = i;
                }
            } else {
                break;
            }
       }
       assert!(in_coeffs.next().is_none(), "Overflow in polynomial construction");
       Polynomial {
           coeffs,
           degree,
       }
    }

    pub fn gen() -> Self {
        Polynomial::new(vec![T::zero(), T::one()].iter())
    }

    fn degree(&self) -> usize {
        self.degree
    }

    fn constant_term(&self) -> T {
        self.coeffs[0]
    }

    fn hightest_term(&self) -> T {
        self.coeffs[self.degree()]
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

impl<T> std::ops::Add for Polynomial<T>
where T : Num + GCD + Copy {
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

impl<T> std::ops::Neg for Polynomial<T>
where T : Num + GCD + Copy + std::ops::Neg<Output=T> {
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

impl<T> std::ops::Sub for Polynomial<T>
where T : Num + GCD + Copy + std::ops::Neg<Output=T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut ans = self.coeffs;
        for i in 0..other.degree()+1 {
            ans[i] = ans[i] - other.coeffs[i];
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

impl<T> std::ops::Mul for Polynomial<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        assert!(self.degree() + other.degree() < MAX_DEGREE, "Attempt to multiply polynomials with overflow");
        let mut ans = [T::zero() ; MAX_DEGREE];
        for i in 0..(self.degree() + other.degree() + 1) {
            for j in 0..i.min(self.degree())+1 {
                if i-j <= other.degree() {
                    ans[i] = ans[i] + self.coeffs[j] * other.coeffs[i-j];
                }
            }
        }
        Polynomial::new(ans.iter())
    }
}

impl<T> std::ops::Mul<T> for Polynomial<T>
where T : Num + GCD + Copy {
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

impl<T> std::ops::Div<T> for Polynomial<T>
where T : Num + GCD + Copy {
    type Output = Self;

    fn div(self, other: T) -> Self {
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

impl<T> std::ops::Div for Polynomial<T>
where T : Num + GCD + Copy + std::fmt::Display + std::ops::Neg<Output=T> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        assert!(!other.is_zero(), "Dividing polynomial by zero");
        if self == Polynomial::zero() {
            return self;
        }
        if self == -other {
            return -Polynomial::one();
        }
        assert!(self.degree() >= other.degree(), format!("Cannot divide polynomials {} {}", self, other));
        if other == Polynomial::one() {
            return self;
        }
        let mut rem = self;
        let mut ans = [T::zero() ; MAX_DEGREE];
        for i in (other.degree()..self.degree()+1).rev() {
            if rem.coeffs[i].is_zero() { continue }
            let m = rem.hightest_term() / other.hightest_term();
            let mut diff = other * m;
            diff.shift(i-other.degree());
            rem = rem - diff;
            ans[i - other.degree] = m;
            assert!(rem.degree() < i, format!("Cannot divide polynomial {} by {}", self,other));
        }
        assert!(rem.is_zero());
        Polynomial::new(ans.iter())
    }
}

impl<T> num::One for Polynomial<T>
where T : Num + GCD + Copy + std::fmt::Display {
    fn one() -> Self {
        Polynomial::new(vec![T::one()].iter())
    }
}

impl<T> num::Zero for Polynomial<T>
where T : Num + GCD + Copy {
    fn zero() -> Self {
        Polynomial::new(vec![T::zero()].iter())
    }

    fn is_zero(&self) -> bool {
        self.degree() == 0 && self.hightest_term().is_zero()
    }
}

impl<T> std::ops::Rem for Polynomial<T>
where T : Num + GCD + Copy + std::fmt::Display + std::ops::Neg<Output=T> {
    type Output = Self;

    fn rem(self, other:Self) -> Self {
        self - self / other
    }
}


impl<T> Num for Polynomial<T>
where T : Num + GCD + Copy + std::fmt::Display + std::ops::Neg<Output=T> {
    type FromStrRadixErr = ();

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl<T> GCD for Polynomial<T>
where T:Num + GCD + Copy + std::ops::Neg<Output=T> + std::fmt::Display {
    fn gcd(&self, other: &Polynomial<T>) -> Polynomial<T> {
        if self.degree() < other.degree() {
            return other.gcd(self);
        }
        if other.is_zero() {
            if self.is_zero() {
                return *self;
            } else {
                let g = self.coeffs[..self.degree()+1].iter().fold(self.hightest_term(), |a, b| a.gcd(b));
                return *self / g;
            }
        }
        if other.degree() == 0 {
            if other.coeffs[0].is_one() || (-other.coeffs[0]).is_one() {
                return Polynomial::one();
            }
        }
        let mut diff = *other;
        diff.shift(self.degree() - other.degree());
        diff = diff * self.hightest_term() - *self * diff.hightest_term();
        let g = diff.coeffs[..diff.degree()+1].iter().fold(diff.hightest_term(), |a,b| a.gcd(b));
        if !g.is_zero() && !g.is_one() {
            diff = diff / g;
        }
        other.gcd(&diff)
    }
}

impl<T> std::fmt::Display for Polynomial<T>
where T: Num + GCD + Copy + std::fmt::Display {
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

impl<T> std::fmt::Debug for Polynomial<T>
where T: Num + GCD + Copy + std::fmt::Display {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl<T> Tex for Polynomial<T>
where T: Num + GCD + Copy + std::fmt::Display {
    fn into_tex(&self) -> String {
        format!("{}", self)
    }
}

pub fn quantum(n : i128) -> Polynomial<i128> {
    if n == 0 {
        return Polynomial::zero();
    }
    fn fact(n : i128) -> i128 {
        (1..n+1).fold(1, |acc, e| acc * e)
    }
    fn choose(n : i128, r: i128) -> i128 {
        fact(n) / (fact(r) * fact(n - r))
    }
    let mut coeffs = Vec::new();
    for _ in 0..n {
        coeffs.push(0)
    }
    for i in 0..(n-1)/2+1 {
        coeffs[(n-1-2*i) as usize] = (1-2*(i%2)) * choose(n-1-i,i);
    }
    Polynomial::new(coeffs.iter())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new () {
        debug_assert_eq!(Polynomial::new(vec![1,2,3,4,0,0,0].iter()).degree(), 3);
        debug_assert_eq!(Polynomial::new(vec![1,2,3,4,0,0,0].iter()) , Polynomial::new(vec![1,2,3,4].iter()));
        debug_assert_eq!(Polynomial::new(vec![0,0,0].iter()).degree() , 0);
    }

    #[test]
    fn gen() {
        debug_assert_eq!(Polynomial::<i128>::gen(), Polynomial::new(vec![0,1].iter()));
        debug_assert_eq!(Polynomial::<i128>::gen().degree(), 1);
        debug_assert_eq!("1X", format!("{}", Polynomial::<i128>::gen()));
        debug_assert_ne!(Polynomial::<i128>::gen(), Polynomial::zero());
        debug_assert_ne!(Polynomial::<i128>::gen(), Polynomial::one());
    }

    #[test]
    fn additive() {
        let a = Polynomial::new(vec![1,2,3].iter());
        let b = Polynomial::new(vec![1,2,3].iter());
        debug_assert_eq!(Polynomial::zero(), a - b);
    }

    #[test]
    fn neg() {
        debug_assert_eq!(-Polynomial::new(vec![1,2,3].iter()), Polynomial::new(vec![-1,-2,-3].iter()));
    }

    #[test]
    fn sub() {
        debug_assert_eq!(Polynomial::new(vec![1,2,3].iter()) - Polynomial::new(vec![2,0,1,5].iter()),
        Polynomial::new(vec![-1,2,2,-5].iter()));
    }

    #[test]
    fn mul() {
        let a = Polynomial::new(vec![1,2,3].iter());
        let b = Polynomial::new(vec![0,4,6].iter());
        debug_assert_eq!(Polynomial::new(vec![0,4,14,24,18].iter()), a * b);
    }

    #[test]
    fn div() {
        let a = Polynomial::new(vec![16384 , - 20480 , - 45056 , 61440 , 45056 , - 70400 , - 19712 , 38720 , 3520 , - 10560 , - 176 , 1320 ,0,  -55 , 0,0].iter().rev());
        let b = Polynomial::new(vec![0,5,0,-20,16].iter());
        debug_assert_eq!(Polynomial::new(
                vec![1024, 0, -2816, 0, 2816, 0, -1232, 0, 220, 0, -11, 0].iter().rev())
                , a / b);
    }

    #[test]
    fn quantum_num() {
        debug_assert_eq!(quantum(0), Polynomial::zero());
        debug_assert_eq!(quantum(1), Polynomial::one());
        debug_assert_eq!(quantum(2), Polynomial::gen());
        debug_assert_eq!(quantum(3), Polynomial::new(vec![-1,0,1].iter()));
        debug_assert_eq!(quantum(8), Polynomial::new(vec![0,-4,0,10,0,-6,0,1].iter()));
    }

    #[test]
    fn gcd() {
        fn mon(x : Polynomial<i128>) -> Polynomial<i128> {
            if x.hightest_term() < 0 {
                - x
            } else {
                x
            }
        }
        debug_assert_eq!(mon(Polynomial::new(vec![1,2,1].iter()).gcd(&Polynomial::new(vec![1,1].iter()))) , Polynomial::new(vec![1,1].iter()));
        debug_assert_eq!(mon(Polynomial::new(vec![2,4,2].iter()).gcd(&Polynomial::new(vec![4,4].iter()))) , Polynomial::new(vec![1,1].iter()));
        debug_assert_eq!(mon(quantum(17).gcd(&quantum(15))) , Polynomial::one());
        debug_assert_eq!(mon(quantum(17).gcd(&quantum(15))) , Polynomial::one());
        debug_assert_eq!(mon(quantum(27).gcd(&quantum(3))) , quantum(3));
        debug_assert_eq!(mon(quantum(27).gcd(&quantum(6))) , quantum(3));
    }
}
