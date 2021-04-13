use num_traits::ToPrimitive;
use std::fmt::{Debug, Display};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::fraction::Fraction;
use crate::gcd::PartialGCD;
use crate::serial::Serialisable;
use crate::structures::{Domain, GCDDomain, Ring, RingOps, Signed, Q};
use crate::tex::Tex;

/// A polynomial in a single variable over a _ring_.
///
/// It is implemented as a vector of coefficients, where the
/// highest coefficient is always zero.  The length is at least
/// two.
#[derive(Clone, PartialEq, Eq)]
pub struct Polynomial<T>(Vec<T>);

impl<R: Ring> Polynomial<R> {
    pub fn new<V>(in_coeffs: &[V]) -> Polynomial<R>
    where
        V: Clone + Into<R>,
    {
        let mut coeffs: Vec<R> = in_coeffs.iter().map(|x| x.clone().into()).collect();
        while let Some(x) = coeffs.last() {
            if x.is_zero() {
                coeffs.pop();
            } else {
                break;
            }
        }
        if coeffs.is_empty() {
            coeffs.push(R::zero());
        }
        coeffs.push(R::zero());
        Polynomial(coeffs)
    }

    pub fn degree(&self) -> usize {
        self.0.len() - 2
    }

    pub fn hightest_term(&self) -> &R {
        &self.0[self.0.len() - 2]
    }

    pub fn coeff(&self, i: usize) -> &R {
        if i > self.degree() {
            self.0.last().unwrap()
        } else {
            self.0.get(i).unwrap()
        }
    }

    fn shift(&mut self, n: usize) {
        let mut coeffs = vec![R::zero(); n];
        coeffs.extend(self.0.clone().into_iter());
        self.0 = coeffs;
    }

    pub fn eval<S>(&self, x: S, one: S) -> S
    where
        for<'r> S: Clone + Add<Output = S> + Mul<Output = S> + Mul<&'r R, Output = S>,
    {
        let mut ans = one * self.coeff(0);
        let mut xp = x.clone();
        for i in 1..self.degree() + 1 {
            ans = ans + xp.clone() * self.coeff(i);
            xp = xp * x.clone();
        }
        ans
    }

    pub fn gen() -> Self {
        Polynomial::new(&[R::zero(), R::one()])
    }
}
impl<R: Ring + Signed> Polynomial<R>
where
    for<'r> &'r R: RingOps<&'r R, R> + Div<&'r R, Output = R>,
{
    fn div_mod(&self, other: &Self) -> (Self, Self) {
        assert!(!other.is_zero(), "Dividing polynomial by zero");
        if self.is_zero() {
            return (Polynomial::zero(), Polynomial::zero());
        }
        if *self == *other {
            return (Polynomial::one(), Polynomial::zero());
        }
        if *self == -other {
            return (-Polynomial::one(), Polynomial::zero());
        }
        if self.degree() < other.degree() {
            return (Polynomial::zero(), self.clone());
        }
        let mut rem = self.clone();
        let mut ans = vec![R::zero(); self.degree() + 1];
        for i in (other.degree()..self.degree() + 1).rev() {
            //TODO(robert) is this tautological?
            if i >= rem.degree() + 1 || rem.coeff(i).is_zero() {
                continue;
            }
            // Because T is guaranteed to be a field, this gives no remainder
            let m = rem.hightest_term() / other.hightest_term();
            let mut diff = other.clone() * m.clone();
            diff.shift(i - other.degree());
            rem = rem - diff;
            ans[i - other.degree()] = m;
            assert!(rem.degree() < i || i == 0);
        }
        (Polynomial::new(&ans), rem)
    }
}

impl<T: Ring> Add<Polynomial<T>> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn add(self, other: Polynomial<T>) -> Polynomial<T> {
        &self + &other
    }
}

impl<T: Ring> Add<&Polynomial<T>> for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn add(self, other: &Polynomial<T>) -> Polynomial<T> {
        let deg = self.degree().max(other.degree());
        let mut ans = vec![T::zero(); deg + 1];
        for i in 0..deg + 1 {
            ans[i] = self.coeff(i) + other.coeff(i)
        }
        Polynomial::new(&ans)
    }
}

impl<T> Neg for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn neg(self) -> Polynomial<T> {
        -&self
    }
}

impl<T> Neg for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn neg(self) -> Polynomial<T> {
        Polynomial(self.0.iter().map(|x| -x).collect())
    }
}

impl<T: Ring> Sub<Polynomial<T>> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn sub(self, other: Polynomial<T>) -> Polynomial<T> {
        &self - &other
    }
}

impl<T: Ring> Sub<&Polynomial<T>> for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn sub(self, other: &Polynomial<T>) -> Polynomial<T> {
        let deg = self.degree().max(other.degree());
        let mut ans = vec![T::zero(); deg + 1];
        for i in 0..self.degree() + 1 {
            ans[i] = self.coeff(i).clone();
        }
        for i in 0..other.degree() + 1 {
            ans[i] = &ans[i] - other.coeff(i)
        }
        Polynomial::new(&ans)
    }
}

impl<T: Ring> Sub<T> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn sub(self, other: T) -> Polynomial<T> {
        self - Polynomial::from(other)
    }
}

impl<T: Ring> Mul<Polynomial<T>> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Polynomial<T> {
        &self * &other
    }
}

impl<T: Ring> Mul<&Polynomial<T>> for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn mul(self, other: &Polynomial<T>) -> Polynomial<T> {
        let max_resulting_degree = self.degree() + other.degree();
        let mut ans = vec![T::zero(); max_resulting_degree + 1];
        for i in 0..max_resulting_degree + 2 {
            for j in i.saturating_sub(other.degree())..i.min(self.degree()) + 1 {
                // In the land of Mordor, where the shadows lie...
                ans[i] = &ans[i] + &(self.coeff(j) * other.coeff(i - j));
            }
        }
        Polynomial::new(&ans)
    }
}

impl<T: Ring> Mul<T> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Polynomial::new(&self.0.iter().map(|x| x * &other).collect::<Vec<_>>())
    }
}

impl<T: Ring> Div<T> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Self;

    fn div(self, other: T) -> Self {
        assert!(!other.is_zero(), "Cannot divide polynomial by zero");
        Polynomial::new(&self.0.iter().map(|x| x / &other).collect::<Vec<_>>())
    }
}

impl<T: Ring + Signed> Div<Polynomial<T>> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn div(self, other: Polynomial<T>) -> Polynomial<T> {
        &self / &other
    }
}

impl<T: Ring + Signed> Div<&Polynomial<T>> for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn div(self, other: &Polynomial<T>) -> Polynomial<T> {
        self.div_mod(other).0
    }
}

impl<T: Ring + Signed> Rem<Polynomial<T>> for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn rem(self, other: Polynomial<T>) -> Polynomial<T> {
        &self % &other
    }
}

impl<T: Ring + Signed> Rem<&Polynomial<T>> for &Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    type Output = Polynomial<T>;

    fn rem(self, other: &Polynomial<T>) -> Polynomial<T> {
        self.div_mod(&other).1
    }
}

impl<T: Ring + Signed> Signed for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
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
}

/// Casting an element to a polynomial in that element is via
/// the usual injection of R into R[X].
impl<T: Ring> From<T> for Polynomial<T> {
    fn from(val: T) -> Polynomial<T> {
        Polynomial::new(&[val])
    }
}

impl PartialGCD for Polynomial<Q> {
    fn partial_gcd(&self, other: &Polynomial<Q>) -> Polynomial<Q> {
        if (self.0.iter().fold(0, |a, x| x.den().abs().max(a)) > 0xffffffff)
            || (self.0.iter().fold(0, |a, x| x.num().abs().max(a)) > 0xffffffff)
            || (other.0.iter().fold(0, |a, x| x.den().abs().max(a)) > 0xfffffffff)
            || (other.0.iter().fold(0, |a, x| x.num().abs().max(a)) > 0xfffffffff)
        {
            return Polynomial::one();
        }
        if self.degree() < other.degree() {
            return other.partial_gcd(self);
        }
        if other.is_zero() {
            return self.clone();
        }
        let div_rem = self.div_mod(other);
        if div_rem.0.is_zero() {
            Polynomial::one()
        } else {
            let mut re = div_rem.1;
            //let c = re.coeffs.iter().fold(re.hightest_term().den(), |x,y | x * y.den() / x.gcd(&y.den()));
            let c = re.hightest_term();
            if !c.is_zero() {
                re = re.clone() / c.clone();
            }
            let sub = other.partial_gcd(&re);
            sub.clone() / sub.clone().hightest_term().clone()
        }
    }

    fn is_small(&self) -> bool {
        self.degree() == 0
    }
}

impl<T: Ring + Display> Display for Polynomial<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut s = false;
        for i in 0..self.degree() + 1 {
            if !self.coeff(i).is_zero() || (i == 0 && self.degree() == 0) {
                if s {
                    write!(f, " + ")?;
                }
                s = true;
                if !self.coeff(i).is_one() || i == 0 {
                    write!(f, "{}", self.coeff(i))?;
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

impl<T: Ring + Display> Debug for Polynomial<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl<T: Ring + Signed + Tex> Tex for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn into_tex(&self) -> String {
        let mut ans = String::new();
        let mut s = false;
        for i in 0..self.degree().max(0) + 1 {
            if !self.coeff(i).is_zero() || (i == 0 && self.degree() == 0) {
                if s {
                    if self.coeff(i).is_positive() {
                        ans += " + ";
                    }
                }
                s = true;
                if !self.coeff(i).abs().is_one() || i == 0 {
                    if self.coeff(i).is_multiterm() {
                        ans += "\\left(";
                    }
                    ans += &self.coeff(i).into_tex();
                    if self.coeff(i).is_multiterm() {
                        ans += "\\right)";
                    }
                }
                if self.coeff(i).is_negative() && (-self.coeff(i)).is_one() && i > 0 {
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
        self.0.iter().filter(|x| !x.is_zero()).count() > 1
    }
}

pub fn quantum<I: ToPrimitive>(n: I) -> Polynomial<Q> {
    let n: i128 = n.to_i64().unwrap() as i128;
    if n.is_zero() {
        return Polynomial::zero();
    }
    if n < 0 {
        return -quantum(-n);
    }
    fn choose(n: i128, r: i128) -> i128 {
        let mut ans = 1;
        for i in n - r + 1..n + 1 {
            ans = ans * i;
        }
        for i in 1..r + 1 {
            ans = ans / i;
        }
        ans
    }
    let mut coeffs = Vec::new();
    for _ in 0..n {
        coeffs.push(Fraction::zero())
    }
    for i in 0..(n - 1) / 2 + 1 {
        coeffs[(n - 1 - 2 * i) as usize] = Fraction::from((1 - 2 * (i % 2)) * choose(n - 1 - i, i));
    }
    Polynomial::new(coeffs.as_slice())
}

impl Serialisable for Polynomial<Q> {
    fn serialise(&self) -> String {
        self.0.serialise()
    }

    fn deserialise(inpt: &str) -> Self {
        let res: Vec<Q> = Vec::deserialise(inpt);
        Polynomial::new(res.as_slice())
    }
}

impl<T: Ring + Signed> RingOps<Polynomial<T>, Polynomial<T>> for Polynomial<T> where
    for<'r> &'r T: RingOps<&'r T, T>
{
}

impl<T: Ring + Signed> Ring for Polynomial<T>
where
    for<'r> &'r T: RingOps<&'r T, T>,
{
    fn zero() -> Self {
        let c: [T; 0] = [];
        Polynomial::new(&c)
    }

    fn is_zero(&self) -> bool {
        self.degree() == 0 && self.hightest_term().is_zero()
    }
    fn one() -> Self {
        Polynomial::new(&[T::one()])
    }
    fn is_one(&self) -> bool {
        self.degree() == 0 && self.hightest_term().is_one()
    }
}

impl<T: Domain + Signed> Domain for Polynomial<T> where for<'r> &'r T: RingOps<&'r T, T> {}

impl GCDDomain for Polynomial<Q> {}

impl<'a, T: 'a + Ring + Signed> RingOps<&'a Polynomial<T>, Polynomial<T>> for &'a Polynomial<T> where
    for<'r> &'r T: RingOps<&'r T, T>
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        assert_eq!(Polynomial::<Q>::new(&[1, 2, 3, 4, 0, 0, 0]).degree(), 3);
        assert_eq!(
            Polynomial::<Q>::new(&[1, 2, 3, 4, 0, 0, 0]),
            Polynomial::new(&[1, 2, 3, 4])
        );
        assert_eq!(Polynomial::<Q>::new(&[0, 0, 0]).degree(), 0);
    }

    #[test]
    fn gen() {
        assert_eq!(
            Polynomial::<Q>::gen(),
            Polynomial::new(&[Q::zero(), Q::one()])
        );
        assert_eq!(Polynomial::<Q>::gen().degree(), 1);
        assert_eq!("X", format!("{}", Polynomial::<Q>::gen()));
        assert_ne!(Polynomial::<Q>::gen(), Polynomial::zero());
        assert_ne!(Polynomial::<Q>::gen(), Polynomial::one());
    }

    #[test]
    fn additive() {
        let a = Polynomial::<Q>::new(&[1, 2, 3]);
        let b = Polynomial::<Q>::new(&[1, 2, 3]);
        assert_eq!(Polynomial::zero(), a - b);
    }

    #[test]
    fn neg() {
        assert_eq!(
            -Polynomial::<Q>::new(&[1, 2, 3]),
            Polynomial::new(&[-1, -2, -3])
        );
    }

    #[test]
    fn sub() {
        assert_eq!(
            Polynomial::<Q>::new(&[1, 2, 3]) - Polynomial::new(&[2, 0, 1, 5]),
            Polynomial::<Q>::new(&[-1, 2, 2, -5])
        );
    }

    #[test]
    fn mul() {
        let a = Polynomial::<Q>::new(&[1, 2, 3]);
        let b = Polynomial::<Q>::new(&[0, 4, 6]);
        assert_eq!(Polynomial::new(&[0, 4, 14, 24, 18]), a * b);
    }

    #[test]
    fn div() {
        let a = Polynomial::<Q>::new(&[
            0, 0, -55, 0, 1320, -176, -10560, 3520, 38720, -19712, -70400, 45056, 61440, -45056,
            -20480, 16384,
        ]);
        let b = Polynomial::new(&[0, 5, 0, -20, 16]);
        assert_eq!(
            Polynomial::<Q>::new(&[0, -11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024]),
            a / b
        );
    }

    #[test]
    fn quantum_num() {
        assert_eq!(quantum(0), Polynomial::zero());
        assert_eq!(quantum(1), Polynomial::one());
        assert_eq!(quantum(2), Polynomial::gen());
        assert_eq!(quantum(3), Polynomial::new(&[-1, 0, 1]));
        assert_eq!(quantum(8), Polynomial::new(&[0, -4, 0, 10, 0, -6, 0, 1]));
    }

    #[test]
    fn gcd() {
        assert_eq!(
            Polynomial::new(&[1, 2, 1]).partial_gcd(&Polynomial::new(&[1, 1])),
            Polynomial::new(&[1, 1])
        );
        assert_eq!(
            Polynomial::new(&[2, 4, 2]).partial_gcd(&Polynomial::new(&[4, 4])),
            Polynomial::new(&[1, 1])
        );
        assert_eq!(quantum(17).partial_gcd(&quantum(15)), Polynomial::one());
        assert_eq!(quantum(17).partial_gcd(&quantum(15)), Polynomial::one());
        assert_eq!(quantum(27).partial_gcd(&quantum(3)), quantum(3));
        assert_eq!(quantum(27).partial_gcd(&quantum(6)), quantum(3));
    }

    #[test]
    fn big_gcd() {
        let a = Polynomial::<Q>::new(&[
            24, 0, -524, 0, 5056, 0, -27152, 0, 84534, 0, -138895, 0, 23826, 0, 445831, 0,
            -1177984, 0, 1742628, 0, -1768127, 0, 1310596, 0, -728166, 0, 305827, 0, -96767, 0,
            22722, 0, -3841, 0, 442, 0, -31, 0, 1,
        ]);

        let b = Polynomial::<Q>::new(&[
            Fraction::new(304, 5),
            Fraction::new(0, 1),
            Fraction::new(21856, -15),
            Fraction::new(0, 1),
            Fraction::new(238652, 15),
            Fraction::new(0, 1),
            Fraction::new(1527248, -15),
            Fraction::new(0, 1),
            Fraction::new(1258064, 3),
            Fraction::new(0, 1),
            Fraction::new(17683936, -15),
            Fraction::new(0, 1),
            Fraction::new(35414599, 15),
            Fraction::new(0, 1),
            Fraction::new(-52134644, 15),
            Fraction::new(0, 1),
            Fraction::new(57703364, 15),
            Fraction::new(0, 1),
            Fraction::new(48757748, -15),
            Fraction::new(0, 1),
            Fraction::new(10579261, 5),
            Fraction::new(0, 1),
            Fraction::new(-1064413, 1),
            Fraction::new(0, 1),
            Fraction::new(6190517, 15),
            Fraction::new(0, 1),
            Fraction::new(-366422, 3),
            Fraction::new(0, 1),
            Fraction::new(406151, 15),
            Fraction::new(0, 1),
            Fraction::new(65266, -15),
            Fraction::new(0, 1),
            Fraction::new(7177, 15),
            Fraction::new(0, 1),
            Fraction::new(161, -5),
            Fraction::new(0, 1),
            Fraction::new(1, 1),
        ]);
        print!("{}", a.partial_gcd(&b));
    }
}
