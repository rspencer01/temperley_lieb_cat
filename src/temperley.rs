use std::collections::HashMap;
use crate::temperley_diagram::TLDiagram;
use crate::poly::{quantum, Polynomial};
use crate::fraction::Fraction;
use crate::num::{Num, Zero};

#[derive(Clone, Debug)]
pub struct TLMorphism<R>
where R : Copy + Clone + Num + std::fmt::Display + std::fmt::Debug {
    coeffs : HashMap<TLDiagram,R>,
    delta : Option<R>,
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> TLMorphism<R> {
    fn new(coeffs : Vec<(TLDiagram, R)>, delta: Option<R>) -> TLMorphism<R> {
        let domain = coeffs.first().unwrap().0.domain();
        let co_domain = coeffs.first().unwrap().0.co_domain();
        let mut ans = HashMap::new();
        for (d,v) in coeffs.iter() {
            assert!(d.domain() == domain);
            assert!(d.co_domain() == co_domain);
            if !v.is_zero() {
                if ans.contains_key(d) {
                    ans.insert(d.clone(), ans[d] + *v);
                } else {
                    ans.insert(d.clone(), *v);
                }
            }
        }
        if ans.is_empty() {
            ans.insert(coeffs[0].0.clone(), coeffs[0].1);
        }
        TLMorphism{
            coeffs:ans,
            delta
        }
    }

    fn ring_compatible(&self, other : &TLMorphism<R>) -> bool {
        other.delta.is_none() || self.delta.is_none() || self.delta == other.delta
    }

    fn ring_point(&self, other : &TLMorphism<R>) -> Option<R> {
        assert!(self.ring_compatible(other));
        match self.delta {
            None => other.delta,
            Some(x) => Some(x),
        }
    }

    fn repoint(&mut self, delta : R) {
        self.delta = Some(delta);
    }

    fn domain(&self) -> usize {
        self.coeffs.keys().next().unwrap().domain()
    }

    fn co_domain(&self) -> usize {
        self.coeffs.keys().next().unwrap().co_domain()
    }

    fn involute(&self) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter().map(|(k, v)| (k.clone().involute(), *v)).collect(),
            self.delta
        )
    }

    fn id(n : usize) -> TLMorphism<R> {
        TLDiagram::id(n).into()
    }

    fn zero(n : usize) -> TLMorphism<R> {
        TLMorphism::id(n) * R::zero()
    }

    fn u(n : usize, i : usize) -> TLMorphism<R> {
        TLDiagram::u(n,i).into()
    }

    fn inject(&self) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter()
            .map(|(k,v)| (k.inject(), *v))
            .collect(),
            self.delta
        )
    }

    pub fn is_jones_wenzl(&self) -> bool {
        self.domain() == self.co_domain() &&
            (self.coeffs.get(&TLDiagram::id(self.domain())).unwrap_or(&R::zero()).is_one()) &&
            (1..self.domain()).all(|i|
            (TLMorphism::<R>::u(self.domain(),i) * self.clone()).is_zero())
    }

    pub fn is_idempotent(&self) -> bool {
        self.clone() * self.clone() == *self
    }

    pub fn support(&self) -> Vec<TLDiagram> {
        self.coeffs.iter()
            .filter(|(_,v)| !v.is_zero())
            .map(|(k,_)| k.clone())
            .collect()
    }
}

pub fn jw(n : usize) -> TLMorphism<Fraction<Polynomial<i128>>> {
    let mut jw = TLMorphism::id(1);
    jw.repoint(Polynomial::gen().into());
    for i in 1..n {
        let jwp = jw.inject();
        jw = jwp.clone() - jwp.clone() * TLMorphism::u(i+1,i) * jwp.clone() * (Fraction::from(quantum(i as i128)) / quantum(i as i128+1));
    }
    jw
}

pub fn jw_nat(n : usize) -> TLMorphism<Fraction<i128>> {
    let mut jw = TLMorphism::id(1);
    jw.repoint(2.into());
    for i in 1..n {
        let jwp = jw.inject();
        jw = jwp.clone() - jwp.clone() * TLMorphism::u(i+1,i) * jwp.clone() * (Fraction::from(i as i128) / (i as i128+1));
    }
    jw
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug > PartialEq for TLMorphism<R> {
    fn eq(&self, other : &TLMorphism<R>) -> bool {
        if self.co_domain() != other.co_domain() {
            return false;
        }
        if self.domain() != other.domain() {
            return false;
        }
        for (k,v) in self.coeffs.iter() {
            if other.coeffs.get(k) != Some(v)  && !v.is_zero() {
                return false;
            }
        }
        for (k,v) in other.coeffs.iter() {
            if self.coeffs.get(k) != Some(v)  && !v.is_zero() {
                return false;
            }
        }
        self.ring_compatible(other)
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Add for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn add(self, other: TLMorphism<R>) -> TLMorphism<R> {
        let mut ans = Vec::new();
        for (k, v) in self.coeffs.iter() {
            ans.push((k.clone(), match other.coeffs.get(k) {
                Some(vp) => {*v + *vp},
                None => {*v},
            }));
        }
        for (k, v) in other.coeffs.iter() {
            match self.coeffs.get(k) {
                Some(_) => {},
                None => {ans.push((k.clone(), *v))},
            };
        }
        TLMorphism::new(ans, self.ring_point(&other))
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Sub for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn sub(self, other: TLMorphism<R>) -> TLMorphism<R> {
        self + (-other)
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Neg for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn neg(self) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter()
            .map(|(k,v)| (k.clone(), R::zero()-*v))
            .collect(),
            self.delta
            )
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Mul for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn mul(self, other: TLMorphism<R>) -> TLMorphism<R> {
        let mut ans = TLMorphism::new(vec![(TLDiagram::any(self.domain(), self.co_domain()), R::zero())], self.delta);
        let delta = self.ring_point(&other).expect("Require ring point to multiply morphisms");
        let mut pows = vec![R::one(); (self.domain() + self.co_domain()) / 2];
        for i in 1.. pows.len() {
            pows[i] = pows[i-1]*delta;
        }
        for (d, v) in self.coeffs.iter() {
            ans = ans + TLMorphism::new(
                other.coeffs.iter().map(|(dp, vp)|{
                    let m = d.clone() * dp.clone();
                    (m.1, pows[m.0] * *v * *vp)
                }).collect(),
                Some(delta)
            )
        }
        ans
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::BitOr for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn bitor(self, other: TLMorphism<R>) -> TLMorphism<R> {
        let mut ans = TLMorphism::new(vec![
            (TLDiagram::any(self.domain() + other.domain(), self.co_domain()+other.co_domain()), R::zero())
        ], self.delta);
        let delta = self.ring_point(&other).expect("Require ring point to multiply morphisms");
        for (d, v) in self.coeffs.iter() {
            ans = ans + TLMorphism::new(
                other.coeffs.iter().map(|(dp, vp)|{
                    (d.clone() | dp.clone(), *v * *vp)
                }).collect(),
                Some(delta)
            )
        }
        ans
    }
}


impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Mul<R> for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn mul(self, other : R) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter()
            .map(|(k, v)| (k.clone(), *v * other))
            .collect(),
            self.delta
        )
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> std::ops::Div<R> for TLMorphism<R> {
    type Output = TLMorphism<R>;

    fn div(self, other : R) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter()
            .map(|(k, v)| (k.clone(), *v / other))
            .collect(),
            self.delta
        )
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> From<TLDiagram> for TLMorphism<R> {
    fn from(diag : TLDiagram) -> TLMorphism<R> {
        TLMorphism::new(
            vec![(diag, R::one())],
            None
        )
    }
}

impl<R:Copy + Clone + Num + std::fmt::Display + std::fmt::Debug> num::Zero for TLMorphism<R> {
    fn zero() -> TLMorphism<R> {
        unimplemented!()
    }

    fn is_zero(&self) -> bool {
        self.coeffs.values().all(|x| x.is_zero())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{quantum, Polynomial};
    use crate::fraction::Fraction;
    use crate::num::One;
    use crate::temperley_link::Link;
    use crate::temperley_site::Site::*;

    #[test]
    fn equality() {
        debug_assert_eq!(TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::from_tableauxs(3, vec![2], vec![2]), -1),
        ], None),
        TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(vec![Link::new(Source(1), Source(2)), Link::new(Target(1), Target(2)), Link::new(Target(3), Source(3))]), -1),
        ], None));
        debug_assert_eq!(TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::from_tableauxs(3,vec![2], vec![2]), -1),
        ], None),
        TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::from_tableauxs(3,vec![2], vec![2]), -1),
            (TLDiagram::from_tableauxs(3,vec![3], vec![2]), 0),
        ], None));
    }

    #[test]
    fn add() {
        debug_assert_eq!(TLMorphism::<i128>::zero(5) + TLMorphism::zero(5), TLMorphism::zero(5));
        debug_assert_ne!(TLMorphism::<i128>::zero(5) + TLMorphism::zero(5), TLMorphism::zero(6));
        debug_assert_eq!(TLMorphism::<i128>::zero(5) + TLMorphism::id(5), TLMorphism::id(5));
        debug_assert_eq!(TLMorphism::<i128>::id(5) + TLMorphism::zero(5), TLMorphism::id(5));
    }

    #[test]
    fn mul() {
        debug_assert_eq!(
            TLMorphism::zero(6),
            TLMorphism::id(6) * (0 as i128)
        );
        type R = Fraction<Polynomial<i128>>;
        let delta = Fraction::from(Polynomial::gen());
        let mut a : TLMorphism<R> = TLMorphism::from(TLDiagram::u(6, 3));
        a.repoint(delta);
        debug_assert_eq!(
            a.clone() * a.clone(),
            a.clone() * delta,
        );
    }

    #[test]
    fn manual_jw3() {
        type R = Fraction<Polynomial<i128>>;
        let delta : R = Polynomial::gen().into();
        let jw3 = TLMorphism::new(vec![
            (TLDiagram::id(3), R::one()),
            (TLDiagram::from_tableauxs(3,vec![2], vec![2]), -R::from(quantum(2)) / quantum(3)),
            (TLDiagram::from_tableauxs(3,vec![3], vec![3]), -R::from(quantum(2)) / quantum(3)),
            (TLDiagram::from_tableauxs(3,vec![2], vec![3]), R::from(quantum(1)) / quantum(3)),
            (TLDiagram::from_tableauxs(3,vec![3], vec![2]), R::from(quantum(1)) / quantum(3)),
        ], Some(Polynomial::gen().into()));
        debug_assert_eq!(jw3, jw3.involute());
        let mut a = TLMorphism::from(TLDiagram::u(3,1));
        let mut b = TLMorphism::from(TLDiagram::u(3,1));
        a.repoint(delta);
        debug_assert_eq!(a.clone() * jw3.clone(), TLMorphism::zero(3));
        debug_assert_eq!(b.clone() * jw3.clone(), TLMorphism::zero(3));
        debug_assert_eq!(jw3.clone() * a.clone(), TLMorphism::zero(3));
        debug_assert_eq!(jw3.clone() * b.clone(), TLMorphism::zero(3));
        debug_assert_eq!(jw3.clone() * jw3.clone(), jw3);
        assert!(jw3.is_jones_wenzl());
    }

    #[test]
    fn general_jw() {
        type R = Fraction<Polynomial<i128>>;
        let mut jw1 = TLMorphism::<R>::id(1);
        jw1.repoint(Polynomial::gen().into());
        assert!(jw(1).is_jones_wenzl());
        assert!(jw(2).is_jones_wenzl());
        assert!(jw(3).is_jones_wenzl());
        assert!(jw(4).is_jones_wenzl());
        assert!(jw(5).is_jones_wenzl());
        assert!(jw(4).inject() * jw(5) == jw(5));
    }

    #[test]
    fn tensored_jw() {
        assert!((jw(4) | jw(3)).is_idempotent());
        assert!((jw(3) | jw(2) | jw(3)).is_idempotent()); 
        assert!(!((jw(3) | jw(2) | jw(3)).is_jones_wenzl()));
        assert!((jw(5) | jw(1)).is_idempotent());
        assert!(!(jw(5) | jw(1)).is_jones_wenzl());
        assert!((jw(4) | jw(2)) * jw(6) == jw(6));
    }
}
