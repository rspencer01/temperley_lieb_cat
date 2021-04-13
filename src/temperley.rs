use std::ops::{Add, Sub, Mul, Div, Neg, BitOr};
use std::collections::HashMap;
use crate::TLDiagram;
use crate::tex::Tex;
use crate::serial::Serialisable;
use crate::structures::{Signed, RingOps, Ring};

#[derive(Clone, Debug)]
pub struct TLMorphism<R : Ring> {
    pub coeffs : HashMap<TLDiagram,R>,
    domain : usize,
    co_domain : usize,
    pub delta : Option<R>,
    pub right_kills : Vec<usize>,
}

impl<R: Ring> TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    pub fn new(coeffs : Vec<(TLDiagram, R)>, delta: Option<R>) -> TLMorphism<R> {
        let representative_diagram = &coeffs.first()
            .expect("Attempt to construct empty morphism")
            .0;
        let domain = representative_diagram.domain();
        let co_domain = representative_diagram.co_domain();
        let mut coeffs_map : HashMap<TLDiagram, R> = HashMap::new();
        for (d,v) in coeffs.iter() {
            assert!(d.domain() == domain);
            assert!(d.co_domain() == co_domain);
            if !v.is_zero() {
                if coeffs_map.contains_key(d) {
                    coeffs_map.insert(*d, &coeffs_map[d] + v);
                } else {
                    coeffs_map.insert(*d, v.clone());
                }
            }
        }
        if coeffs_map.is_empty() {
            coeffs_map.insert(representative_diagram.clone(), R::zero());
        }

        let mut right_kills = Vec::new();
        if let Some(delta) = delta.clone() {
            let mut pows = vec![R::one(); co_domain];
            for i in 1.. pows.len() {
                pows[i] = &pows[i-1] * &delta;
            }
            for i in 1..co_domain{
                if TLMorphism::new(
                    coeffs_map.iter().map(|(dp, vp)|{
                        let m = dp * &TLDiagram::u(co_domain, i);
                        (m.1, &pows[m.0] * vp)
                    }).collect(),
                    None
                ).is_zero() {
                    right_kills.push(i);
                }
            }
        }

        TLMorphism{
            coeffs:coeffs_map,
            delta,
            right_kills,
            domain,
            co_domain,
        }
    }

    /// Checks if the possibly pointed morphisms can correspond to the same pointed ring
    ///
    /// Since morphisms can be possibly pointed, this checks that if they are both pointed
    /// they do not point to different deltas.
    fn ring_compatible(&self, other : &TLMorphism<R>) -> bool {
        other.delta.is_none() || self.delta.is_none() || self.delta == other.delta
    }

    pub fn domain(&self) -> usize {
        self.domain
    }

    pub fn co_domain(&self) -> usize {
        self.co_domain
    }

    fn ring_point(&self, other : &TLMorphism<R>) -> Option<R> {
        assert!(self.ring_compatible(other));
        match self.delta.clone() {
            None => other.delta.clone(),
            Some(x) => Some(x),
        }
    }

    pub fn repoint(&mut self, delta : Option<R>) {
        self.delta = delta.clone();
        let co_domain = self.co_domain();
        let mut right_kills = Vec::new();
        if delta.is_some() {
            let mut pows = vec![R::one(); co_domain];
            for i in 1.. pows.len() {
                pows[i] = &pows[i-1]*delta.as_ref().unwrap();
            }
            for i in 1..co_domain{
                let cap = TLDiagram::u(co_domain, i);
                if TLMorphism::new(
                    self.coeffs.iter().map(|(dp, vp)|{
                        let m = dp * &cap;
                        (m.1, &pows[m.0] * vp)
                    }).collect(),
                    None
                ).is_zero() {
                    right_kills.push(i);
                }
            }
        }
        self.right_kills = right_kills;
    }


    pub fn involute(&self) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter().map(|(k, v)| (k.involute(), v.clone())).collect(),
            self.delta.clone()
        )
    }

    pub fn flip(&self) -> TLMorphism<R> {
        TLMorphism::new(
            self.coeffs.iter().map(|(k, v)| (k.flip(), v.clone())).collect(),
            self.delta.clone()
        )
    }

    pub fn id(n : usize) -> TLMorphism<R> {
        TLDiagram::id(n).into()
    }

    pub fn id_zero(n : usize) -> TLMorphism<R> {
        TLMorphism::id(n) * R::zero()
    }

    pub fn u(n : usize, i : usize) -> TLMorphism<R> {
        TLDiagram::u(n,i).into()
    }

    pub fn big_u(n : usize, i : usize, j : usize) -> TLMorphism<R> {
        TLDiagram::big_u(n,i, j).into()
    }

    /// Check if this is the Jones Wenzl idempotent
    ///
    /// The Jones-Wenzl idempotent can be characterised by
    ///   * Being a map n -> n
    ///   * Having unit coefficient for the identity diagram
    ///   * Being killed by all caps on the right
    /// The other features (idempotent, invariant under involution)
    /// follow from these.
    pub fn is_jones_wenzl(&self) -> bool {
        self.domain == self.co_domain &&
            (self.coeffs.get(&TLDiagram::id(self.domain)).unwrap_or(&R::zero()).is_one()) &&
            self.right_kills.len() == self.co_domain - 1
    }

    pub fn is_idempotent(&self) -> bool {
        self * self == *self
    }

    pub fn support(&self) -> Vec<&TLDiagram> {
        self.coeffs.iter()
            .filter(|(_,v)| !v.is_zero())
            .map(|(k,_)| k)
            .collect()
    }

    pub fn map_on_diagrams<F>(&self, f : F) -> TLMorphism<R>
    where F : Fn(&TLDiagram) -> TLDiagram {
        TLMorphism::new(
            self.coeffs.iter()
                .map(|(k,v)| (f(k), v.clone()))
                .collect(),
            self.delta.clone()
        )
    }

    pub fn map_on_coeffs<F>(&self, f : F) -> TLMorphism<R>
    where F : Fn(&R) -> R {
        TLMorphism::new(
            self.coeffs.iter()
                .map(|(k,v)| (*k, f(v)))
                .collect(),
            self.delta.clone()
        )
    }

    pub fn turn_down(&self, i : isize) -> TLMorphism<R> {
        self.map_on_diagrams(|d| d.turn_down(i))
    }

    pub fn turn_up(&self, i : isize) -> TLMorphism<R> {
        self.map_on_diagrams(|d| d.turn_up(i))
    }

    pub fn rotate(&self, i : isize) -> TLMorphism<R> {
        self.map_on_diagrams(|d| d.rotate(i))
    }

    pub fn trace(&self) -> TLMorphism<R> {
        let u  = TLMorphism::from(TLDiagram::cup(self.domain() + 1, self.domain()));
        u * self.turn_down(1)
    }
}

impl<R: Ring> PartialEq for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
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

impl<R: Ring> Eq for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {}

impl<R: Ring> Add<TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn add(self, other: TLMorphism<R>) -> TLMorphism<R> {
        &self + &other
    }
}

impl<R: Ring> Add<&TLMorphism<R>> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn add(self, other: &TLMorphism<R>) -> TLMorphism<R> {
        let delta = self.ring_point(&other);
        let mut ans = Vec::new();
        ans.extend(self.coeffs.iter().map(|(k,v)| (*k, v.clone())));
        ans.extend(other.coeffs.iter().map(|(k,v)| (*k, v.clone())));
        TLMorphism::new(ans, delta)
    }
}

impl<R: Ring> Sub<TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn sub(self, other: TLMorphism<R>) -> TLMorphism<R> {
        &self - &other
    }
}

impl<R: Ring> Sub<&TLMorphism<R>> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn sub(self, other: &TLMorphism<R>) -> TLMorphism<R> {
        let delta = self.ring_point(&other);
        let mut ans = Vec::new();
        ans.extend(self.coeffs.iter().map(|(k,v)| (*k, v.clone())));
        ans.extend(other.coeffs.iter().map(|(k,v)|(*k, -v.clone())));
        TLMorphism::new(ans, delta)
    }
}

impl<R: Ring> Neg for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn neg(self) -> TLMorphism<R> {
        self.map_on_coeffs(|v| -v.clone())
    }
}

impl<R: Ring> Mul<TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other: TLMorphism<R>) -> TLMorphism<R> {
        &self * &other
    }
}

impl<R: Ring> Mul<&TLMorphism<R>> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other: &TLMorphism<R>) -> TLMorphism<R> {
        let mut ans_vec = Vec::new();
        let delta = self.ring_point(&other).expect("Require ring point to multiply morphisms");
        let mut pows = vec![R::one(); (self.domain() + self.co_domain()) / 2 + 1];
        for i in 1.. pows.len() {
            pows[i] = &pows[i-1]*&delta;
        }
        for (d, v) in other.coeffs.iter() {
            let mut dont = false;
            for j in d.simple_links().iter() {
                if self.right_kills.contains(j) {
                    dont = true;
                    break;
                }
            }
            if dont {continue}
            ans_vec.extend(
                self.coeffs.iter().map(|(dp, vp)|{
                    let m = dp * d;
                    // In the land of Mordor, where the shadows lie
                    (m.1, &(&pows[m.0] * vp) * v)
                })
            );
        }
        if ans_vec.is_empty() {
            TLMorphism::new(
                vec![(TLDiagram::any(self.domain, other.co_domain()), R::zero())],
                Some(delta)
            )
        } else {
            TLMorphism::new(
                ans_vec,
                Some(delta),
            )
        }
    }
}

impl<R: Ring> BitOr<TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    /// Bitwise or, `a | b`, is the tensor product
    fn bitor(self, other: TLMorphism<R>) -> TLMorphism<R> {
        &self | &other
    }
}

impl<R: Ring> BitOr<&TLMorphism<R>> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    /// Bitwise or, `a | b`, is the tensor product
    fn bitor(self, other: &TLMorphism<R>) -> TLMorphism<R> {
        let mut ans_vec = Vec::new();
        let delta = self.ring_point(other);
        for (d, v) in self.coeffs.iter() {
            ans_vec.extend(
                other.coeffs.iter().map(|(dp, vp)|{
                    (*d | dp, v * vp)
                })
            );
        }
        TLMorphism::new(
            ans_vec,
            delta
        )
    }
}

impl<R: Ring> Mul<&R> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other : &R) -> TLMorphism<R> {
        self.map_on_coeffs(|x| x * other)
    }
}

impl<R: Ring> Mul<R> for &TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other : R) -> TLMorphism<R> {
        self * &other
    }
}


impl<R: Ring> Mul<&R> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other : &R) -> TLMorphism<R> {
        &self * other
    }
}

impl<R: Ring> Mul<R> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    type Output = TLMorphism<R>;

    fn mul(self, other : R) -> TLMorphism<R> {
        self * &other
    }
}

impl<R: Ring> Div<R> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> + Div<&'r R, Output=R> {
    type Output = TLMorphism<R>;

    fn div(self, other : R) -> TLMorphism<R> {
        self.map_on_coeffs(|x| x / &other)
    }
}

impl<R: Ring> Div<TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> + Div<&'r R, Output=R> {
    type Output = TLMorphism<R>;

    fn div(self, other : TLMorphism<R>) -> TLMorphism<R> {
        unimplemented!()
    }
}

impl<R: Ring> From<TLDiagram> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    fn from(diag : TLDiagram) -> TLMorphism<R> {
        TLMorphism::new(
            vec![(diag, R::one())],
            None
        )
    }
}

impl<R: Ring> RingOps<TLMorphism<R>, TLMorphism<R>> for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {}

impl<R: Ring> Ring for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    fn zero() -> TLMorphism<R> {
        unimplemented!()
    }

    fn is_zero(&self) -> bool {
        self.coeffs.values().all(|x| x.is_zero())
    }

    fn one() -> TLMorphism<R> {
        unimplemented!()
    }

    fn is_one(&self) -> bool {
        unimplemented!()
    }
}

impl<R: Ring + Tex + Signed> Tex for TLMorphism<R> {
    fn into_tex(&self) -> String {
        let mut ans = String::new();
        for (k,v) in self.coeffs.iter() {
            if v.is_zero() { continue; }
            if !ans.is_empty() {
                if v.is_positive() {
                    ans += " + ";
                } else {
                    ans += " - ";
                }
            }
            if ans.is_empty() && v.is_negative() {
                ans += "-";
            }
            if v.is_multiterm() {
                ans += "\\left(";
            }
            if !v.abs().is_one() {
                ans += &v.abs().into_tex();
            }
            if v.is_multiterm() {
                ans += "\\right)";
            }
            ans += "\\,";
            ans += &k.into_tex();
        }
        ans
    }

    fn is_multiterm(&self) -> bool {
        self.coeffs.values().filter(|x| !x.is_zero()).count() > 1
    }
}

impl<R: Ring + Signed + Serialisable> Serialisable for TLMorphism<R>
where for<'r> &'r R : RingOps<&'r R, R> {
    fn serialise(&self) -> String {
        let mut ans = String::new();
        for (k,v) in self.coeffs.iter() {
            ans += &format!("{}|{}",k.serialise(), v.serialise());
            ans += "\n";
        }
        ans
    }

    fn deserialise(inpt : &str) -> Self {
        let mut ans = Vec::new();
        for line in inpt.split("\n") {
            if line.len() == 0 {
                break;
            }
            let parts : Vec<_> = line.split("|").collect();
            assert!(parts.len()==2, "Cannot parse TLMorphism");
            ans.push((TLDiagram::deserialise(parts[0]), R::deserialise(parts[1])));
        }
        TLMorphism::new(ans, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Polynomial;
    use crate::fraction::Fraction;
    use crate::structures::Q;

    #[test]
    fn equality() {
        assert_eq!(TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(3, 3, 0b0100, 0b0100), -1),
        ], None),
        TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(3, 3, 0b0100, 0b0100), -1),
        ], None));
        assert_eq!(TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(3,3,0b100, 0b100), -1),
        ], None),
        TLMorphism::new(vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(3,3,0b100, 0b100), -1),
            (TLDiagram::new(3,3,0b1000, 0b100), 0),
        ], None));
    }

    #[test]
    fn add() {
        assert_eq!(TLMorphism::<i128>::id_zero(5) + TLMorphism::id_zero(5), TLMorphism::id_zero(5));
        assert_ne!(TLMorphism::<i128>::id_zero(5) + TLMorphism::id_zero(5), TLMorphism::id_zero(6));
        assert_eq!(TLMorphism::<i128>::id_zero(5) + TLMorphism::id(5), TLMorphism::id(5));
        assert_eq!(TLMorphism::<i128>::id(5) + TLMorphism::id_zero(5), TLMorphism::id(5));
    }

    #[test]
    fn mul() {
        assert_eq!(
            TLMorphism::id_zero(6),
            TLMorphism::id(6) * (0 as i128)
        );
        type R = Fraction<Polynomial<Q>>;
        let delta = Fraction::from(Polynomial::gen());
        let mut a : TLMorphism<R> = TLMorphism::from(TLDiagram::u(6, 3));
        a.repoint(Some(delta.clone()));
        assert_eq!(
            &a * &a,
            a * delta,
        );
    }
}
