use crate::poly::{quantum, Polynomial};
use crate::temperley::TLMorphism;
use crate::temperley_diagram::TLDiagram;
use crate::fraction::Fraction;
use crate::structures::Q;
use std::collections::HashMap;
use num::{One, Zero};
use std::sync::Mutex;

lazy_static! {
    static ref JW_COEFFS : Mutex<HashMap<TLDiagram, Fraction<Polynomial<Q>>>> = {
        let mut m = HashMap::new();
        m.insert(TLDiagram::new(0,0,0,0), Fraction::one());
        m.insert(TLDiagram::new(1,1,0,0), Fraction::one());
        Mutex::new(m)
    };
    static ref JW_COMPUTED : Mutex<usize> = Mutex::new(1);
}

/// Construct the Jones-Wenzl element on a certain number of strands
/// 
/// This function uses Morrison's formula for coefficients and caches
/// the results so succesive calls will be faster
pub fn jw(n : usize)-> TLMorphism<Fraction<Polynomial<Q>>>  {
    let x = Polynomial::gen().into();
    for m in (*JW_COMPUTED.lock().unwrap()+1)..(n+1) {
        for d in TLDiagram::all(m,m) {
            let mut ans = Fraction::zero();
            let dp = d.turn_down(1);
            for i in dp.simple_links().into_iter() {
                let coeff = if (i + m) % 2 == 0 {
                    Fraction::one()
                } else {
                    -Fraction::one()
                };
                ans = ans + &JW_COEFFS.lock().unwrap()[&dp.remove_cap(i)]
                    * &Fraction::new(quantum(-(i as isize)) , quantum(-(m as isize))) * &coeff;
            }
            JW_COEFFS.lock().unwrap().insert(d,ans);
        }
    }
    if n > *JW_COMPUTED.lock().unwrap() {
        *JW_COMPUTED.lock().unwrap() = n;
    }
    let mut ans = TLMorphism::new(TLDiagram::all(n,n).iter()
                    .map(|d| (*d, JW_COEFFS.lock().unwrap()[d].clone()))
                    .collect::<Vec<_>>(),
                    None);
    ans.right_kills.extend(1..n);
    ans.delta = Some(x);
    ans
}

/// Construct the Jones-Wenzl element on a certain number of strands
/// 
/// This function computes the idempotents using the usual iterative
/// form.
pub fn jw_old(n : usize) -> TLMorphism<Fraction<Polynomial<Q>>> {
    let mut jw = TLMorphism::id(1);
    jw.repoint(Some(Polynomial::gen().into()));
    for i in 1..n {
        let jwa = jw.clone() | TLMorphism::id(1);
        let jwb = jw.turn_down(1) * jw.turn_down(1).involute();
        let jwb = jwb * Fraction::new(quantum(i), quantum(i+1));
        jw = jwa - jwb;
    }
    jw
}

fn l_p_adic_digits(n : usize, l : usize, p : usize) -> Vec<usize> {
    let mut ans = Vec::new();
    if n > 0 {
        ans.push(n % l);
        ans.extend(l_p_adic_digits(n / l, p, p));
    }
    ans
}

fn jw2lp(n : usize, l : usize, p : usize) -> TLMorphism<Fraction<Polynomial<Q>>> {
    let jwn = l_p_reduce_mod(jw(n), l, p);
    let mut t = TLMorphism::id(2 * n + 1) * Fraction::zero();
    let n = n as isize;
    for i in -n..n+1 {
        let coeff = if i %2 == 0 { Fraction::one() } else { - Fraction::one() };
        t = t + (jwn.turn_down(i) | TLMorphism::id(1) | jwn.turn_up(-i) * coeff);
    }
    t
}

fn reduce_mod(a : TLMorphism<Fraction<Polynomial<Q>>>, p : Polynomial<Q>) -> TLMorphism<Fraction<Polynomial<Q>>> {
    let foo =
        a.coeffs.into_iter()
        .map(|(k,v)|
            (k, Fraction::new(v.num() % &p, v.den() % &p))
        )
        .collect::<Vec<_>>();
    TLMorphism::new(
        foo,
        a.delta,
    )
}

fn l_p_reduce_mod(a : TLMorphism<Fraction<Polynomial<Q>>>, l : usize, p : usize) -> TLMorphism<Fraction<Polynomial<Q>>> {
    reduce_mod(a, quantum(l))
}


pub fn jwlp(n : usize, l : usize, p : usize) -> TLMorphism<Fraction<Polynomial<Q>>> {
    let digits = l_p_adic_digits(n + 1, l, p);
    if digits.iter().filter(|x| **x != 0).count() == 1 {
        // We are eve
        if digits.iter().sum::<usize>() == 2 {
            jw2lp((n + 1) / 2 - 1, l, p)
        } else {
            l_p_reduce_mod(jw(n), l, p)
        }
    } else {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use num::Zero;

    #[test]
    fn general_jw() {
        for i in 1..7 {
            let mut jw = jw(i);
            jw.repoint(jw.delta.clone());
            assert!(jw.is_jones_wenzl());
        }
        assert!((jw(4) | TLMorphism::id(1)) * jw(5) == jw(5));
        assert!(jw(6).flip() == jw(6));
        assert!((jw(4) | jw(2)).flip() == (jw(2) | jw(4)));
    }

    #[test]
    fn tensored_jw() {
        assert!((jw(4) | jw(3)).is_idempotent());
        assert_eq!((jw(4) | jw(3)).right_kills, vec![1,2,3,5,6]);
        assert!((jw(3) | jw(2) | jw(3)).is_idempotent());
        assert!(!((jw(3) | jw(2) | jw(3)).is_jones_wenzl()));
        assert!((jw(5) | jw(1)).is_idempotent());
        assert!(!(jw(5) | jw(1)).is_jones_wenzl());
        assert!((jw(4) | jw(2)) * jw(6) == jw(6));
    }

    #[test]
    fn jw5over3() {
        for p in &[2, 3, 5, 7, 11] {
            let t = jwlp(5, 3, *p);
            for i in 1..5 {
                let s = &TLMorphism::u(5,i) * &t;
                for c in s.coeffs.values() {
                    assert!((c.num() % &quantum(3)).is_zero());
                }
            }
        }
    }

    #[test]
    fn jw11over2_3() {
        let t = jwlp(11, 2, 3);
        assert!(!t.is_zero());
        for i in 1..11 {
            let s = &TLMorphism::u(11,i) * &t;
            for c in s.coeffs.values() {
                let d = c.num() % &quantum(2);
                for i in 0..d.degree()+1 {
                    let e = d.coeff(i);
                    assert!(e.is_integral());
                    assert!(d.coeff(i).num() % 3 == 0);
                }
            }
        }
    }

    #[test]
    fn jw9over5() {
        let jw6 = jw(6) | TLMorphism::id(1);
        let jw7 = &jw6 * &TLMorphism::u(7,6);
        let jw7 = &jw7 * &jw6;
//        println!("{}", jw7.serialise());
        let jw7 = jw7 * Fraction::new(quantum(6), quantum(7));
//        println!("{}", jw7.serialise());
//        let mut t = TLMorphism::id(9) * Fraction::zero();
//        for i in -4..5 {
//            let coeff = if i %2 == 0 { Fraction::one() } else { - Fraction::one() };
//            t = t + (jw4.turn_down(i) | TLMorphism::id(1) | jw4.turn_up(-i) * coeff);
//        }
//        for i in 1..9 {
//            let s = TLMorphism::u(9,i) * t.clone();
//            for c in s.coeffs.values() {
//                assert!((c.num() % quantum(5)).is_zero());
//            }
//        }
    }

    #[test]
    fn p_adic() {
        assert_eq!(vec![2,3,0,1usize], l_p_adic_digits(2+3*5+1*125, 5, 5));
    }

    #[test]
    fn morrison() {
        assert_eq!(jw(3), jw_old(3));
        assert_eq!(jw(7), jw_old(7));
    }
}
