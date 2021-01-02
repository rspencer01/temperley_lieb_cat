/// If $n = a + b$ and $e = JW_a \otimes JW_b$ then the algebra $e TL_n e$
/// is semisimple and a polynomial algebra.
use temperley_lieb_cat::fraction::Fraction;
use temperley_lieb_cat::poly::{quantum, Polynomial};
use temperley_lieb_cat::{TLMorphism, jw};
use temperley_lieb_cat::structures::Q;
use num::Zero;

type R = Fraction<Polynomial<Q>>;

#[test]
fn two_jw_semisimple_poly() {
    let n1 : i128 = 3;
    let n2 : i128 = 3;
    let e : TLMorphism<R> = jw(n1 as usize) | jw(n2 as usize);
    let mut phi = Vec::new();
    for i in 0..n1 as usize+1 {
        phi.push(
            e.clone() *
            TLMorphism::big_u((n1+n2) as usize, n1 as usize, i) *
            e.clone() *
            R::from(quantum(n1) * quantum(n2)));
    }

    for k in 0..n1 as usize {
        assert_eq!(
        phi[k+1].clone() * Fraction::from(quantum(n1-k as i128) * quantum(n2-k as i128)),
        (phi[1].clone() - TLMorphism::<R>::id((n1+n2) as usize) * R::from(quantum(k as i128) * quantum(n1+n2-k as i128+1))) * phi[k].clone());
    }
    assert!(
        ((phi[1].clone() - TLMorphism::<R>::id((n1+n2) as usize) * R::from(quantum(n1) * quantum(n2+1))) * phi[n1 as usize].clone()).is_zero()
    );

    let mut t = TLMorphism::id((n1+n2) as usize);
    for i in 0..n1+1 {
        t = t * (phi[1].clone() - TLMorphism::<R>::id((n1+n2) as usize) * R::from(quantum(i) * quantum(n1+n2-i+1)));
    }
    assert!(t.is_zero());
}
