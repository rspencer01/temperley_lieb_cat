/// If $n = a + b$ and $e = JW_a \otimes JW_b$ then the algebra $e TL_n e$
/// is semisimple and a polynomial algebra.
use temperley_lieb_cat::fraction::Fraction;
use temperley_lieb_cat::poly::{quantum, Polynomial};
use temperley_lieb_cat::{TLMorphism, jw};
use temperley_lieb_cat::structures::{Ring, Q};

type R = Fraction<Polynomial<Q>>;

#[test]
fn two_jw_semisimple_poly() {
    let n1 : usize = 2;
    let n2 : usize = 4;
    let id = TLMorphism::<R>::id(n1+n2);
    let e : TLMorphism<R> = jw(n1) | jw(n2);

    let factor = R::from(quantum(n1) * quantum(n2));
    let phi = (0..n1+1).map(|i|
            &(
            &e *
            &TLMorphism::big_u(n1+n2, n1, i)) *
            &e *
            &factor
        ).collect::<Vec<_>>();

    for k in 0..n1 as usize {
        assert_eq!(
        &phi[k+1] * R::from(quantum(n1-k) * quantum(n2-k)),
        &(phi[1].clone() - &id * R::from(quantum(k) * quantum(n1+n2-k+1))) * &phi[k]);
    }
    assert!(
        (&(phi[1].clone() - &id * R::from(quantum(n1) * quantum(n2+1))) * &phi[n1]).is_zero()
    );

    let mut min_poly = Polynomial::<R>::one();
    for i in 0..n1+1 {
        min_poly = min_poly * (Polynomial::gen() - R::from(quantum(i) * quantum(n1+n2-i+1)));
    }
    assert!(
        min_poly.eval(phi[1].clone(), id).is_zero()
    );

}
