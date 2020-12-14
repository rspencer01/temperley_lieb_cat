use temperley_lieb_cat::fraction::Fraction;
use temperley_lieb_cat::poly::{quantum, Polynomial};
use temperley_lieb_cat::temperley::{TLMorphism, jw};
use num::Zero;

#[test]
fn two_jw_semisimple_poly() {
    type R = Fraction<Polynomial<i128>>;
    let n1 : i128 = 3;
    let n2 : i128 = 3;
    let e = jw(n1 as usize) | jw(n2 as usize);
    let mut phi = Vec::new();
    for i in 0..n1 as usize+1 {
        println!("Constructing \\phi_{}", i);
        phi.push(e.clone() *TLMorphism::U((n1+n2) as usize, n1 as usize, i) * e.clone() * (R::from(quantum(n1)) * quantum(n2)));
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
        println!("Adding term {}", i);
        t = t * (phi[1].clone() - TLMorphism::<R>::id((n1+n2) as usize) * R::from(quantum(i) * quantum(n1+n2-i+1)));
    }
    assert!(t.is_zero());
}
