/// Manually construct $JW_3$ and confirm its properties

use temperley_lieb_cat::fraction::Fraction;
use temperley_lieb_cat::poly::{quantum, Polynomial};
use temperley_lieb_cat::{TLMorphism, TLDiagram};
use temperley_lieb_cat::structures::Q;
use num::One;

#[test]
fn manual_jw3() {
    type R = Fraction<Polynomial<Q>>;
    let delta : R = Polynomial::gen().into();
    let jw3 = TLMorphism::new(vec![
        (TLDiagram::id(3), R::one()),
        (TLDiagram::new(3,3, vec![2], vec![2]), -R::from(quantum(2)) / Fraction::from(quantum(3))),
        (TLDiagram::new(3,3, vec![3], vec![3]), -R::from(quantum(2)) / Fraction::from(quantum(3))),
        (TLDiagram::new(3,3, vec![2], vec![3]), R::from(quantum(1)) / Fraction::from(quantum(3))),
        (TLDiagram::new(3,3, vec![3], vec![2]), R::from(quantum(1)) / Fraction::from(quantum(3))),
    ], Some(Polynomial::gen().into()));
    assert_eq!(jw3, jw3.involute());
    let mut a = TLMorphism::from(TLDiagram::u(3,1));
    let b = TLMorphism::from(TLDiagram::u(3,1));
    a.repoint(Some(delta));
    assert_eq!(a.clone() * jw3.clone(), TLMorphism::id_zero(3));
    assert_eq!(b.clone() * jw3.clone(), TLMorphism::id_zero(3));
    assert_eq!(jw3.clone() * a.clone(), TLMorphism::id_zero(3));
    assert_eq!(jw3.clone() * b.clone(), TLMorphism::id_zero(3));
    assert_eq!(jw3.clone() * jw3.clone(), jw3);
    assert!(jw3.is_jones_wenzl());
}

