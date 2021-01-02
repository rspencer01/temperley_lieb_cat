/// Confirm the form of jw3 over characteristic 2
///
/// This can be done over integers

use temperley_lieb_cat::{TLMorphism, TLDiagram};
use num::Zero;

#[test]
fn test_jw3over2() {
    let jw3lp = TLMorphism::<i128>::new(vec![
        (TLDiagram::id(3), 1),
        (TLDiagram::from_tableauxs(3, vec![2].into_iter(), vec![3].into_iter()), -1),
        (TLDiagram::from_tableauxs(3, vec![3].into_iter(), vec![2].into_iter()), -1)
    ], Some(0));
    assert!((&TLMorphism::u(3,1) * &jw3lp).is_zero());
    assert!((&TLMorphism::u(3,2) * &jw3lp).is_zero());
    assert!(jw3lp == jw3lp.involute());
}

