use temperley_lieb_cat::structures::Ring;
/// Confirm the form of jw3 over characteristic 2
///
/// This can be done over integers
use temperley_lieb_cat::{TLDiagram, TLMorphism};

#[test]
fn test_jw3over2() {
    let jw3lp = TLMorphism::<i128>::new(
        vec![
            (TLDiagram::id(3), 1),
            (TLDiagram::new(3, 3, 0b0100, 0b1000), -1),
            (TLDiagram::new(3, 3, 0b1000, 0b0100), -1),
        ],
        Some(0),
    );
    assert!((&TLMorphism::u(3, 1) * &jw3lp).is_zero());
    assert!((&TLMorphism::u(3, 2) * &jw3lp).is_zero());
    assert!(jw3lp == jw3lp.involute());
}
