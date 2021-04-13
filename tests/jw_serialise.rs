use temperley_lieb_cat::jw;
use temperley_lieb_cat::Serialisable;
/// Test serialisation of large elements
use temperley_lieb_cat::TLMorphism;

#[test]
fn jw8_serialise() {
    let jw = jw(8);
    let s = jw.serialise();
    let jwp = TLMorphism::deserialise(&s);
    assert_eq!(jw, jwp);
}
