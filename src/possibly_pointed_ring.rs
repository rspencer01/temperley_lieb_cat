use num::Num;

#[derive(Copy, Clone)]
pub struct PossiblyPointedRing<R>
where R : Num + Copy + std::fmt::Debug + std::fmt::Display {
    val : R,
    delta : Option<R>
}
