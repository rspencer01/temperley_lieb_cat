use temperley_lieb_cat::jw;
use temperley_lieb_cat::quantum;
use temperley_lieb_cat::structures::Ring;
use temperley_lieb_cat::structures::{RingOps, Q};
use temperley_lieb_cat::Fraction;
use temperley_lieb_cat::Polynomial;
use temperley_lieb_cat::TLDiagram;
use temperley_lieb_cat::TLMorphism;
use temperley_lieb_cat::Tex;

type R = Fraction<Polynomial<Q>>;

fn ladder<F>(path: &[i128], g: F) -> TLMorphism<R>
where
    F: Fn(usize) -> TLMorphism<R>,
{
    let mut ans = jw(0);
    let mut r = 0;
    for n in 0..path.len() {
        ans = ans | TLMorphism::id(i128::abs(path[n]) as usize);
        if path[n] < 0 {
            ans = &ans
                * &TLMorphism::from(TLDiagram::big_cap(
                    ans.co_domain(),
                    r,
                    i128::abs(path[n]) as usize,
                ));
            r -= i128::abs(path[n]) as usize;
        } else {
            r += i128::abs(path[n]) as usize;
        }
        println!("{} {} {}->{}", n, r, ans.domain(), ans.co_domain());
        assert!(r == ans.co_domain());
        ans = ans * g(r);
    }
    ans
}

fn quotient_non_monic(f: TLMorphism<R>) -> TLMorphism<R> {
    TLMorphism::new(
        f.coeffs
            .iter()
            .filter(|(d, _)| d.propagation() == d.co_domain())
            .map(|(k, v)| (*k, v.clone()))
            .collect::<Vec<_>>(),
        f.delta,
    )
}

fn main() {
    jw(3).trace().render_tex();
    return;
    //ladder(vec![1, 1], |n| {
    //    TLMorphism::id(n)
    //})
    //.render_tex();
    let d = ladder(vec![3, -1], |n| jw(n));
    d.render_tex();
    quotient_non_monic(d).render_tex();
}
