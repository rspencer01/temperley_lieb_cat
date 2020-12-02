extern crate num;

mod fraction;
mod gcd;
mod poly;
mod temperley;
mod temperley_diagram;
mod temperley_site;
mod temperley_link;
mod tex;

use temperley::jw;
use temperley::TLMorphism;
use poly::{quantum, Polynomial};
use fraction::Fraction;

use tex::Tex;

fn main() {
    type R = Fraction<Polynomial<i128>>;
    let n1 : i128 = 3;
    let n2 : i128 = 4;
    let e = jw(n1 as usize) | jw(n2 as usize);
    let mut phi = Vec::new();
    for i in 0..n1 as usize+1 {
        println!("Constructing \\phi_{}", i);
        phi.push(e.clone() *TLMorphism::U((n1+n2) as usize, n1 as usize, i) * e.clone() * (Fraction::from(quantum(n1)) * quantum(n2)));
    }

    for k in 0..n1 as usize {
        assert_eq!(
        phi[k+1].clone() * Fraction::from(quantum(n1-k as i128) * quantum(n2-k as i128)),
        (phi[1].clone() - TLMorphism::<R>::id((n1+n2) as usize) * Fraction::from(quantum(k as i128) * quantum(n1+n2-k as i128+1))) * phi[k].clone());
    }
}
