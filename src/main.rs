#![allow(dead_code)]
extern crate num;


#[macro_use]
mod macros;

mod fraction;
mod gcd;
mod poly;
mod temperley;
mod temperley_diagram;
mod temperley_site;
mod temperley_link;
mod tex;
mod serial;
mod structures;

use temperley::jw;
use temperley::TLMorphism;
use temperley_diagram::TLDiagram;
use poly::{Polynomial,quantum};
use fraction::Fraction;

use structures::{NumOps, Q};

use tex::Tex;
use num::{One,Zero};
use temperley_site::Site::*;
use serial::Serialisable;

fn reduce_mod<R>(a : TLMorphism<Fraction<R>>,p : R) -> TLMorphism<Fraction<R>> 
where R : gcd::PartialGCD + num::Signed + Clone  + NumOps<R, R> {
    TLMorphism::new(
        a.coeffs.iter()
        .map(|(k,v)| (k.clone(), Fraction::new(v.num() % p, v.den() % p)))
        .collect::<Vec<_>>(),
        a.delta,
    )
}

fn jw5over3() -> TLMorphism<Fraction<Polynomial<Q>>> {
    let jw2 =jw(2);
    let mut t = TLMorphism::id(5) * Fraction::zero();
    for i in -2..3 {
        let coeff = if i %2 == 0 { Fraction::one() } else { - Fraction::one() };
        t = t + (jw2.turn_down(i) | TLMorphism::id(1) | jw2.turn_up(-i) * coeff);
    }
    t
}

fn jw7over4() -> TLMorphism<Fraction<Polynomial<Q>>>{
    let jw3 =jw(3);
    let mut t = TLMorphism::id_zero(7);
    for i in -3..4 {
        let coeff = if i % 2 == 0 { Fraction::one() } else { - Fraction::one() };
        t = t + (jw3.turn_down(i) | TLMorphism::id(1) | jw3.turn_up(-i) * coeff);
    }
    t
}

fn make_and_save_jw(n : usize) {
    let dest = std::path::PathBuf::from(format!("./cache/jw{}",n));
    let jw = jw(n);
    assert!(jw.is_jones_wenzl());
    jw.to_file(&dest);
}

fn main() {
    foo2();
    return
    let mut jw10 : TLMorphism<Fraction<Polynomial<Q>>> = TLMorphism::from_file(&std::path::PathBuf::from("./cache/jw10"));
    println!(".");
    jw10.repoint(Some(Polynomial::gen().into()));
    assert!(jw10.is_jones_wenzl());
    println!(".");
    jw10 = jw10 | TLMorphism::id(1);
    println!(".");
    let jw11 = &jw10 * &TLMorphism::u(10,10);
    println!(".");
    let jw11bar = TLMorphism::new(
        jw10.coeffs.iter()
        .filter(|(k,_)|
                !(k.link(Source(1)) == Source(2))
            )
        .filter(|(k,_)|
                !(k.link(Source(2)) == Source(3))
            )
        .filter(|(k,_)|
                !(k.link(Source(3)) == Source(4))
            )
        .filter(|(k,_)|
                !(k.link(Source(4)) == Source(5))
            )
        .filter(|(k,_)|
                !(k.link(Source(5)) == Source(6))
            )
        .filter(|(k,_)|
                !(k.link(Source(6)) == Source(7))
            )
        .filter(|(k,_)|
                !(k.link(Source(7)) == Source(8))
            )
        .map(|(k,v)| (k.clone(), v.clone()))
        .collect::<Vec<_>>()
        ,
        Some(Polynomial::gen().into())
    );
    println!("{}", jw11bar.support().len());
    let jw11bar = jw11 * jw11bar;
    println!(".");
    let jw11 = &jw10 - &(jw11bar * Fraction::new(quantum(10), quantum(10)));
    println!(".");
    println!("{}", jw11.serialise());
    assert!(jw11.is_jones_wenzl());
}
fn foo2() {
    let jw3lp = TLMorphism::new(vec![
        (TLDiagram::id(3), Fraction::one()),
        (TLDiagram::from_tableauxs(3, vec![2].into_iter(), vec![3].into_iter()), -Fraction::one()),
        (TLDiagram::from_tableauxs(3, vec![3].into_iter(), vec![2].into_iter()), -Fraction::one())
    ], Some(Fraction::from(quantum(2))));
    let jw5 = TLMorphism::from_file(&std::path::PathBuf::from("./cache/jw5"));
    let mut jwlp = reduce_mod(jw5, quantum(2));
    println!("jwlp has {} terms", jwlp.support().len());
    let half : Fraction<Polynomial<Q>> = Fraction::new(
        Polynomial::one(),
        Polynomial::<Q>::one() * Q::from(2)
    );
    let mut candidate = TLMorphism::id_zero(5);
    for i in -1..2 {
        let c = if i%2==0 {Fraction::one()} else {-Fraction::one()};
        candidate = candidate + (
            jw(1).turn_down(i) |
            TLMorphism::id(1) |
            jw3lp.turn_up(-i)
        ) * c;
        candidate = candidate + (
            jw3lp.turn_down(-i) |
            TLMorphism::id(1) |
            jw(1).turn_up(i)
        ) * c;
    }
    for i in 0..2 {
        let a =
            TLMorphism::id(i)|
            (
                TLMorphism::id(1) |
                jw(1).rotate(1+i as isize)
            ).turn_down(2) |
            TLMorphism::id(1-i);
        for j in 0..2 {
            let c = if (i+j)%2==0 {Fraction::one()} else {-Fraction::one()};
            let b =
                TLMorphism::id(j)|
                (
                    jw(1).rotate(j as isize)|
                    TLMorphism::id(1)
                ).turn_up(-2) |
                TLMorphism::id(1-j);
            candidate = candidate + a.clone() * jw(1) * b * c;
        }
    }
    (TLMorphism::u(5,2) * candidate.clone()).render_tex();
    jwlp = jwlp - candidate * half;
    println!("after candidate subtracted, has {} terms", jwlp.support().len());
    jwlp.render_tex();
}

/*
fn foo4() {
    //p^r = 4
    let jw3 = jw(3);
    let start = std::time::Instant::now();
    let jw11 = jw(10);
    println!("Done 11 in {}s", start.elapsed().as_secs());
    return;
    let delta = jw3.delta;

    let factor : Fraction<Polynomial<Q>> = Fraction::new(
        Polynomial::<Q>::one() * Q::from(2),
        Polynomial::one()
    );

    let jw7 = reduce_mod(jw7over4(), quantum(4));
    println!("{}", jw7.support().len());

    let half : Fraction<Polynomial<Q>> = Fraction::new(
        Polynomial::one(),
        Polynomial::<Q>::one() * Q::from(2)
    );
    let mut candidate : TLMorphism<Fraction<Polynomial<Q>>> = TLMorphism::id_zero(11);

    for i in -3..4 {
        let c = if i%2==0 {Fraction::one()} else {-Fraction::one()};
        candidate = candidate + (
            jw3.turn_down(-i) |
            TLMorphism::id(1) |
            jw7.turn_up(i)
        ) * c;
        candidate = candidate + (
            jw7.turn_down(i) |
            TLMorphism::id(1) |
            jw3.turn_up(-i)
        ) * c;
        print!(".");
    }

    for i in 0..4 {
        let a =
            TLMorphism::id(i)|
            (
                TLMorphism::id(1) |
                jw3.rotate(1+i as isize)
            ).turn_down(4) |
            TLMorphism::id(3-i);
        for j in 0..4 {
            let c = if (i+j)%2==0 {Fraction::one()} else {-Fraction::one()};
            let b =
                TLMorphism::id(j)|
                (
                    jw3.rotate(j as isize)|
                    TLMorphism::id(1)
                ).turn_up(-4) |
                TLMorphism::id(3-j);
            candidate = candidate + a.clone() * jw3.clone() * b * c;
        }
        print!(".");
    }
    print!("\r           ");
    let temp = reduce_mod(TLMorphism::u(11,3) * candidate.clone(),quantum(3));
    println!("{}", temp.support().len());
    temp.render_tex();

}
*/
/*
fn foo3() {
    let jw2 = jw(2);
    let jw8 = jw(8);
    let mut jwlp = reduce_mod(jw8,quantum(3));
    println!("{}", jwlp.support().len());
    let delta = jwlp.delta;

//    for k in -2..3 {
//        for i in (-2+k).max(-2)..(3+k).min(3) {
//            let coeff = if (i + k + 4) % 2 == 1 {Fraction::one()} else {-Fraction::one()};
//            let x = jw2.clone().turn_down(k) |
//                TLMorphism::id(1) | 
//                jw2.clone().turn_down(-k).turn_up(i) | 
//                TLMorphism::id(1) | 
//                jw2.clone().turn_up(-i);
//            jwlp = jwlp + x * coeff;
//        }
//        jwlp = reduce_mod(jwlp, quantum(3));
//        println!("{}",jwlp.support().len());
//    }

//    jwlp = TLMorphism::new(
//        jwlp.coeffs.iter()
//        .filter(|(k,v)| {
//            k.link(Source(4)) == Target(2)
//        })
//        .map(|(k,v)| (k.clone(), *v))
//        .collect::<Vec<_>>(),
//        delta
//    );
//
//    let cappy_boi = jw(2).turn_down(2);
//    let cappy_cap_boi = (TLMorphism::id(1) | cappy_boi.clone()).turn_down(1);
//    assert_eq!(cappy_cap_boi.domain(), 6);
//    assert_eq!(cappy_cap_boi.co_domain(), 0);
//
//    let cuppy_boi = jw(2).turn_down(-2);
//    let cuppy_cup_boi = (TLMorphism::id(1) | cuppy_boi.clone()).turn_down(-1);
//    assert_eq!(cuppy_cup_boi.co_domain(), 6);
//    assert_eq!(cuppy_cup_boi.domain(), 0);
//
//    let coeff = Fraction::new(quantum(1), quantum(1)*2);
//    let x = jw2.clone().turn_down(1) |
//        TLMorphism::id(1) |
//        cuppy_cup_boi |
//        cappy_boi;
//    jwlp = jwlp + x * coeff;
//    jwlp = reduce_mod(jwlp, quantum(3));
//    jwlp.render_tex();
//    println!("{}",jwlp.support().len());

//    let mut ans = TLMorphism::new(vec![(TLDiagram::any(8,8),Fraction::zero())], None);
//    for i in 0..16 {
//        for j in i+6..16 {
//            if (j+5)%16 == i { break;}
//            if (3 <= i) && (i<8) && (11<=j) { continue; }
//            println!("(i,j) = {:?}", (i,j));
//            let mut done = vec![false ; 16];
//            let mut links = Vec::new();
//            links.push(Link::new(Source(i+1), Source((i+5)%16 + 1)));
//            links.push(Link::new(Source((i+1)%16+1), Source((i+4)%16 + 1)));
//            links.push(Link::new(Source((i+2)%16+1), Source((i+3)%16 + 1)));
//            for k in 0..6 {
//                done[(i+k)%16] = true;
//            }
//
//
//            links.push(Link::new(Source(j+1), Source((j+5)%16 + 1)));
//            links.push(Link::new(Source((j+1)%16+1), Source((j+4)%16 + 1)));
//            links.push(Link::new(Source((j+2)%16+1), Source((j+3)%16 + 1)));
//            for k in 0..6 {
//                done[(j+k)%16] = true;
//            }
//            let mut stack = Vec::new();
//            let mut k = (j+6)%16;
//            while stack.len() < 3 - (i%2)  {
//                if !done[k] {
//                    stack.push(k);
//                }
//                k = (k+1)%16;
//            }
//            while stack.len() > 1-(i%2) {
//                if !done[k] {
//                    let s = stack.pop().unwrap();
//                    links.push(Link::new(Source(k+1), Source(s+1)));
//                }
//                k = (k+1)%16;
//            }
//            let diag = TLDiagram::new(links).turn_down(-8);
//            diag.render_tex();
//            let correction = if (1+i+j)%2==1 {-quantum(1)*2} else {quantum(1)*2};
//            println!("{}", jwlp.coeffs[&diag]*correction);
//            println!();
//            ans = ans + TLMorphism::new(vec![(diag.clone(),jwlp.coeffs[&diag]*correction)], None);
//        }
//    }
//    ans.render_tex();
//    println!("{}", ans.support().len());
    let factor : Fraction<Polynomial<Q>> = Fraction::new(
        Polynomial::<Q>::one() * Q::from(2),
        Polynomial::one()
    );

    let jw5 = jw5over3();

    let half : Fraction<Polynomial<Q>> = Fraction::new(
        Polynomial::one(),
        Polynomial::<Q>::one() * Q::from(2)
    );
    let mut candidate : TLMorphism<Fraction<Polynomial<Q>>> = TLMorphism::id_zero(8);

    //println!("Before, jwlp has {} terms", jwlp.support().len());

    for i in -2..3 {
        let c = if i%2==0 {Fraction::one()} else {-Fraction::one()};
        candidate = candidate + (
            jw2.turn_down(-i) |
            TLMorphism::id(1) |
            jw5.turn_up(i)
        ) * c;
        candidate = candidate + (
            jw5.turn_down(i) |
            TLMorphism::id(1) |
            jw2.turn_up(-i)
        ) * c;
    }

    for i in 0..3 {
        let a =
            TLMorphism::id(i)|
            (
                TLMorphism::id(1) |
                jw2.rotate(1+i as isize)
            ).turn_down(3) |
            TLMorphism::id(2-i);
        for j in 0..3 {
            let c = if (i+j)%2==0 {-Fraction::one()} else {Fraction::one()};
            let b =
                TLMorphism::id(j)|
                (
                    jw2.rotate(j as isize)|
                    TLMorphism::id(1)
                ).turn_up(-3) |
                TLMorphism::id(2-j);
            candidate = candidate + a.clone() * jw2.clone() * b * c;
        }
    }
    let x = reduce_mod(TLMorphism::u(8,3) * candidate.clone(),quantum(3));
    x.render_tex();

    jwlp = jwlp - candidate * half;

    jwlp = reduce_mod(jwlp, quantum(3));

    println!("After, jwlp has {} terms", jwlp.support().len());
    jwlp = TLMorphism::new(
        jwlp.coeffs.iter()
        .filter(|(k,v)| {
            k.link(Source(8)) == Source(5) ||
            (k.link(Source(8)) == Source(7) &&
             k.link(Source(6)) == Source(5))
        })
        .map(|(k,v)| (k.clone(), *v))
        .collect::<Vec<_>>(),
        delta
    );
    println!("{}", jwlp.support().len());

    for i in 1..8 {
        println!("{} {}",i,
                 reduce_mod(TLMorphism::u(8,i) * jwlp.clone(), quantum(3)).support().len());
    }
    //jwlp = reduce_mod(TLMorphism::u(8,5) * jwlp, quantum(3));

    jwlp.render_tex();
}*/
