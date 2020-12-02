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

fn main() {
    let j = jw(4);
    println!("{:?}", j);
    println!("{}", j.support().len());
    println!("{}", j.is_jones_wenzl());
}
