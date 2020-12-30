extern crate partitions;

use crate::temperley_site::{Site, Site::*};
use crate::temperley_link::Link;
use crate::tex::Tex;
use crate::serial::Serialisable;

#[derive(Clone, PartialOrd, Ord)]
pub struct TLDiagram(Vec<Link>);

impl std::fmt::Debug for TLDiagram {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl std::fmt::Display for TLDiagram {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f,"{:?}", self.0)
    }
}

impl std::hash::Hash for TLDiagram {
    fn hash<H:std::hash::Hasher>(&self, state : &mut H) {
        for i in self.0.iter() {
            i.hash(state)
        }
    }
}

impl TLDiagram {
    pub fn new(links : Vec<Link>) -> TLDiagram {
        let mut lnk = links.clone();
        lnk.sort();
        TLDiagram(lnk)
    }

    pub fn from_tableaux<I>(n : usize, tab : I) -> TLDiagram
    where I : Iterator<Item=usize> {
        let mut stack = Vec::new();
        let mut links = Vec::new();
        let tab = tab.collect::<Vec<_>>();
        for i in 1..n+1 {
            if tab.contains(&i) {
                links.push(Link::new(Source(stack.pop().unwrap()), Source(i)));
            } else {
                stack.push(i);
            }
        }
        for i in 0..stack.len() {
            links.push(Link::new(Source(stack[i]), Target(i+1)));
        }
        TLDiagram::new(links)
    }

    pub fn from_tableauxs<I>(n : usize, tab1 : I, tab2 : I) -> TLDiagram
    where I : Iterator<Item=usize> {
        (TLDiagram::from_tableaux(n, tab1) * TLDiagram::from_tableaux(n, tab2).involute()).1
    }

    pub fn u(n :usize, i : usize) -> TLDiagram {
        TLDiagram::from_tableauxs(n, i+1..i+2, i+1..i+2)
    }

    pub fn U(n :usize, i : usize, j : usize) -> TLDiagram {
        TLDiagram::from_tableauxs(n, i+1..i+1+j, i+1..i+1 + j)
    }

    pub fn domain(&self) -> usize {
        self.0.iter()
            .map(|x|
            match x {
                Link(Source(_), Source(_)) => 2,
                Link(Source(_), Target(_)) => 1,
                Link(Target(_), Source(_)) => 1,
                Link(Target(_), Target(_)) => 0,
            }
            ).sum()
    }

    pub fn co_domain(&self) -> usize {
        self.0.len() * 2 - self.domain()
    }

    pub fn propagation(&self) -> usize {
        self.0.iter()
            .map(|x|
            match x {
                Link(Source(_), Source(_)) => 0,
                Link(Source(_), Target(_)) => 1,
                Link(Target(_), Source(_)) => 1,
                Link(Target(_), Target(_)) => 0,
            }
            ).sum()
    }

    pub fn involute(self) -> TLDiagram {
        TLDiagram::new(
            self.0.into_iter().map(|x| {
                x.involute()
            }).collect()
        )
    }

    pub fn flip(&self) -> TLDiagram {
        TLDiagram::new(
            self.0.iter().map(|x|
                match x {
                    Link(Source(a), Source(b)) => Link::new(Source(self.domain() + 1 - a), Source(self.domain() + 1 - b)),
                    Link(Source(a), Target(b)) => Link::new(Source(self.domain() + 1 - a), Target(self.co_domain() + 1 - b)),
                    Link(Target(a), Target(b)) => Link::new(Target(self.co_domain() + 1 - a), Target(self.co_domain() + 1 - b)),
                    _ => panic!("Unexpected form of diagram"),
                }
            ).collect()
        )
    }

    pub fn cap(n: usize, i: usize) -> TLDiagram {
        TLDiagram::from_tableaux(n, i+1..i+2)
    }

    pub fn cup(n: usize, i: usize) -> TLDiagram {
        TLDiagram::cap(n,i).involute()
    }

    pub fn id(n: usize) -> TLDiagram {
        TLDiagram::new((1..n+1).map(|j| Link::new(Source(j), Target(j))).collect())
    }

    pub fn inject(&self) -> TLDiagram {
        let mut links = self.0.clone();
        links.push(Link::new(Source(self.domain()+1), Target(self.co_domain()+1)));
        TLDiagram::new(links)
    }

    pub fn any(n : usize, m : usize) -> TLDiagram {
        if n < m {
            TLDiagram::any(m,n).involute()
        } else {
            let mut v = Vec::new();
            for i in 0..(n - m)/2 {
                v.push(i*2 + 2);
            }
            TLDiagram::from_tableaux(n,v.into_iter())
        }
    }

    pub fn link(&self, site: Site) -> Site {
        for Link(a, b) in self.0.iter() {
            if *a == site { return *b }
            if *b == site { return *a }
        }
        panic!("Site {:?} not found",site)
    }

    pub fn simple_links(&self) -> Vec<usize> {
        (1..self.domain())
        .filter(|i| self.link(Source(*i)) == Source(i+1))
        .collect()
    }


    pub fn turn_down(&self, i : isize) -> TLDiagram {
        if i == 0 {
            return self.clone();
        }
        let n = self.domain();
        let m = self.co_domain();
        let links =
            if i > 0 {
                assert!(m > 0,"Cannot turn element");
                self.0.iter()
                    .map(|x|
                         if x.0 == Target(m) {
                             Link::new(x.1, Source(n+1))
                         } else if x.1 == Target(m) {
                             Link::new(x.0, Source(n+1))
                         } else {
                             *x
                         }
                    ).collect::<Vec<_>>()
            } else {
                assert!(n > 0, "Cannot turn element");
                self.0.iter()
                    .map(|x|
                         if x.0 == Source(n) {
                             Link::new(x.1, Target(m+1))
                         } else if x.1 == Source(n) {
                             Link::new(x.0, Target(m+1))
                         } else {
                             *x
                         }
                    ).collect::<Vec<_>>()
            };
        let new = if i < 0 { i + 1 } else { i - 1 };
        TLDiagram::new(links).turn_down(new)
    }

    pub fn turn_up(&self, i: isize) -> TLDiagram {
        self.flip().turn_down(i).flip()
    }

    /// Rotate the diagram anticlockwise
    pub fn rotate(&self, i : isize) -> TLDiagram {
        if i == 0 {
            self.clone()
        } else {
            if i > 0 {
                self.turn_up(1).turn_down(-1).rotate(i-1)
            } else {
                self.turn_up(1).turn_down(-1).rotate(i+1)
            }
        }
    }
}

impl PartialEq for TLDiagram {
    fn eq(&self, other : &TLDiagram) -> bool {
        if (self.domain(), self.co_domain()) != (other.domain(), other.co_domain()) {
            return false
        }
        for i in 1..self.domain()+1 {
            if self.link(Source(i)) != other.link(Source(i)) {
                return false;
            }
        }
        for i in 1..self.co_domain()+1 {
            if self.link(Target(i)) != other.link(Target(i)) {
                return false;
            }
        }
        true
    }
}

impl Eq for TLDiagram {
}

impl std::ops::Mul for TLDiagram {
    type Output = (usize, TLDiagram);

    fn mul(self, other: TLDiagram) -> (usize, TLDiagram) {
        assert_eq!(self.co_domain(), other.domain(), "Cannot multiply diagrams if their (co)domains don't match");
        // Inlcusive
        // Sources 1   .. a
        // Middle  a+1 .. a+b
        // Target  a+b .. a+b+c
        let left  = |x : usize| {x};
        let mid   = |x : usize| {x + self.domain()};
        let right = |x : usize| {x + self.domain() + self.co_domain()};
        let mut uf : partitions::PartitionVec<usize> = partitions::PartitionVec::new();
        for i in 0..self.domain() + self.co_domain() + other.co_domain() + 1 {
            uf.push(i);
        }
        for i in 1..self.domain()+1 {
            match self.link(Source(i)) {
                Source(j) => {uf.union(left(i), left(j))},
                Target(j) => {uf.union(left(i), mid(j))},
            }
        }
        for i in 1..self.co_domain() + 1 {
            match self.link(Target(i)) {
                Target(j) => {uf.union(mid(i), mid(j))},
                _ => {},
            }
            match other.link(Source(i)) {
                Source(j) => {uf.union(mid(i), mid(j))},
                Target(j) => {uf.union(mid(i), right(j))},
            }
        }
        for i in 1..other.co_domain()+1 {
            match other.link(Target(i)) {
                Target(j) => {uf.union(right(i),right(j))},
                _ => {},
            }
        }

        let mut links : Vec<Link> = Vec::new();
        for i in 1..self.domain()+1 {
            let targ : Site = match uf.set(i).find(|x| *x.1 > self.domain() + self.co_domain()) {
                Some(j) => Target(j.1 - self.domain() - self.co_domain()),
                None => Source(
                    *uf.set(i).find(|x| *x.1 <= self.domain() && i != *x.1).unwrap().1
                ),
            };
            let link = Link::new(Source(i), targ);
            if links.iter().find(|x| **x == link).is_none() {
                links.push(link);
            }
        }
        for i in 1..other.co_domain()+1 {
            if let Some(j) = uf.set(right(i)).find(|x| *x.1 > i + self.domain() + self.co_domain()) {
                links.push(Link::new(Target(i), Target(j.1 - self.domain() - self.co_domain())));
            }
        }
        (
            uf.amount_of_sets() -
                (self.domain() + other.co_domain())/2 - 1,
            TLDiagram::new(links)
        )
    }
}

impl std::ops::BitOr for TLDiagram {
    type Output = TLDiagram;

    fn bitor(mut self, other: TLDiagram) -> TLDiagram {
        let n = self.domain();
        let m = self.co_domain();
        for lnk in other.0.iter() {
            self.0.push(lnk.shift(n,m));
        }
        self.0.sort();
        self
    }
}

impl Tex for TLDiagram {
    fn into_tex(&self) -> String {
        let width = self.domain().max(self.co_domain()) / 2;
        let mut ans = String::from("\\vcenter{\\hbox{\\begin{tikzpicture}[scale=0.3]");
        ans += &format!("\\draw[thin] (0,.5) -- (0, {});",self.domain() as f32 + 0.5);
        ans += &format!("\\draw[thin] ({0},.5) -- ({0}, {1});",width, self.co_domain() as f32 + 0.5);
        for i in 1..self.domain()+1 {
            ans += &format!("\\fill (0,{}) circle(0.1);",i);
        }
        for i in 1..self.co_domain()+1 {
            ans += &format!("\\fill ({0},{1}) circle(0.1);",width, i);
        }
        for lnk in self.0.iter() {
            if let Link(Source(i), Source(j)) = lnk {
                ans += &format!(
                    "\\draw[very thick] (0,{}) edge[out=0, in=0] (0,{});",
                    self.domain()-i+1,
                    self.domain()-j+1);
            }
            if let Link(Target(i), Target(j)) = lnk {
                ans += &format!(
                    "\\draw[very thick] ({0},{1}) edge[out=180, in=180] ({0},{2});",
                    width,
                    self.co_domain()-i+1,
                    self.co_domain()-j+1);
            }
            if let Link(Source(i), Target(j)) = lnk {
                ans += &format!(
                    "\\draw[very thick] (0,{1}) edge[out=0, in=180] ({0},{2});",
                    width,
                    self.domain()-i+1,
                    self.co_domain()-j+1);
            }
        }
        ans += "\\end{tikzpicture} } }";
        ans
    }

    fn is_multiterm(&self) -> bool {
        false
    }
}

impl Serialisable for TLDiagram {
    fn serialise(&self) -> String {
        self.0.serialise()
    }

    fn deserialise(inpt : &str) -> Self {
        TLDiagram::new(
            Vec::deserialise(inpt)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn caps_cups() {
        debug_assert_eq!(TLDiagram::cap(10,3).domain(), 10);
        debug_assert_eq!(TLDiagram::cap(10,3).co_domain(), 8);
        debug_assert_eq!(TLDiagram::cup(10,3).domain(), 8);
        debug_assert_eq!(TLDiagram::cup(10,3).co_domain(), 10);
        debug_assert_eq!(TLDiagram::cup(4,3).inject(), TLDiagram::cup(5,3));
        debug_assert_eq!(TLDiagram::cap(2,1), TLDiagram::new(vec![Link::new(Source(1), Source(2))]));
    }

    #[test]
    fn from_tableaux() {
        debug_assert_eq!(TLDiagram::cap(10,3), TLDiagram::from_tableaux(10, vec![4].into_iter()));
    }

    #[test]
    fn composition() {
        let a = TLDiagram::cap(4,1) * TLDiagram::cup(4,2);
        debug_assert_eq!(a , (0, TLDiagram::new(vec![Link::new(Source(1), Source(2)),Link::new(Source(3), Target(1)),Link::new(Source(4), Target(4)),Link::new(Target(3), Target(2))])));
        let a = TLDiagram::cap(4,1) * TLDiagram::cup(4,1);
        debug_assert_eq!(a , (0, TLDiagram::new(vec![Link::new(Source(1), Source(2)),Link::new(Source(3), Target(3)),Link::new(Source(4), Target(4)),Link::new(Target(1), Target(2))])));
        let a = TLDiagram::cup(4,1) * TLDiagram::cap(4,1);
        debug_assert_eq!(a , (1, TLDiagram::new(vec![Link::new(Source(1), Target(1)),Link::new(Source(2), Target(2))])));
    }

    #[test]
    fn inject() {
        let a = (TLDiagram::cap(4,1) * TLDiagram::cup(4,2)).1;
        debug_assert_eq!(a.inject(), a | TLDiagram::id(1));
    }

    #[test]
    fn any() {
        debug_assert_eq!(TLDiagram::any(5,3), TLDiagram::from_tableaux(5, vec![2].into_iter()));
        debug_assert_eq!(TLDiagram::any(15,5).domain(), 15);
        debug_assert_eq!(TLDiagram::any(15,5).co_domain(), 5);
        debug_assert_eq!(TLDiagram::any(16,12).domain(), 16);
        debug_assert_eq!(TLDiagram::any(16,12).co_domain(), 12);
    }

    #[test]
    fn turn() {
        debug_assert_eq!(TLDiagram::id(3).turn_down(1), TLDiagram::from_tableaux(4, vec![4].into_iter()));
        debug_assert_eq!(TLDiagram::any(14,4).turn_down(3).turn_down(-3), TLDiagram::any(14,4));
        debug_assert_eq!(TLDiagram::id(5).turn_up(-2).turn_down(2), TLDiagram::from_tableauxs(5,vec![4,5].into_iter(), vec![3,4].into_iter()));
    }
}
