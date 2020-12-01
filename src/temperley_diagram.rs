extern crate partitions;

#[derive(Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Site {
    Source(usize),
    Target(usize),
}

impl std::fmt::Debug for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl std::fmt::Display for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Source(i) => write!(f, "{}·", i),
            Target(i) => write!(f, "·{}", i),
        }
    }
}

impl Site {
    fn involute(self) -> Site {
        match self {
            Source(i) => Target(i),
            Target(i) => Source(i),
        }
    }
}

use Site::*;

#[derive(Copy, Clone, Hash, Ord, PartialOrd, Eq, PartialEq)]
pub struct Link(pub Site, pub Site);

impl std::fmt::Debug for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl std::fmt::Display for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f,"{}—{}", self.0, self.1)
    }
}

impl Link {
    pub fn new(site_a : Site, site_b : Site) -> Link {
        match (site_a, site_b) {
            (Source(i), Source(j)) => Link(Source(usize::min(i, j)), Source(usize::max(i, j))),
            (Target(i), Target(j)) => Link(Target(usize::min(i, j)), Target(usize::max(i, j))),
            (Source(i), Target(j)) => Link(Source(i), Target(j)),
            (Target(i), Source(j)) => Link(Source(j), Target(i)),
        }
    }

    fn involute(self) -> Link {
        Link::new(self.0.involute(), self.1.involute())
    }
}

#[derive(Clone)]
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

    pub fn from_tableaux(n : usize, tab : Vec<usize>) -> TLDiagram {
        let mut stack = Vec::new();
        let mut links = Vec::new();
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

    pub fn from_tableauxs(n : usize, tab1 : Vec<usize>, tab2 : Vec<usize>) -> TLDiagram {
        (TLDiagram::from_tableaux(n, tab1) * TLDiagram::from_tableaux(n, tab2).involute()).1
    }

    pub fn u(n :usize, i : usize) -> TLDiagram {
        TLDiagram::from_tableauxs(n, vec![i+1], vec![i+1])
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

    pub fn involute(self) -> TLDiagram {
        TLDiagram::new(
            self.0.into_iter().map(|x| {
                x.involute()
            }).collect()
        )
    }

    pub fn cap(n: usize, i: usize) -> TLDiagram {
        let mut ans : Vec<Link> = (1..(i))
                    .map(|j| Link::new(Source(j), Target(j)))
                    .chain((i+2..n+1)
                        .map(|j| Link::new(Source(j), Target(j-2)))
                    )
                    .collect();
        ans.push(Link::new(Source(i), Source(i+1)));
        TLDiagram::new(ans)
    }

    pub fn cup(n: usize, i: usize) -> TLDiagram {
        TLDiagram::cap(n,i).involute()
    }

    pub fn identity(n: usize) -> TLDiagram {
        TLDiagram::new((1..n+1).map(|j| Link::new(Source(j), Target(j))).collect())
    }

    pub fn inject(&self) -> TLDiagram {
        let mut links = self.0.clone();
        links.push(Link::new(Source(self.domain()+1), Target(self.co_domain()+1)));
        TLDiagram::new(links)
    }

    pub fn link(&self, site: Site) -> Site {
        for Link(a, b) in self.0.iter() {
            if *a == site { return *b }
            if *b == site { return *a }
        }
        panic!("Site {:?} not found",site)
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
        assert_eq!(self.co_domain(), other.domain());
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
        debug_assert_eq!(TLDiagram::cap(10,3), TLDiagram::from_tableaux(10, vec![4]));
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
}
