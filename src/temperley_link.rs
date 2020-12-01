use crate::temperley_site::{Site, Site::*};

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

    pub fn involute(self) -> Link {
        Link::new(self.0.involute(), self.1.involute())
    }

    pub fn shift(self, n : usize, m : usize) -> Link {
        match self {
            Link(Source(i), Source(j)) => {
                Link(Source(i+n), Source(j+n))
            },
            Link(Source(i), Target(j)) => {
                Link(Source(i+n), Target(j+m))
            },
            Link(Target(i), Target(j)) => {
                Link(Target(i+m), Target(j+m))
            },
            Link(Target(_), Source(_)) => {
                panic!("Link should not have a target before a source")
            },
        }
    }
}

