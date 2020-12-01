/// A site in a Temperley-Lieb diagram.
///
/// May reside on the left/bottom (source) or the top/right (target)
/// Sites are stored context-less, which is to say they are agnostic
/// to the number of other sites in the diagram.
#[derive(Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Site {
    /// Left / bottom nodes
    Source(usize),
    /// Right / top nodes
    Target(usize),
}

use Site::{Source, Target};

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
    pub fn involute(self) -> Site {
        match self {
            Source(i) => Target(i),
            Target(i) => Source(i),
        }
    }
}
