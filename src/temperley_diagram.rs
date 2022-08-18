extern crate itertools;
extern crate partitions;

use crate::serial::Serialisable;
use crate::tex::Tex;
use itertools::Itertools;
use std::fmt::Write as _;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TLDiagram {
    domain: usize,
    co_domain: usize,
    left_tab: u64,
    right_tab: u64,
}

impl std::fmt::Debug for TLDiagram {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "TL {{ {} ↦ {} {:b} {:b} }}",
            self.domain, self.co_domain, self.left_tab, self.right_tab
        )
    }
}

impl std::fmt::Display for TLDiagram {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl TLDiagram {
    pub fn new(domain: usize, co_domain: usize, left_tab: u64, right_tab: u64) -> TLDiagram {
        assert!(
            domain as u32 - 2 * left_tab.count_ones()
                == co_domain as u32 - 2 * right_tab.count_ones()
        );
        TLDiagram {
            domain,
            co_domain,
            left_tab,
            right_tab,
        }
    }

    pub fn u(n: usize, i: usize) -> TLDiagram {
        TLDiagram::new(n, n, 1 << (i + 1), 1 << (i + 1))
    }

    // IS this correct???
    pub fn big_u(n: usize, i: usize, j: usize) -> TLDiagram {
        TLDiagram::new(
            n,
            n,
            (1 << (i + j + 1)) - (1 << i),
            (1 << (i + j + 1)) - (1 << i),
        )
    }

    pub fn domain(&self) -> usize {
        self.domain
    }

    pub fn co_domain(&self) -> usize {
        self.co_domain
    }

    pub fn propagation(&self) -> usize {
        self.domain - 2 * self.left_tab.count_ones() as usize
    }

    pub fn involute(self) -> TLDiagram {
        TLDiagram::new(self.co_domain, self.domain, self.right_tab, self.left_tab)
    }

    pub fn flip(&self) -> TLDiagram {
        let mut new_left_tab = 0;
        let mut stack = Vec::new();
        for i in 1..self.domain + 1 {
            if (self.left_tab & (1 << i)) != 0 {
                new_left_tab |= 1 << (self.domain + 1 - stack.pop().unwrap()) as u64;
            } else {
                stack.push(i);
            }
        }
        let mut new_right_tab = 0;
        let mut stack = Vec::new();
        for i in 1..self.co_domain + 1 {
            if (self.right_tab & (1 << i)) != 0 {
                new_right_tab |= 1 << (self.co_domain + 1 - stack.pop().unwrap()) as u64;
            } else {
                stack.push(i);
            }
        }
        TLDiagram::new(self.domain, self.co_domain, new_left_tab, new_right_tab)
    }

    pub fn cap(n: usize, i: usize) -> TLDiagram {
        TLDiagram::new(n, n - 2, 1 << (i + 1), 0)
    }

    pub fn cup(n: usize, i: usize) -> TLDiagram {
        TLDiagram::cap(n, i).involute()
    }

    pub fn big_cap(n: usize, i: usize, j: usize) -> TLDiagram {
        TLDiagram::new(n, n - 2 * j, (1 << (i + j + 1)) - (1 << (i + 1)), 0)
    }

    pub fn id(n: usize) -> TLDiagram {
        TLDiagram::new(n, n, 0, 0)
    }

    pub fn inject(&self) -> TLDiagram {
        TLDiagram::new(
            self.domain + 1,
            self.co_domain + 1,
            self.left_tab,
            self.right_tab,
        )
    }

    pub fn any(n: usize, m: usize) -> TLDiagram {
        if n < m {
            TLDiagram::any(m, n).involute()
        } else {
            let mut v = 0;
            for i in 0..(n - m) / 2 {
                v |= 1 << (i * 2 + 2);
            }
            TLDiagram::new(n, m, v, 0)
        }
    }
    pub fn remove_cap(&self, i: usize) -> TLDiagram {
        TLDiagram {
            left_tab: (self.left_tab & ((1 << i) - 1))
                | ((self.left_tab & !((1 << (i + 2)) - 1)) >> 2),
            domain: self.domain - 2,
            ..*self
        }
    }

    pub fn all(n: usize, m: usize) -> Vec<TLDiagram> {
        fn ok(tab: &Vec<usize>) -> bool {
            for (c, i) in tab.iter().enumerate() {
                if *i < 2 * (c + 1) {
                    return false;
                }
            }
            true
        }
        if n < m {
            TLDiagram::all(m, n)
                .into_iter()
                .map(|d| d.involute())
                .collect()
        } else {
            let mut ans = Vec::new();
            for rm in 0..(m / 2 + 1) {
                let prop = m - 2 * rm;
                let rn = (n - prop) / 2;
                for lt in (1..(n + 1)).combinations(rn).filter(ok) {
                    for rt in (1..(m + 1)).combinations(rm).filter(ok) {
                        ans.push(TLDiagram::new(
                            n,
                            m,
                            lt.iter().fold(0, |l, n| l | (1 << n)),
                            rt.iter().fold(0, |l, n| l | (1 << n)),
                        ));
                    }
                }
            }
            ans
        }
    }

    pub fn simple_links(&self) -> Vec<usize> {
        (1..self.domain() + 1)
            .filter(|i| (self.left_tab & (1 << i) != 0) & (self.left_tab & (1 << (i - 1)) == 0))
            .map(|i| i - 1)
            .collect()
    }

    pub fn turn_down(&self, i: isize) -> TLDiagram {
        if i == 0 {
            return *self;
        }
        if i > 0 {
            let i = i as usize;
            assert!(self.co_domain >= i, "Cannot turn element");
            let mut left_tab = self.left_tab;
            let mut right_tab = self.right_tab;
            for j in 0..i {
                if right_tab & (1 << (self.co_domain - j)) != 0 {
                    right_tab ^= 1 << (self.co_domain - j);
                } else {
                    left_tab |= 1 << (self.domain + 1 + j);
                }
            }
            TLDiagram::new(self.domain + i, self.co_domain - i, left_tab, right_tab)
        } else {
            self.involute().turn_down(-i).involute()
        }
    }

    pub fn turn_up(&self, i: isize) -> TLDiagram {
        self.flip().turn_down(i).flip()
    }

    /// Rotate the diagram anticlockwise
    pub fn rotate(&self, i: isize) -> TLDiagram {
        match i {
            0 => *self,
            i if i > 0 => self.turn_up(1).turn_down(-1).rotate(i - 1),
            _ => self.turn_up(1).turn_down(-1).rotate(i + 1),
        }
    }
}

impl std::ops::Mul for TLDiagram {
    type Output = (usize, TLDiagram);

    fn mul(self, other: TLDiagram) -> (usize, TLDiagram) {
        &self * &other
    }
}

impl std::ops::Mul for &TLDiagram {
    type Output = (usize, TLDiagram);

    fn mul(self, other: &TLDiagram) -> (usize, TLDiagram) {
        assert_eq!(
            self.co_domain, other.domain,
            "Cannot multiply diagrams if their (co)domains don't match"
        );
        // Inlcusive
        // Sources 1   .. a
        // Middle  a+1 .. a+b
        // Target  a+b .. a+b+c
        let left = |x: usize| x;
        let mid = |x: usize| x + self.domain;
        let right = |x: usize| x + self.domain + self.co_domain;
        let mut uf: partitions::PartitionVec<usize> = partitions::PartitionVec::new();
        for i in 0..self.domain + self.co_domain + other.co_domain + 1 {
            uf.push(i);
        }

        let mut left_stack = Vec::new();
        for i in 1..self.domain() + 1 {
            if self.left_tab & (1 << i) != 0 {
                uf.union(left(left_stack.pop().unwrap()), left(i));
            } else {
                left_stack.push(i);
            }
        }
        let mut right_stack = Vec::new();
        for i in 1..self.co_domain + 1 {
            if self.right_tab & (1 << i) != 0 {
                uf.union(mid(right_stack.pop().unwrap()), mid(i));
            } else {
                right_stack.push(i);
            }
        }
        assert!(left_stack.len() == right_stack.len());
        for i in 0..left_stack.len() {
            uf.union(left(left_stack[i]), mid(right_stack[i]));
        }

        let mut left_stack = Vec::new();
        for i in 1..other.domain() + 1 {
            if other.left_tab & (1 << i) != 0 {
                uf.union(mid(left_stack.pop().unwrap()), mid(i));
            } else {
                left_stack.push(i);
            }
        }
        let mut right_stack = Vec::new();
        for i in 1..other.co_domain + 1 {
            if other.right_tab & (1 << i) != 0 {
                uf.union(right(right_stack.pop().unwrap()), right(i));
            } else {
                right_stack.push(i);
            }
        }
        assert!(left_stack.len() == right_stack.len());
        for i in 0..left_stack.len() {
            uf.union(mid(left_stack[i]), right(right_stack[i]));
        }

        let mut left_tab = 0;
        for i in 1..self.domain + 1 {
            if let Some(d) = uf.set(left(i)).find(|x| i < *x.1 && *x.1 <= self.domain()) {
                left_tab |= 1 << (*d.1);
            }
        }
        let mut right_tab = 0;
        for i in 1..other.co_domain + 1 {
            if let Some(d) = uf
                .set(right(i))
                .find(|x| self.domain + self.co_domain + i < *x.1)
            {
                right_tab |= 1 << (d.1 - self.domain - self.co_domain);
            }
        }
        (
            uf.amount_of_sets() - (self.domain() + other.co_domain()) / 2 - 1,
            TLDiagram::new(self.domain, other.co_domain, left_tab, right_tab),
        )
    }
}

impl std::ops::BitOr for TLDiagram {
    type Output = TLDiagram;

    fn bitor(mut self, other: TLDiagram) -> TLDiagram {
        let n = self.domain;
        let m = self.co_domain;
        self.left_tab |= other.left_tab << n;
        self.right_tab |= other.right_tab << m;
        self.domain += other.domain;
        self.co_domain += other.co_domain;
        self
    }
}

impl std::ops::BitOr<&TLDiagram> for TLDiagram {
    type Output = TLDiagram;

    fn bitor(mut self, other: &TLDiagram) -> TLDiagram {
        let n = self.domain;
        let m = self.co_domain;
        self.left_tab |= other.left_tab << n;
        self.right_tab |= other.right_tab << m;
        self.domain += other.domain;
        self.co_domain += other.co_domain;
        self
    }
}

impl Tex for TLDiagram {
    fn into_tex(&self) -> String {
        let width = self.domain().max(self.co_domain()) / 2;
        let mut ans = String::from("\\vcenter{\\hbox{\\begin{tikzpicture}[scale=0.3]");
        let _ = write!(
            ans,
            "\\draw[thin] (0,.5) -- (0, {});",
            self.domain() as f32 + 0.5
        );
        let _ = write!(
            ans,
            "\\draw[thin] ({0},.5) -- ({0}, {1});",
            width,
            self.co_domain() as f32 + 0.5
        );
        for i in 1..self.domain() + 1 {
            let _ = write!(ans, "\\fill (0,{}) circle(0.1);", i);
        }
        for i in 1..self.co_domain() + 1 {
            let _ = write!(ans, "\\fill ({0},{1}) circle(0.1);", width, i);
        }
        let mut left_stack = Vec::new();
        let mut right_stack = Vec::new();
        for i in 1..self.domain + 1 {
            if self.left_tab & (1 << i) != 0 {
                let _ = write!(
                    ans,
                    "\\draw[very thick] (0,{}) edge[out=0, in=0] (0,{});",
                    self.domain - left_stack.pop().unwrap() + 1,
                    self.domain - i + 1
                );
            } else {
                left_stack.push(i);
            }
        }
        for i in 1..self.co_domain + 1 {
            if self.right_tab & (1 << i) != 0 {
                let _ = write!(
                    ans,
                    "\\draw[very thick] ({0},{1}) edge[out=180, in=180] ({0},{2});",
                    width,
                    self.co_domain - right_stack.pop().unwrap() + 1,
                    self.co_domain - i + 1
                );
            } else {
                right_stack.push(i)
            }
        }
        for i in 0..left_stack.len() {
            let _ = write!(
                ans,
                "\\draw[very thick] (0,{1}) edge[out=0, in=180] ({0},{2});",
                width,
                self.domain() - left_stack[i] + 1,
                self.co_domain() - right_stack[i] + 1
            );
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
        format!(
            "<{}.{}.{}.{}>",
            self.domain.serialise(),
            self.co_domain.serialise(),
            self.left_tab.serialise(),
            self.right_tab.serialise()
        )
    }

    fn deserialise(inpt: &str) -> Self {
        assert!(inpt.starts_with('<'));
        assert!(inpt.ends_with('>'));
        let mut items = inpt[1..inpt.len() - 1].split('.');
        TLDiagram {
            domain: usize::deserialise(items.next().unwrap()),
            co_domain: usize::deserialise(items.next().unwrap()),
            left_tab: u64::deserialise(items.next().unwrap()),
            right_tab: u64::deserialise(items.next().unwrap()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn caps_cups() {
        assert_eq!(TLDiagram::cap(10, 3).domain(), 10);
        assert_eq!(TLDiagram::cap(10, 3).co_domain(), 8);
        assert_eq!(TLDiagram::cup(10, 3).domain(), 8);
        assert_eq!(TLDiagram::cup(10, 3).co_domain(), 10);
        assert_eq!(TLDiagram::cup(4, 3).inject(), TLDiagram::cup(5, 3));
        assert_eq!(TLDiagram::cap(2, 1), TLDiagram::new(2, 0, 0b100, 0));
    }
    #[test]
    fn cap_from_tableaux() {
        assert_eq!(
            TLDiagram::cap(10, 3),
            TLDiagram::new(10, 8, 0b00000010000, 0)
        );
    }

    #[test]
    fn composition() {
        let a = TLDiagram::cap(4, 1) * TLDiagram::cup(4, 2);
        assert_eq!(a, (0, TLDiagram::new(4, 4, 0b00100, 0b01000)));
        let a = TLDiagram::cap(4, 1) * TLDiagram::cup(4, 1);
        assert_eq!(a, (0, TLDiagram::new(4, 4, 0b00100, 0b00100)));
        let a = TLDiagram::cup(4, 1) * TLDiagram::cap(4, 1);
        debug_assert_eq!(a, (1, TLDiagram::id(2)));
    }

    #[test]
    fn inject() {
        let a = TLDiagram::u(4, 2);
        assert_eq!(a.inject(), a | TLDiagram::id(1));
    }

    #[test]
    fn any() {
        assert_eq!(TLDiagram::any(15, 5).domain(), 15);
        assert_eq!(TLDiagram::any(15, 5).co_domain(), 5);
        assert_eq!(TLDiagram::any(16, 12).domain(), 16);
        assert_eq!(TLDiagram::any(16, 12).co_domain(), 12);
    }

    #[test]
    fn turn() {
        assert_eq!(TLDiagram::id(3).turn_down(1), TLDiagram::cap(4, 3));
        assert_eq!(
            TLDiagram::id(3).turn_down(2),
            TLDiagram::new(5, 1, 0b110000, 0)
        );
        assert_eq!(
            TLDiagram::any(14, 4).turn_down(3).turn_down(-3),
            TLDiagram::any(14, 4)
        );
        assert_eq!(
            TLDiagram::id(5).turn_up(-2).turn_down(2),
            TLDiagram::new(5, 5, 0b110000, 0b011000)
        );
    }

    #[test]
    fn all() {
        println!("{:?}", TLDiagram::all(3, 1));
        use std::collections::HashSet;
        assert_eq!(TLDiagram::all(3, 1).iter().collect::<HashSet<_>>().len(), 2);
        assert_eq!(
            TLDiagram::all(9, 5).iter().collect::<HashSet<_>>().len(),
            42 * 5 + 48 * 4 + 27 * 1
        );
    }

    #[test]
    fn remove_cap() {
        assert_eq!(TLDiagram::cap(8, 3).remove_cap(3), TLDiagram::id(6));
        assert_eq!(
            TLDiagram::new(8, 6, 0b01100100, 0b001100).remove_cap(1),
            TLDiagram::new(6, 6, 0b00011000, 0b001100)
        );
        assert_eq!(
            TLDiagram::new(8, 6, 0b01100100, 0b001100).remove_cap(4),
            TLDiagram::new(6, 6, 0b00010100, 0b001100)
        );
    }
}
