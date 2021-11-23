/// Types that can calculate greatest common factors.
pub trait GCD {
    /// Find a greatest common factor of self and other
    fn gcd(&self, other: &Self) -> Self;
}

/// Types that can calculate common factors.
pub trait PartialGCD {
    /// Find a (large) common factor of self and other
    fn partial_gcd(&self, other: &Self) -> Self;

    /// Is the number ``small''
    fn is_small(&self) -> bool;
}

impl<T: GCD> PartialGCD for T {
    #[inline(always)]
    fn partial_gcd(&self, other: &Self) -> Self {
        self.gcd(other)
    }

    #[inline(always)]
    fn is_small(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gcd() {
        debug_assert_eq!(15.gcd(&5), 5);
        debug_assert_eq!(5.gcd(&15), 5);
        debug_assert_eq!(5.gcd(&19), 1);
        debug_assert_eq!(19.gcd(&19), 19);
        debug_assert_eq!(0.gcd(&17), 17);
        debug_assert_eq!(1.gcd(&-1), 1);
        debug_assert_eq!((-15).gcd(&5), 5);
    }
}
