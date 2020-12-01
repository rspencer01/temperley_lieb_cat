pub trait GCD {
    fn gcd(&self, other: &Self) -> Self;
}

impl GCD for i128 {
    fn gcd(&self, other: &i128) -> i128 {
        if *self < *other {
            other.gcd(self)
        } else {
            let mut a = *self;
            let mut b = *other;
            while b != 0 {
                let c = a % b;
                a = b;
                b = c;
            }
            a
        }
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
    }
}
