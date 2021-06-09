//! Operations for dealing with expanding numbers in terms of their digits

/// Turn a list of digits mod `p` into a number
///
/// The digits are not constrained to the range `[0, p)` but
/// can be negative or larger than `p-1` if required.
pub fn anti_p_adic(digs: &[i128], p: u128) -> i128 {
    let mut ans = 0;
    let mut exp = 1;
    for dig in digs {
        ans += exp * dig;
        exp *= p as i128;
    }
    ans
}

/// Express a number base `p`
///
/// All digits will be within the range `[0, p)`.
pub fn p_adic(n: u128, p: u128) -> Vec<i128> {
    if n == 0 {
        return vec![0];
    }
    let dig = (n % p) as i128;
    let n = n / p;
    if n > 0 {
        let mut ans = p_adic(n, p);
        ans.insert(0, dig);
        ans
    } else {
        vec![dig]
    }
}

/// Turn a list of digits mod p into a number
///
/// The digits are not constrained to and can be arbitrary integers
pub fn anti_l_p_adic(digs: &[i128], l: u128, p: u128) -> i128 {
    if digs.len() == 0 {
        0
    } else {
        digs[0] + l as i128 * anti_p_adic(&digs[1..], p)
    }
}

/// Express a number base (`l`, `p`)
///
/// The first digit will be in the range `[0,l)` and all subsequents in the range `[0, p)`.
pub fn l_p_adic(n: u128, l: u128, p: u128) -> Vec<i128> {
    let dig = (n % l) as i128;
    let n = n / l;
    if n > 0 {
        let mut ans = p_adic(n, p);
        ans.insert(0, dig);
        ans
    } else {
        vec![dig]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn p_digits() {
        assert_eq!(p_adic(42, 2), [0, 1, 0, 1, 0, 1]);
        assert_eq!(p_adic(43, 2), [1, 1, 0, 1, 0, 1]);
        assert_eq!(p_adic(0, 2), [0]);
        assert_eq!(p_adic(1, 2), [1]);
        assert_eq!(p_adic(8, 2), [0, 0, 0, 1]);

        assert_eq!(p_adic(501, 10), [1, 0, 5]);

        assert_eq!(p_adic(100, 3), [1, 0, 2, 0, 1]);
        assert_eq!(p_adic(0, 3), [0]);
    }

    #[test]
    fn anti_p_digits() {
        assert_eq!(anti_p_adic(&[0, 1, 0, 1, 0, 1], 2), 42);
        assert_eq!(anti_p_adic(&[1, 1, 0, 1, 0, 1], 2), 43);
        assert_eq!(anti_p_adic(&[0], 2), 0);
        assert_eq!(anti_p_adic(&[1], 2), 1);
        assert_eq!(anti_p_adic(&[0, 0, 0, 1], 2), 8);

        assert_eq!(anti_p_adic(&[1, 0, 5], 10), 501);

        assert_eq!(anti_p_adic(&[1, 0, 2, 0, 1], 3), 100);
        assert_eq!(anti_p_adic(&[0], 3), 0);
    }

    #[test]
    fn l_p_digits() {
        assert_eq!(l_p_adic(42, 5, 2), [2, 0, 0, 0, 1]);
        assert_eq!(l_p_adic(43, 5, 2), [3, 0, 0, 0, 1]);
        assert_eq!(l_p_adic(0, 5, 2), [0]);
        assert_eq!(l_p_adic(1, 5, 2), [1]);
        assert_eq!(l_p_adic(10, 5, 2), [0, 0, 1]);

        assert_eq!(l_p_adic(501, 2, 10), [1, 0, 5, 2]);
    }

    #[test]
    fn anti_l_p_digits() {
        assert_eq!(anti_l_p_adic(&[2, 0, 0, 0, 1], 5, 2), 42);
        assert_eq!(anti_l_p_adic(&[3, 0, 0, 0, 1], 5, 2), 43);
        assert_eq!(anti_l_p_adic(&[0], 5, 2), 0);
        assert_eq!(anti_l_p_adic(&[1], 5, 2), 1);
        assert_eq!(anti_l_p_adic(&[0, 0, 1], 5, 2), 10);

        assert_eq!(anti_l_p_adic(&[1, 0, 5, 2], 2, 10), 501);
    }
}
