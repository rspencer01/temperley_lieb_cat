//! Operations for dealing with expanding numbers in terms of their digits
//!
//! The representation theory of the Temperley-Lieb category is closely tied
//! to the two torsion parameters, $p$ and $\ell$, of the underlying ring.
//! This module contains various helper functions for dealing with numbers in
//! their $(\ell, p)$-adic expansions.
//!
//! Let $p^{(i)} = \ell p^{i-1}$ where we understand $p^{(0)} = 1$.
//! Then the $(\ell, p)$-adic expansion of $n$ is written
//! $n = [n\_r, \ldots, n\_0]\_{p, \ell}$ where
//! $n = \sum_{i = 0}^{r} n\_ip^{(i)}$ and $0\le n_0 < \ell$ and 
//! $0 \le n_i < p$ for $0< i$.

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

/// Calculates the support (cousins) of $n$ in $(\ell,p)$ digits
///
/// The support of $n$ is the set
/// $$
/// \\big\\{[n\_r, \pm n\_{r-1},\ldots, \pm n\_0]\_{p, \ell} - 1\\big\\}
/// $$
/// where $n + 1 = [n_r, n_{r-1},\ldots, n_0]_{p, \ell}$.
///
/// See also [`support_digits`] for a function to obtain the individual
/// elements by their $(\ell, p)$-digits.
///
/// ### Example
///
/// ```rust
/// # use temperley_lieb_cat::digits::support;
/// assert_eq!(
///     support(10, 2, 3),
///     [10, 8, 2, 0]
/// );
/// ```
pub fn support(n: u128, l: u128, p: u128) -> Vec<u128> {
    support_digits(n, l, p)
        .iter()
        .map(|x| anti_l_p_adic(x, l, p) as u128 - 1)
        .collect::<Vec<_>>()
}

/// Calculates the digits of the support (cousins) of $n$ in $(\ell,p)$ digits
///
/// The support of $n$ is the set
/// $$
/// \\big\\{[n\_r, \pm n\_{r-1},\ldots, \pm n\_0]\_{p, \ell} - 1\\big\\}
/// $$
/// where $n + 1 = [n_r, n_{r-1},\ldots, n_0]_{p, \ell}$.  This function returns
/// a list of the $(\ell, p)$ digits of the support elements.
///
/// See also [`support`] for a function to obtain the numerical elements.
///
/// ### Note
/// The digits are presented from least to most significant digit.  That is, the
/// first element is $\pm n\_0$, the second $\pm n\_1$ and so on.
///
/// ### Example
///
/// ```rust
/// # use temperley_lieb_cat::digits::support_digits;
/// assert_eq!(
///     support_digits(9, 2, 3),
///     [[0,2,1], [0,-2,1]]
/// );
/// ```
pub fn support_digits(n: u128, l: u128, p: u128) -> Vec<Vec<i128>> {
    let mut ans = vec![];
    let digs = l_p_adic(n + 1, l, p);
    let indices: Vec<usize> = digs
        .iter()
        .enumerate()
        .take(digs.len() - 1)
        .filter(|(_, d)| **d != 0)
        .map(|(i, _)| i)
        .collect();

    for t in 0..(1 << indices.len()) {
        let mut these_digs = digs.clone();
        for i in 0..indices.len() {
            if t & (1 << i) != 0 {
                these_digs[indices[i]] *= -1;
            }
        }
        ans.push(these_digs);
    }
    ans
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

    #[test]
    fn test_support() {
        assert_eq!(
            support(100, 2, 3),
            [100, 98, 92, 90, 88, 86, 80, 78, 28, 26, 20, 18, 16, 14, 8, 6]
        );
        for p in [2, 3, 5] {
            for l in 2..10 {
                assert_eq!(support(0, l, p), [0]);
                assert_eq!(support(l - 1, l, p), [l - 1]);
                assert_eq!(
                    support((p - 1) * l * p * p - 1, l, p),
                    [(p - 1) * l * p * p - 1]
                );
                assert_eq!(
                    support((p - 1) * l * p * p, l, p),
                    [(p - 1) * l * p * p, (p - 1) * l * p * p - 2]
                );
            }
        }
    }
}
