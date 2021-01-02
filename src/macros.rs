#[macro_export]
macro_rules! val_op_by_ref {
    ($op:tt, $op_fn:tt, $t:ty, <$N:ident $(: $b0:ident $(+$b:ident$(<$c:ident=$d:ty>)?)* )?>) => (
        impl<$N $(: $b0 $(+$b$(<$c=$d>)?)* )?> $op<$t> for $t {
            type Output = $t;

            fn $op_fn(self, other: $t) -> $t {
                (&self).$op_fn(&other)
            }
        }

        impl<$N $(: $b0 $(+$b$(<$c=$d>)?)* )?> $op<&$t> for $t {
            type Output = $t;

            fn $op_fn(self, other: &$t) -> $t {
                (&self).$op_fn(other)
            }
        }

        impl<$N $(: $b0 $(+$b$(<$c=$d>)?)* )?> $op<$t> for &$t {
            type Output = $t;

            fn $op_fn(self, other: $t) -> $t {
                self.$op_fn(&other)
            }
        }
    )
}

