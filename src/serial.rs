use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

/// Traits that can be serialised or deserialised to strings
pub trait Serialisable {
    fn serialise(&self) -> String;

    fn deserialise(inpt: &str) -> Self;

    fn to_file(&self, path: &Path) {
        File::create(path)
            .expect("Cannot open file for writing")
            .write(self.serialise().as_bytes())
            .expect("Could not write object");
    }

    fn from_file(path: &Path) -> Self
    where
        Self: Sized,
    {
        let mut buff = String::new();
        File::open(path)
            .expect("Cannot open file for reading")
            .read_to_string(&mut buff)
            .expect("Could not read file contents");
        Self::deserialise(buff.as_str())
    }
}

impl Serialisable for i128 {
    fn serialise(&self) -> String {
        format!("{}", self)
    }

    fn deserialise(inpt: &str) -> Self {
        inpt.parse::<i128>().expect("Could not parse int")
    }
}

impl Serialisable for usize {
    fn serialise(&self) -> String {
        format!("{}", self)
    }

    fn deserialise(inpt: &str) -> Self {
        inpt.parse::<usize>().expect("Could not parse int")
    }
}

impl Serialisable for u64 {
    fn serialise(&self) -> String {
        format!("{}", self)
    }

    fn deserialise(inpt: &str) -> Self {
        inpt.parse::<u64>().expect("Could not parse int")
    }
}

impl<T> Serialisable for Vec<T>
where
    T: Serialisable,
{
    fn serialise(&self) -> String {
        self.iter().fold(String::from("["), |acc, v| {
            acc + v.serialise().as_str() + ","
        }) + "]"
    }

    fn deserialise(inpt: &str) -> Self {
        assert!(
            inpt.chars().nth(0) == Some('['),
            "No open bracket in vec parsing"
        );
        assert!(
            inpt.chars().nth(inpt.len() - 1) == Some(']'),
            "No close bracket in vec parsing"
        );
        let mut ans = Vec::new();
        for item in inpt[1..inpt.len() - 1].split(",") {
            if item.len() == 0 {
                break;
            }
            ans.push(T::deserialise(item));
        }
        ans
    }
}

mod test {
    #[cfg(test)]
    use crate::serial::Serialisable;

    #[test]
    fn vec_serial() {
        assert_eq!(Vec::<i128>::new().serialise(), "[]");
        assert!(Vec::<i128>::deserialise("[]").is_empty());
        assert_eq!(vec![1i128, 2, 3].serialise(), "[1,2,3,]");
        assert_eq!(Vec::<i128>::deserialise("[1,2,3,]"), vec![1, 2, 3]);
    }
}
