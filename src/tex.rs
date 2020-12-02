use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{Command, Stdio};

pub trait Tex {
    fn into_tex(&self) -> String;

    fn is_multiterm(&self) -> bool;

    fn render_tex(&self) {
        let document = format!(
    "\\documentclass[crop=true]{{standalone}}
    \\standaloneconfig{{border=1em}}
    \\usepackage{{amsmath}}
    \\usepackage{{amssymb}}
    \\usepackage{{tikz}}
    \\begin{{document}}

    $
    {}
    $

    \\end{{document}}", self.into_tex());
        let mut command = Command::new("/home/robert/bin/pdflatex")
            .stdin(Stdio::piped())
            .arg("-output-directory=/tmp")
            .arg("-shell-escape")
            .spawn()
            .unwrap();
        command.stdin
            .as_mut()
            .unwrap()
            .write_all(&document.into_bytes())
            .expect("Can't write tex to pdflatex");
        command.wait().expect("Cannot await task");
    }
}


impl Tex for i128 {
    fn into_tex(&self) -> String {
        format!("{}",self)
    }

    fn is_multiterm(&self) -> bool {
        false
    }
}
