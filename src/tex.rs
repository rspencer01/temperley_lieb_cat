//! Traits for writing out TeX files

use std::io::Write;
use std::process::{Command, Stdio};
extern crate which;
use which::which;

/// A trait to render the type into (La)TeX
pub trait Tex {
    fn into_tex(&self) -> String;

    /// Indicate if this element should be surrounded by braces to eliminate ambiguity
    ///
    /// This should be true for items separated by `+`, such as `x + y`, but false
    /// if there is no ambiguity.  The function `Tex::into_tex` should use this
    /// to determine if sub-expressions need to be places in parenthesis.
    fn is_multiterm(&self) -> bool;

    /// Renders this object to a PDF file for displaying
    ///
    /// The output is written to the file `texput.pdf` in the current directory.
    ///
    /// This can be very useful if your PDF viewer automatically updates its display
    /// when the underlying file changes.
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

    \\end{{document}}",
            self.into_tex()
        );
        let mut command = Command::new(which("pdflatex").expect("Could not find latex compiler"))
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .arg("-shell-escape")
            .spawn()
            .unwrap();
        command
            .stdin
            .as_mut()
            .unwrap()
            .write_all(&document.into_bytes())
            .expect("Can't write tex to pdflatex");
        command.wait().expect("Cannot await task");
    }
}

impl Tex for i128 {
    fn into_tex(&self) -> String {
        format!("{}", self)
    }

    fn is_multiterm(&self) -> bool {
        false
    }
}
