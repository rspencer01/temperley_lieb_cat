//! Traits for writing out TeX files

use std::io::Write;
use std::process::{Command, Stdio};
extern crate which;
use which::which;

/// A trait to render the type into (La)TeX
///
/// Often it is useful to view the results of computations in a visual
/// format.  Here this is done by rendering to TeX and then to PDF.  See [Tex::render_tex] for details.
///
/// To see which types can be rendered to TeX, check which types
/// implement this trait.
/// To add this functionality to an existing type, you will need to provide [Tex::into_tex] and
/// [Tex::is_multiterm].
///
/// The document framework that the strucutre is rendered into is
/// ```tex
/// \documentclass[crop=true]{standalone}
/// \standaloneconfig{border=1em}
/// \usepackage{amsmath}
/// \usepackage{amssymb}
/// \usepackage{tikz}
/// \begin{document}
///
/// $
/// {}
/// $
///
/// \end{document}
/// ```
///
/// ## Example
/// ```rust
/// # use temperley_lieb_cat::Tex;
/// # use temperley_lieb_cat::Polynomial;
/// Polynomial::<i128>::new(&[1,2,0,0,1]).render_tex();
/// ```
pub trait Tex {
    /// Write out the (math-mode) TeX to represent this object
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

impl<T: Tex> Tex for [T] {
    fn into_tex(&self) -> String {
        let mut ans = String::new();
        ans += "\\left(";
        let mut first = true;
        for i in self {
            if !first {
                ans += ", ";
            } else {
                first = false;
            }
            ans += &i.into_tex();
        }
        ans += "\\right)";
        ans
    }

    fn is_multiterm(&self) -> bool {
        false
    }
}
