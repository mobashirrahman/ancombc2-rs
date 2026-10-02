//! Formula parsing, matching the subset of R's formula language that
//! `ancombc2()`'s `fix_formula` argument actually uses.
//!
//! Supported: an optional leading `~`, an intercept term, `+`-separated terms,
//! `name:other` and `name*other` interactions (which R expands to
//! `name + other + name:other`), numeric literals, and the `-` removal operator.
//! Not supported, and rejected rather than ignored: `|`, `/`, `^`, `.`, `I()`,
//! functions, and transformation syntax like `log(x)`. Silently dropping an
//! unrecognised operator would build a design that does not match the oracle's,
//! which is the one failure mode the parity contract cannot detect.

use std::fmt;

use super::IoError;

/// One term of a formula.
#[derive(Debug, Clone, PartialEq)]
pub enum Term {
    Intercept,
    Variable(String),
    Interaction(String, String),
}

impl Term {
    /// The column names this term contributes, given the metadata's variables.
    pub fn label(&self) -> String {
        match self {
            Term::Intercept => "1".into(),
            Term::Variable(v) => v.clone(),
            Term::Interaction(a, b) => format!("{a}:{b}"),
        }
    }
}

/// A parsed formula: its terms, and whether the intercept was requested.
#[derive(Debug, Clone, PartialEq)]
pub struct Formula {
    pub terms: Vec<Term>,
    pub intercept: bool,
}

impl Formula {
    /// The variable names the formula references, in first-appearance order.
    pub fn variables(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for t in &self.terms {
            for v in match t {
                Term::Intercept => continue,
                Term::Variable(v) => std::slice::from_ref(v),
                Term::Interaction(a, b) => {
                    let pair = [a.clone(), b.clone()];
                    for v in pair {
                        if !out.contains(&v) {
                            out.push(v);
                        }
                    }
                    continue;
                }
            } {
                if !out.contains(v) {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    /// The formula as R would print it, for the result metadata.
    pub fn to_r_string(&self) -> String {
        let mut parts = Vec::new();
        if self.intercept {
            parts.push("1".to_string());
        }
        for t in &self.terms {
            parts.push(t.label());
        }
        format!("~ {}", parts.join(" + "))
    }
}

impl fmt::Display for Formula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_r_string())
    }
}

/// Parse a formula. `~` is optional and stripped, as `terms(formula(x))` does.
pub fn parse(src: &str) -> Result<Formula, IoError> {
    let mut s = src.trim();
    if let Some(rest) = s.strip_prefix('~') {
        s = rest.trim();
    }
    if s.is_empty() {
        return Err(IoError::Formula("the formula is empty".into()));
    }
    for (op, why) in [
        ('|', "the conditioning operator `|` is not supported"),
        ('/', "the `/` nesting operator is not supported"),
        ('^', "`^` is not supported"),
        ('(', "function calls and `I()` are not supported"),
        ('[', "subsetting is not supported"),
        ('$', "`$` is not supported"),
    ] {
        if s.contains(op) {
            return Err(IoError::Formula(format!("{why} (found {op:?} in {src:?})")));
        }
    }
    if s.contains('.')
        && !s
            .split_whitespace()
            .all(|w| w.chars().all(|c| c.is_ascii_digit()))
    {
        return Err(IoError::Formula(format!(
            "the `.` placeholder is not supported (found in {src:?}); name the \
             covariates explicitly"
        )));
    }

    let mut terms: Vec<Term> = Vec::new();
    // R includes the intercept unless the formula removes it with `-1` or `+0`,
    // so the default is on and only an explicit removal turns it off.
    let mut intercept = true;
    // Removed terms are remembered so a later `+` of the same term does not
    // silently re-add it, which is what `terms()` does.
    let mut removed: Vec<Term> = Vec::new();
    for (sign, raw) in split_top(s) {
        let raw = raw.trim();
        if raw.is_empty() {
            continue;
        }
        if raw == "1" {
            intercept = sign > 0;
            continue;
        }
        if raw == "0" && sign > 0 {
            // `+ 0` is the documented way to write "no intercept".
            intercept = false;
            continue;
        }
        if raw == "-1" {
            continue;
        }
        for t in expand(raw) {
            if t == Term::Intercept {
                intercept = sign > 0;
                continue;
            }
            if sign > 0 {
                removed.retain(|r| *r != t);
                if !terms.contains(&t) {
                    terms.push(t);
                }
            } else {
                terms.retain(|x| *x != t);
                if !removed.contains(&t) {
                    removed.push(t);
                }
            }
        }
    }
    if terms.is_empty() && !intercept {
        return Err(IoError::Formula(format!(
            "{src:?} has no terms after removal"
        )));
    }
    Ok(Formula { terms, intercept })
}

/// Split on the top-level `+` and `-` operators *and* on whitespace, returning
/// each piece with its sign (`+1` or `-1`).
///
/// R treats `+`, `-` and plain whitespace at the top level of a formula
/// identically apart from the sign, so `~ a + b - a` is `b` and `~ a b` is
/// `a + b`. A `-` inside a name (`x-1`) is not an operator, and neither is one
/// inside an interaction, so the scan tracks nesting. Whitespace is *not* a
/// separator inside parentheses, because R has no bracketed term list this
/// parser accepts anyway, and treating it as one would mangle a rejected
/// expression's error message.
fn split_top(s: &str) -> Vec<(i32, String)> {
    let mut out = Vec::new();
    let mut sign = 1i32;
    let mut cur = String::new();
    let mut depth = 0i32;
    for ch in s.chars() {
        match ch {
            '(' | '[' => {
                depth += 1;
                cur.push(ch);
            }
            ')' | ']' => {
                depth -= 1;
                cur.push(ch);
            }
            '+' | '-' if depth == 0 => {
                if !cur.trim().is_empty() {
                    out.push((sign, std::mem::take(&mut cur)));
                } else {
                    cur.clear();
                }
                sign = if ch == '+' { 1 } else { -1 };
            }
            c if c.is_whitespace() && depth == 0 => {
                if !cur.is_empty() {
                    out.push((sign, std::mem::take(&mut cur)));
                }
            }
            _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push((sign, cur));
    }
    out
}

/// Expand one term into the terms R's `terms()` would produce: `a*b` becomes
/// `a + b + a:b`, and `a:b` stays as it is.
fn expand(raw: &str) -> Vec<Term> {
    if let Some((a, b)) = raw.split_once('*') {
        let a = a.trim();
        let b = b.trim();
        let mut out = vec![Term::Variable(a.to_string()), Term::Variable(b.to_string())];
        if !a.is_empty() && !b.is_empty() {
            out.push(Term::Interaction(a.to_string(), b.to_string()));
        }
        return out;
    }
    if let Some((a, b)) = raw.split_once(':') {
        let a = a.trim();
        let b = b.trim();
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        return vec![Term::Interaction(a.to_string(), b.to_string())];
    }
    if raw.contains(|c: char| c.is_ascii_digit()) && raw.parse::<f64>().is_err() {
        // `x1` is a variable name, not a product; only a bare number is a number.
    }
    if raw.parse::<f64>().is_ok() {
        return vec![Term::Intercept];
    }
    vec![Term::Variable(raw.to_string())]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_plain_sum() {
        let f = parse("~ group + x1 + x2").unwrap();
        assert!(f.intercept);
        assert_eq!(
            f.terms,
            vec![
                Term::Variable("group".into()),
                Term::Variable("x1".into()),
                Term::Variable("x2".into())
            ]
        );
    }

    #[test]
    fn expands_star_into_three_terms() {
        let f = parse("~ group*x1").unwrap();
        assert_eq!(
            f.terms,
            vec![
                Term::Variable("group".into()),
                Term::Variable("x1".into()),
                Term::Interaction("group".into(), "x1".into()),
            ]
        );
    }

    #[test]
    fn keeps_an_explicit_interaction() {
        let f = parse("~ a:b").unwrap();
        assert_eq!(f.terms, vec![Term::Interaction("a".into(), "b".into())]);
    }

    #[test]
    fn no_intercept_form() {
        let f = parse("~ 0 + x1").unwrap();
        assert!(!f.intercept);
        assert_eq!(f.terms, vec![Term::Variable("x1".into())]);
        let f = parse("~ x1 - 1").unwrap();
        assert!(!f.intercept);
    }

    #[test]
    fn removal_drops_the_term() {
        let f = parse("~ a + b - a").unwrap();
        assert_eq!(f.terms, vec![Term::Variable("b".into())]);
    }

    #[test]
    fn rejects_what_it_cannot_represent() {
        for src in [
            "~ log(x)", "~ a | b", "~ a/b", "~ . + a", "~ x^2", "~ (a+b)",
        ] {
            assert!(parse(src).is_err(), "{src} should be rejected");
        }
    }

    #[test]
    fn variables_are_reported_in_order_without_repeats() {
        let f = parse("~ g + x1 + g:x2").unwrap();
        assert_eq!(f.variables(), vec!["g", "x1", "x2"]);
    }

    #[test]
    fn whitespace_separates_terms_like_r() {
        let f = parse("group x1 x2").unwrap();
        assert_eq!(
            f.terms,
            vec![
                Term::Variable("group".into()),
                Term::Variable("x1".into()),
                Term::Variable("x2".into())
            ]
        );
        let f = parse("~ a + b").unwrap();
        assert_eq!(f.terms, parse("a b").unwrap().terms);
    }

    #[test]
    fn a_digit_in_a_name_is_not_a_product() {
        let f = parse("~ x1 + x2").unwrap();
        assert_eq!(f.terms.len(), 2);
    }
}
