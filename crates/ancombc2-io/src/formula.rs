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
    /// An interaction of two or more factors, each already a single variable
    /// name. `*` and `:` are the same operator in R, so `a:b * c` is
    /// `Interaction(["a", "b", "c"])` and its column label is `a:b:c`.
    ///
    /// The flag records whether this term came out of a `*` expansion rather than
    /// being written with `:`, because `model.matrix` codes a factor differently in
    /// the two cases -- see [`Formula::terms`].
    Interaction(Vec<String>, bool),
}

impl Term {
    /// The column names this term contributes, given the metadata's variables.
    pub fn label(&self) -> String {
        match self {
            Term::Intercept => "1".into(),
            Term::Variable(v) => v.clone(),
            Term::Interaction(fs, _) => fs.join(":"),
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
                Term::Interaction(fs, _) => {
                    for v in fs {
                        if !out.contains(v) {
                            out.push(v.clone());
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
        if raw.trim_start_matches([' ', '\t']).starts_with(['*', ':'])
            || raw.trim_end_matches([' ', '\t']).ends_with(['*', ':'])
        {
            return Err(IoError::Formula(format!(
                "{raw:?} is not a term: an interaction operator (`*` or `:`) needs a \
                 variable on both sides"
            )));
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
    // A cross-term's **label** lists its variables in the order they first appear
    // in the formula, not sorted and not written. Both of these are the oracle:
    //
    // ```text
    // > colnames(model.matrix(~ group + x1 + ... + x9 + x10 * x1, m))[-1]
    // [1] "x1:x10"      # x1 is written before x10, so it is labelled first
    // > colnames(model.matrix(~ x10 * x1, m))[-1]
    // [1] "x10:x1"      # here x10 is written first, and wins
    // ```
    //
    // The label is part of the Level A contract, so this is not cosmetic, and
    // "sort the factors" gets the first case right by luck and the second wrong.
    let mut appearance: Vec<String> = Vec::new();
    for t in &terms {
        for v in match t {
            Term::Intercept => continue,
            Term::Variable(v) => std::slice::from_ref(v),
            Term::Interaction(fs, _) => fs.as_slice(),
        } {
            if !appearance.contains(v) {
                appearance.push(v.clone());
            }
        }
    }
    let rank = |v: &String| appearance.iter().position(|a| a == v).unwrap_or(usize::MAX);
    for t in &mut terms {
        if let Term::Interaction(fs, _) = t {
            fs.sort_by_key(rank);
        }
    }

    // R's `terms()` emits terms grouped by interaction order, ascending, and in
    // generation order within an order. The sort is stable so that order is
    // preserved inside each group.
    //
    // This is what puts `group` before `x1:x2` in `~ x1:x2 * group`, and `c`
    // before `a:b` in `~ a:b * c` -- both written the other way round, and both
    // observed on the oracle:
    //
    // ```text
    // > colnames(model.matrix(~ x1:x2 * group, m))
    // [1] "(Intercept)" "group2" "x1:x2" "x1:x2:group2"
    // ```
    //
    // Column order is part of the Level A contract, so this cannot be skipped.
    terms.sort_by_key(interaction_order);
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
    let chars: Vec<char> = s.chars().collect();
    let mut out = Vec::new();
    let mut sign = 1i32;
    let mut cur = String::new();
    let mut depth = 0i32;
    for (i, &ch) in chars.iter().enumerate() {
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
                // Whitespace separates terms, but `*` and `:` bind to their
                // neighbours, so `a * b` is **one** term and not three.
                //
                // Splitting on the whitespace instead made `*` a term of its own,
                // and `expand("*")` then produced `Term::Variable("")` -- twice.
                // The CLI rejected the whole formula with `the formula uses ""`,
                // which names neither the operator nor the term that caused it,
                // and it rejected the exact fixed formula the fixture matrix
                // ships: `group + x1 + ... + x9 + x10 * x1`.
                let next_binds = chars[i + 1..]
                    .iter()
                    .find(|c| !c.is_whitespace())
                    .is_some_and(|c| matches!(c, '*' | ':'));
                let prev_binds = cur.trim_end().ends_with(['*', ':']);
                if next_binds || prev_binds {
                    cur.push(' ');
                } else if !cur.is_empty() {
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

/// Expand one term into the terms R's `terms()` would produce.
///
/// `*` and `:` are **one** operator in R -- `?` `Interaction` -- and `a*b*c` is
/// every variable plus every pairwise product plus the triple. So the piece is
/// first split into `*`-separated factors, each of which is a `:`-chain of
/// variable names, and then every non-empty subset of those factors is added.
///
/// Two rules here were wrong before and are worth stating:
///
/// * `a:b * c` is `c + a:b + a:b:c` -- **not** `a + b + a:b + c + ...`. R keeps
///   `a:b` as a single two-column term, because the formula never asked for the
///   main effects `a` and `b`. Verified against `model.matrix`:
///
///   ```text
///   > colnames(model.matrix(~ x1:x2 * group, m))
///   [1] "(Intercept)" "group2" "x1:x2" "x1:x2:group2"
///   ```
///
/// * the terms come out **grouped by interaction order**, ascending, and in
///   generation order within an order. That is why `group` precedes `x1:x2` in
///   the line above even though it is written second. Generating factors first
///   and cross-products afterwards would put `x1:x2` first and diverge from the
///   oracle's column order, which is part of the Level A contract.
fn expand(raw: &str) -> Vec<Term> {
    let raw = raw.trim();
    if raw.is_empty() {
        return Vec::new();
    }
    // A term written with `:` alone is not a `*` expansion, and `model.matrix`
    // codes a factor inside it differently -- every level survives, where a `*`
    // expansion drops the base level in both the main effect and the interaction.
    let from_star = raw.contains('*');
    let factors: Vec<Vec<String>> = raw
        .split('*')
        .map(|f| {
            f.split(':')
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
                .collect()
        })
        .filter(|f: &Vec<String>| !f.is_empty())
        .collect::<Vec<Vec<String>>>();
    if factors.is_empty() {
        return Vec::new();
    }

    // `from_star` is false for a term written with `:` and true for one this
    // function generated from a `*`, which is what `build_design` needs in order to
    // decide whether a factor inside it drops its base level.
    let as_term = |fs: Vec<String>| -> Term {
        if fs.len() == 1 {
            Term::Variable(fs[0].clone())
        } else {
            Term::Interaction(fs, from_star)
        }
    };

    // Every non-empty *subset* of the factors, smallest first. `a*b*c` is not just
    // the pairs -- R includes the triple:
    //
    // ```text
    // > colnames(model.matrix(~ x1 * x2 * group, m))[-1]
    // [1] "x1" "x2" "group2" "x1:x2" "x1:group2" "x2:group2" "x1:x2:group2"
    // ```
    //
    // `parse` does the stable sort by interaction order, so generating smallest
    // first only fixes the tie-break within an order.
    let k = factors.len();
    let mut out: Vec<Term> = Vec::new();
    for size in 1..=k {
        // Combinations of `size` indices from `0..k`, in lexicographic order:
        // for k = 3 that is {0} {1} {2}, then {0,1} {0,2} {1,2}, then {0,1,2}.
        let mut idx: Vec<usize> = (0..size).collect();
        loop {
            let fs: Vec<String> = idx.iter().flat_map(|&i| factors[i].clone()).collect();
            // `x1 * x1` has no `x1:x1` term: R's `terms()` drops an interaction
            // that repeats a variable, and `model.matrix(~ x1 * x1)` is just `x1`.
            let mut seen: Vec<&String> = Vec::new();
            if fs.iter().any(|v| {
                let dup = seen.contains(&v);
                seen.push(v);
                dup
            }) {
                // keep walking for the next combination rather than emitting
            } else {
                out.push(as_term(fs));
            }
            // Step to the next combination, or stop when this was the last one.
            // The flag is separate from the cursor: `p` legitimately reaches 0 on a
            // successful advance too, and folding the two cases together ended the
            // walk after the *first* combination.
            let mut advanced = false;
            let mut p = size;
            while p > 0 {
                p -= 1;
                if idx[p] != p + k - size {
                    idx[p] += 1;
                    for q in (p + 1)..size {
                        idx[q] = idx[q - 1] + 1;
                    }
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
        }
    }
    out
}

/// The number of variables a term involves, which is what R's `terms()` orders
/// by. Cross-terms are labelled with their factors **sorted**, because R writes
/// `x10:x1` as `x1:x10`: the label is part of the Level A contract.
fn interaction_order(t: &Term) -> usize {
    match t {
        Term::Intercept => 0,
        Term::Variable(_) => 1,
        Term::Interaction(fs, _) => fs.len(),
    }
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

    /// `a * b` with spaces around the operator is **one** term, not three.
    ///
    /// `split_top` also splits on whitespace, which used to leave the `*` as a
    /// term of its own and `expand("*")` then produced `Term::Variable("")` --
    /// twice. The CLI rejected the formula with `the formula uses ""`, naming
    /// neither the operator nor the term, and the exact fixed formula the fixture
    /// matrix ships (`group + x1 + ... + x9 + x10 * x1`) could not be run from the
    /// command line at all.
    #[test]
    fn a_spaced_interaction_is_one_term() {
        let f = parse("~ group + x9 + x10 * x1").unwrap();
        assert_eq!(
            f.variables(),
            vec!["group".to_string(), "x9".into(), "x10".into(), "x1".into()],
            "no variable may come back empty, and `x10` must survive the split"
        );
        assert!(
            f.terms
                .contains(&Term::Interaction(vec!["x10".into(), "x1".into()], true)),
            "the interaction itself is a term, labelled in formula appearance \
             order -- here `x10` is written first; got {f:?}"
        );
        // The same formula spelled without spaces must give the same answer.
        assert_eq!(
            parse("~ group+x9+x10*x1").unwrap().terms,
            f.terms,
            "spacing around `*` is not semantic"
        );
        // The interaction's *label* lists its variables in formula appearance
        // order, so the same interaction is labelled differently depending on
        // which of its variables the formula mentions first. Both are the oracle:
        // `~ x10 * x1` gives `x10:x1`, while the `covariates-10-interaction` cell's
        // `... + x9 + x10 * x1` gives `x1:x10` because `x1` is written first.
        assert_eq!(
            parse("~ group + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 * x1")
                .unwrap()
                .terms
                .last(),
            Some(&Term::Interaction(vec!["x1".into(), "x10".into()], true)),
            "the golden's fix_eff ends in `x1:x10`, so sorting the factors would \
             be right here by luck and wrong for `~ x10 * x1`"
        );
    }

    /// A spaced nested interaction, `a:b * c`, which is the shape R users write.
    #[test]
    fn a_spaced_nested_interaction_keeps_both_colons() {
        let f = parse("~ x1 : x2 * group").unwrap();
        assert_eq!(
            f.terms,
            vec![
                Term::Variable("group".into()),
                Term::Interaction(vec!["x1".into(), "x2".into()], true),
                Term::Interaction(vec!["x1".into(), "x2".into(), "group".into()], true),
            ],
            "`x1:x2 * group` is `group + x1:x2 + x1:x2:group`: `group` first \
             because it is the only order-1 term, and R keeps `x1:x2` as one \
             two-column term rather than asking for the main effects; got {f:?}"
        );
    }

    /// A dangling operator is reported as what it is, rather than as a variable
    /// named `""`.
    #[test]
    fn a_dangling_interaction_operator_names_itself() {
        let e = parse("~ group * ").unwrap_err().to_string();
        assert!(
            e.contains("interaction operator"),
            "the message must name the problem, not an empty variable; got {e:?}"
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
                Term::Interaction(vec!["group".into(), "x1".into()], true),
            ]
        );
    }

    #[test]
    fn keeps_an_explicit_interaction() {
        let f = parse("~ a:b").unwrap();
        assert_eq!(
            f.terms,
            vec![Term::Interaction(vec!["a".into(), "b".into()], false)]
        );
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
