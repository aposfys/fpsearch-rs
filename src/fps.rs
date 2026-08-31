//! The plain-text fingerprint interchange format.
//!
//! One record per line, `id<TAB>hex`, with `#` comments and blank lines ignored:
//!
//! ```text
//! # n_bits=2048
//! 1<TAB>00a3f1...
//! 2<TAB>4400b2...
//! ```
//!
//! Deliberately boring. Anything that can emit hex — an RDKit script, a shell pipeline —
//! can feed the index builder without linking against this crate, which is what keeps the
//! Rust side from having to own a chemistry toolkit.

use std::fmt;

/// A malformed line, reported with its line number so a 10-million-line file can be fixed.
#[derive(Debug, PartialEq, Eq)]
pub struct ParseError {
    pub line: usize,
    pub kind: ParseErrorKind,
}

#[derive(Debug, PartialEq, Eq)]
pub enum ParseErrorKind {
    /// No tab separating identifier from payload.
    MissingSeparator,
    /// The identifier is not a `u64`.
    BadId(String),
    /// The payload contains a character that is not a hex digit.
    BadHex(char),
    /// The payload has an odd number of hex digits, so it is not a whole number of bytes.
    OddLength(usize),
    /// The payload is a different width from the first record in the file.
    ///
    /// A file mixing widths cannot produce a coherent index, so it is refused here rather
    /// than at query time.
    InconsistentWidth { expected: usize, actual: usize },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: ", self.line)?;
        match &self.kind {
            ParseErrorKind::MissingSeparator => write!(f, "expected 'id<TAB>hex'"),
            ParseErrorKind::BadId(s) => write!(f, "identifier {s:?} is not a u64"),
            ParseErrorKind::BadHex(c) => write!(f, "{c:?} is not a hex digit"),
            ParseErrorKind::OddLength(n) => {
                write!(
                    f,
                    "hex payload has {n} digits, which is not a whole number of bytes"
                )
            }
            ParseErrorKind::InconsistentWidth { expected, actual } => write!(
                f,
                "fingerprint is {actual} bytes but the file started at {expected}; \
                 mixed widths cannot share an index"
            ),
        }
    }
}

impl std::error::Error for ParseError {}

/// One parsed record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub id: u64,
    pub fingerprint: Vec<u64>,
}

fn hex_value(c: char, line: usize) -> Result<u8, ParseError> {
    c.to_digit(16).map(|d| d as u8).ok_or(ParseError {
        line,
        kind: ParseErrorKind::BadHex(c),
    })
}

/// Decode a hex fingerprint into packed `u64` words.
///
/// The first hex digit pair is the first byte of the first word, little-endian within the
/// word, so the same byte order round-trips through [`to_hex`].
pub fn from_hex(hex: &str, line: usize) -> Result<Vec<u64>, ParseError> {
    if hex.len() % 2 != 0 {
        return Err(ParseError {
            line,
            kind: ParseErrorKind::OddLength(hex.len()),
        });
    }
    let mut bytes = Vec::with_capacity(hex.len() / 2);
    let chars: Vec<char> = hex.chars().collect();
    for pair in chars.chunks(2) {
        let hi = hex_value(pair[0], line)?;
        let lo = hex_value(pair[1], line)?;
        bytes.push((hi << 4) | lo);
    }
    // Pad the final partial word with zero bytes, so a width that is not a multiple of
    // eight bytes still packs cleanly.
    while bytes.len() % 8 != 0 {
        bytes.push(0);
    }
    Ok(bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect())
}

/// Encode packed words back to hex.
pub fn to_hex(fingerprint: &[u64]) -> String {
    let mut out = String::with_capacity(fingerprint.len() * 16);
    for word in fingerprint {
        for byte in word.to_le_bytes() {
            out.push_str(&format!("{byte:02x}"));
        }
    }
    out
}

/// Parse a whole `.fps` document.
///
/// Returns the records and the fingerprint width in bits, which the index builder needs and
/// which is checked to be the same on every line.
pub fn parse(text: &str) -> Result<(Vec<Record>, u32), ParseError> {
    let mut records = Vec::new();
    let mut hex_width: Option<usize> = None;

    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let (id_part, hex_part) = trimmed.split_once('\t').ok_or(ParseError {
            line,
            kind: ParseErrorKind::MissingSeparator,
        })?;
        let id: u64 = id_part.trim().parse().map_err(|_| ParseError {
            line,
            kind: ParseErrorKind::BadId(id_part.trim().to_string()),
        })?;
        let hex = hex_part.trim();
        match hex_width {
            None => hex_width = Some(hex.len()),
            Some(expected) if expected != hex.len() => {
                return Err(ParseError {
                    line,
                    kind: ParseErrorKind::InconsistentWidth {
                        expected: expected / 2,
                        actual: hex.len() / 2,
                    },
                })
            }
            _ => {}
        }
        records.push(Record {
            id,
            fingerprint: from_hex(hex, line)?,
        });
    }

    let n_bits = hex_width.unwrap_or(0) as u32 * 4;
    Ok((records, n_bits))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_round_trips() {
        let fp = vec![0x0123_4567_89ab_cdefu64, 0xfedc_ba98_7654_3210u64];
        let hex = to_hex(&fp);
        assert_eq!(from_hex(&hex, 1).unwrap(), fp);
    }

    #[test]
    fn parses_a_small_document() {
        let text = "# a comment\n1\t0f00000000000000\n\n2\tff00000000000000\n";
        let (records, n_bits) = parse(text).unwrap();
        assert_eq!(n_bits, 64);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, 1);
        assert_eq!(crate::popcount(&records[0].fingerprint), 4);
        assert_eq!(crate::popcount(&records[1].fingerprint), 8);
    }

    #[test]
    fn mixed_widths_are_refused() {
        let text = "1\t0f00\n2\tff0000\n";
        let err = parse(text).unwrap_err();
        assert_eq!(err.line, 2);
        assert!(matches!(
            err.kind,
            ParseErrorKind::InconsistentWidth {
                expected: 2,
                actual: 3
            }
        ));
    }

    #[test]
    fn a_bad_hex_digit_names_itself_and_its_line() {
        let text = "1\t0f00\n2\tzz00\n";
        let err = parse(text).unwrap_err();
        assert_eq!(err.line, 2);
        assert_eq!(err.kind, ParseErrorKind::BadHex('z'));
    }

    #[test]
    fn a_missing_tab_is_an_error_not_a_silent_skip() {
        let err = parse("1 0f00\n").unwrap_err();
        assert_eq!(err.kind, ParseErrorKind::MissingSeparator);
    }
}
