//! Generates the precomputed constants of `sqisign-verify` and `sqisign-rs`.
//!
//! Prime-only constants (the deterministic basis of `E0[2^f]`, the action
//! matrices of the `O0` generators on it, and the quaternion-layer sizes)
//! are derived from `p` alone. The SQIsign round-3 level constants (response
//! length, `D_mix`, `RICofactor`, `EIBox`, encoded sizes) are derived from
//! `(p, lambda)`. See PRECOMP.md for the inventory and the sources of each
//! formula.
//!
//! Usage:
//!   gen-precomp --write ROOT    regenerate every output file
//!   gen-precomp --check ROOT    regenerate and diff (CI)
//!   gen-precomp --dump PRIME              print the values for one prime

mod curve;
mod e0;
mod emit;
mod field;
mod numtheory;

use sqisign_verify::params::{P324_3, P500_27, P664_17};
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.as_slice() {
        [flag, root] if flag == "--write" || flag == "--check" => {
            let root = Path::new(root);
            let outputs = all_outputs();
            let mut stale = 0;
            for (rel, content) in &outputs {
                let path = root.join(rel);
                if flag == "--write" {
                    std::fs::write(&path, content).expect("write output");
                    eprintln!("wrote {}", path.display());
                } else {
                    let have = std::fs::read_to_string(&path).unwrap_or_default();
                    if &have != content {
                        eprintln!("stale: {}", path.display());
                        stale += 1;
                    }
                }
            }
            if stale > 0 {
                std::process::exit(1);
            }
            if flag == "--check" {
                eprintln!("all {} precomp files up to date", outputs.len());
            }
        }
        [flag, name] if flag == "--dump" => {
            let text = match name.as_str() {
                "p324_3" => emit::prime_dump::<P324_3>(),
                "p500_27" => emit::prime_dump::<P500_27>(),
                "p664_17" => emit::prime_dump::<P664_17>(),
                _ => {
                    eprintln!("unknown prime {name}");
                    std::process::exit(2);
                }
            };
            print!("{text}");
        }
        _ => {
            eprintln!("usage: gen-precomp (--write <root> | --check <root> | --dump PRIME)");
            std::process::exit(2);
        }
    }
}

/// Every generated file, as (path relative to the workspace root, content).
fn all_outputs() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let verify_dir = "crates/verify/src/precomp";
    let sign_dir = "crates/sqisign-rs/src/precomp";
    macro_rules! prime_files {
        ($t:ty, $name:literal) => {{
            let data = emit::PrimeData::<$t>::compute($name);
            out.push((
                format!("{verify_dir}/{}.rs", $name),
                emit::rustfmt(data.verify_source()),
            ));
            out.push((
                format!("{sign_dir}/{}.rs", $name),
                emit::rustfmt(data.sign_source()),
            ));
        }};
    }
    prime_files!(P324_3, "p324_3");
    prime_files!(P500_27, "p500_27");
    prime_files!(P664_17, "p664_17");
    out.push((
        "crates/verify/src/params/sqisign_v3.rs".to_string(),
        emit::rustfmt(emit::sqisign_v3_levels()),
    ));
    out
}
