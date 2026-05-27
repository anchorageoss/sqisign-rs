use sha2::{Digest, Sha256};
use sqisign_rs::{formats::AnySignature, generate, PublicKey, SigningKey, Verifier};
use std::path::PathBuf;

fn usage() -> ! {
    eprintln!("usage: project-sign <command> [options]");
    eprintln!();
    eprintln!("commands:");
    eprintln!("  generate                          Generate a keypair (hex)");
    eprintln!("  hash                              Compute deterministic source hash");
    eprintln!("  sign    --secret-key <hex> --message <hex>");
    eprintln!("  verify  --public-key <hex> --message <hex> --signature <hex>");
    eprintln!("  update-readme --secret-key <hex>  Hash, sign, and update README.md");
    std::process::exit(1);
}

fn parse_flag<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.windows(2).find_map(|w| {
        if w[0] == flag {
            Some(w[1].as_str())
        } else {
            None
        }
    })
}

fn require_flag<'a>(args: &'a [String], flag: &str) -> &'a str {
    parse_flag(args, flag).unwrap_or_else(|| {
        eprintln!("error: missing required flag {flag}");
        std::process::exit(1);
    })
}

fn cmd_generate() {
    let mut rng = rand::rngs::OsRng;
    let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
    let pk_hex = hex::encode(pk.to_bytes());
    let sk_hex = hex::encode(sk.to_bytes().expect("encoding signing key"));
    println!("PUBLIC_KEY={pk_hex}");
    println!("SECRET_KEY={sk_hex}");
}

fn collect_source_files() -> Vec<PathBuf> {
    let mut files = Vec::new();
    walk_dir("crates".as_ref(), &mut files);
    files.sort();
    files
}

fn walk_dir(dir: &std::path::Path, out: &mut Vec<PathBuf>) {
    let entries = match std::fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            walk_dir(&path, out);
        } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
            out.push(path);
        }
    }
}

fn cmd_hash() {
    println!("{}", compute_hash());
}

fn cmd_sign(args: &[String]) {
    let sk_hex = require_flag(args, "--secret-key");
    let msg_hex = require_flag(args, "--message");

    let sk_bytes = hex::decode(sk_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid secret key hex: {e}");
        std::process::exit(1);
    });
    let msg_bytes = hex::decode(msg_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid message hex: {e}");
        std::process::exit(1);
    });

    let sk: SigningKey = SigningKey::from_bytes(&sk_bytes).unwrap_or_else(|e| {
        eprintln!("error: failed to parse signing key: {e}");
        std::process::exit(1);
    });

    let mut rng = rand::rngs::OsRng;
    let sig = sk.sign(&msg_bytes, &mut rng).unwrap_or_else(|e| {
        eprintln!("error: signing failed: {e}");
        std::process::exit(1);
    });

    let compressed = sig.compress();
    println!("{}", hex::encode(compressed.to_bytes()));
}

fn cmd_verify(args: &[String]) {
    let pk_hex = require_flag(args, "--public-key");
    let msg_hex = require_flag(args, "--message");
    let sig_hex = require_flag(args, "--signature");

    let pk_bytes = hex::decode(pk_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid public key hex: {e}");
        std::process::exit(1);
    });
    let msg_bytes = hex::decode(msg_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid message hex: {e}");
        std::process::exit(1);
    });
    let sig_bytes = hex::decode(sig_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid signature hex: {e}");
        std::process::exit(1);
    });

    let pk: PublicKey = PublicKey::from_bytes(&pk_bytes).unwrap_or_else(|e| {
        eprintln!("error: failed to parse public key: {e}");
        std::process::exit(1);
    });

    let sig = AnySignature::from_bytes(&sig_bytes).unwrap_or_else(|e| {
        eprintln!("error: failed to parse signature: {e}");
        std::process::exit(1);
    });

    match pk.verify(&msg_bytes, &sig) {
        Ok(()) => println!("OK"),
        Err(e) => {
            eprintln!("FAILED: {e}");
            std::process::exit(1);
        }
    }
}

fn wrap_hex(hex: &str, width: usize) -> String {
    hex.as_bytes()
        .chunks(width)
        .map(|c| std::str::from_utf8(c).unwrap())
        .collect::<Vec<_>>()
        .join("\n")
}

fn compute_hash() -> String {
    let files = collect_source_files();
    let mut outer = Sha256::new();
    for path in &files {
        let contents = std::fs::read(path).unwrap_or_else(|e| {
            eprintln!("error reading {}: {e}", path.display());
            std::process::exit(1);
        });
        let hash = Sha256::digest(&contents);
        let line = format!("{}  {}\n", hex::encode(hash), path.display());
        outer.update(line.as_bytes());
    }
    hex::encode(outer.finalize())
}

fn skip_block<'a>(lines: &mut impl Iterator<Item = &'a str>) {
    let mut in_fence = false;
    for l in lines.by_ref() {
        if l.starts_with("```") {
            if in_fence {
                return;
            }
            in_fence = true;
        } else if !in_fence {
            continue;
        }
    }
}

fn cmd_update_readme(args: &[String]) {
    let sk_hex = require_flag(args, "--secret-key");

    let sk_bytes = hex::decode(sk_hex).unwrap_or_else(|e| {
        eprintln!("error: invalid secret key hex: {e}");
        std::process::exit(1);
    });
    let sk: SigningKey = SigningKey::from_bytes(&sk_bytes).unwrap_or_else(|e| {
        eprintln!("error: failed to parse signing key: {e}");
        std::process::exit(1);
    });

    let hash_hex = compute_hash();
    let msg_bytes = hex::decode(&hash_hex).unwrap();

    let mut rng = rand::rngs::OsRng;
    let sig = sk.sign(&msg_bytes, &mut rng).unwrap_or_else(|e| {
        eprintln!("error: signing failed: {e}");
        std::process::exit(1);
    });
    let sig_hex = hex::encode(sig.compress().to_bytes());

    let pk = sk.public_key();
    let sig_raw = hex::decode(&sig_hex).unwrap();
    let any = AnySignature::from_bytes(&sig_raw).unwrap();
    pk.verify(&msg_bytes, &any).unwrap_or_else(|e| {
        eprintln!("error: self-verification failed: {e}");
        std::process::exit(1);
    });

    let readme = std::fs::read_to_string("README.md").unwrap_or_else(|e| {
        eprintln!("error: cannot read README.md: {e}");
        std::process::exit(1);
    });

    let hash_block = format!("**Source hash** (SHA-256):\n```\n{}\n```", hash_hex);
    let sig_block = format!(
        "**Signature** (129 bytes, SQIsign Level 1 compressed):\n```\n{}\n```",
        wrap_hex(&sig_hex, 46)
    );

    let mut out = String::new();
    let mut lines = readme.lines().peekable();
    while let Some(line) = lines.next() {
        if line.starts_with("**Source hash**") {
            out.push_str(&hash_block);
            out.push('\n');
            skip_block(&mut lines);
        } else if line.starts_with("**Signature**") {
            out.push_str(&sig_block);
            out.push('\n');
            skip_block(&mut lines);
        } else {
            out.push_str(line);
            out.push('\n');
        }
    }

    std::fs::write("README.md", &out).unwrap_or_else(|e| {
        eprintln!("error: cannot write README.md: {e}");
        std::process::exit(1);
    });

    eprintln!("Updated README.md");
    eprintln!("  hash: {hash_hex}");
    eprintln!("  sig:  {sig_hex}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }

    match args[1].as_str() {
        "generate" => cmd_generate(),
        "hash" => cmd_hash(),
        "sign" => cmd_sign(&args[2..]),
        "verify" => cmd_verify(&args[2..]),
        "update-readme" => cmd_update_readme(&args[2..]),
        _ => usage(),
    }
}
