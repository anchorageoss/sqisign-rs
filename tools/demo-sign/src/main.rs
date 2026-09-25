//! Sign the command-line arguments under a fresh level I key and verify
//! them, printing the encodings.
use sqisign_rs::{generate, Level1, PublicKey, SigningKey};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: demo-sign <message> [message ...]");
        std::process::exit(1);
    }
    let mut rng = rand::rngs::OsRng;
    eprintln!("Generating a level I key pair...");
    let (pk, sk): (PublicKey<Level1>, SigningKey<Level1>) = generate(&mut rng);
    println!("pk = {}", hex::encode(pk.to_bytes()));
    println!();
    for msg in &args {
        eprintln!("Signing {msg:?}...");
        let sig = sk.sign(msg.as_bytes(), &mut rng).expect("signing");
        let valid = pk.verify(msg.as_bytes(), &sig).is_ok();
        println!("msg = {msg:?}");
        println!("sig = {}", hex::encode(sig.to_bytes()));
        println!("verified = {valid}");
        println!();
    }
}
