use sqisign_rs::keygen::keypair;
use sqisign_rs::sign::sign;
use sqisign_rs::{Level1, Verifier};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: demo-sign <message> [message ...]");
        std::process::exit(1);
    }

    let mut rng = rand::thread_rng();

    eprintln!("Generating keypair...");
    let (pk, sk) = keypair::<Level1>(&mut rng);
    let pk_bytes = pk.to_bytes();

    println!("pk = {}", hex::encode(pk_bytes));
    println!();

    for msg in &args {
        eprintln!("Signing \"{}\"...", msg);
        let sig = sign::<Level1>(&sk, &pk, msg.as_bytes(), &mut rng).expect("signing must succeed");
        let sig_bytes = sig.to_bytes();

        let valid = pk.verify(msg.as_bytes(), &sig).is_ok();

        println!("msg = \"{}\"", msg);
        println!("sig = {}", hex::encode(sig_bytes));
        println!("verified = {}", valid);
        println!();
    }
}
