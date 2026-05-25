use crypto_bigint::{Odd, U1024, U128, U256, U512};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use rand::Rng;
use std::hint::black_box;
use std::time::{Duration, Instant};

const ITERS: u64 = 1_000_000;

fn ns_per_op(d: Duration) -> f64 {
    d.as_nanos() as f64 / ITERS as f64
}

fn ratio(nb: Duration, cb: Duration) -> f64 {
    cb.as_nanos() as f64 / nb.as_nanos() as f64
}

fn random_bytes(rng: &mut impl Rng, n: usize) -> Vec<u8> {
    let mut buf = vec![0u8; n];
    rng.fill(&mut buf[..]);
    buf[0] |= 0x80; // ensure high bit set for full width
    buf
}

fn make_nb(bytes: &[u8]) -> BigInt {
    BigInt::from_bytes_be(Sign::Plus, bytes)
}

fn make_u1024(bytes: &[u8]) -> U1024 {
    let mut padded = [0u8; 128];
    let offset = 128 - bytes.len();
    padded[offset..].copy_from_slice(bytes);
    U1024::from_be_slice(&padded)
}

fn make_u512(bytes: &[u8]) -> U512 {
    let mut padded = [0u8; 64];
    let offset = 64 - bytes.len();
    padded[offset..].copy_from_slice(bytes);
    U512::from_be_slice(&padded)
}

fn make_u256(bytes: &[u8]) -> U256 {
    let mut padded = [0u8; 32];
    let offset = 32 - bytes.len();
    padded[offset..].copy_from_slice(bytes);
    U256::from_be_slice(&padded)
}

fn bench_multiply(rng: &mut impl Rng) {
    println!("\n--- Multiply ---");

    // 64-bit
    {
        let a_nb = BigInt::from(123456789i64);
        let b_nb = BigInt::from(987654321i64);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb * &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = U128::from(123456789u64);
        let b_cb = U128::from(987654321u64);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_mul(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "multiply       64-bit    {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb * &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let b_cb = make_u512(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_mul(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "multiply       271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb * &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_mul(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "multiply       1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }
}

fn bench_add_sub(rng: &mut impl Rng) {
    println!("\n--- Add ---");

    // 64-bit
    {
        let a_nb = BigInt::from(123456789i64);
        let b_nb = BigInt::from(987654321i64);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb + &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = U128::from(123456789u64);
        let b_cb = U128::from(987654321u64);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_add(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "add            64-bit    {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb + &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let b_cb = make_u512(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_add(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "add            271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb + &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_add(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "add            1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    println!("\n--- Subtract ---");

    // 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb - &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let b_cb = make_u512(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_sub(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "subtract       271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb - &b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.wrapping_sub(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "subtract       1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }
}

fn bench_divide(rng: &mut impl Rng) {
    println!("\n--- Divide (div_rem) ---");

    // 271-bit / 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_nb.div_rem(&b_nb));
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let b_cb = make_u512(&bytes_b);
        let b_nz = crypto_bigint::NonZero::new(b_cb).unwrap();
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.div_rem(&b_nz));
        }
        let cb = t0.elapsed();

        println!(
            "div_rem        271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit / 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_nb.div_rem(&b_nb));
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let b_nz = crypto_bigint::NonZero::new(b_cb).unwrap();
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.div_rem(&b_nz));
        }
        let cb = t0.elapsed();

        println!(
            "div_rem        1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }
}

fn bench_gcd(rng: &mut impl Rng) {
    println!("\n--- GCD ---");

    // 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_nb.gcd(&b_nb));
        }
        let nb = t0.elapsed();

        let a_cb = make_u256(&bytes_a[..32]);
        let b_cb = make_u256(&bytes_b[..32]);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.gcd(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "gcd            256-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let iters_gcd = 100_000u64;
        let t0 = Instant::now();
        for _ in 0..iters_gcd {
            black_box(a_nb.gcd(&b_nb));
        }
        let nb = t0.elapsed();
        let nb_ns = nb.as_nanos() as f64 / iters_gcd as f64;

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..iters_gcd {
            black_box(a_cb.gcd(&b_cb));
        }
        let cb = t0.elapsed();
        let cb_ns = cb.as_nanos() as f64 / iters_gcd as f64;

        println!(
            "gcd            1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            nb_ns,
            cb_ns,
            cb_ns / nb_ns
        );
    }
}

fn bench_compare(rng: &mut impl Rng) {
    println!("\n--- Compare ---");

    // 271-bit
    {
        let bytes_a = random_bytes(rng, 34);
        let bytes_b = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_nb > b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let b_cb = make_u512(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.gt(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "compare        271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit
    {
        let bytes_a = random_bytes(rng, 128);
        let bytes_b = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let b_nb = make_nb(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_nb > b_nb);
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let b_cb = make_u1024(&bytes_b);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.gt(&b_cb));
        }
        let cb = t0.elapsed();

        println!(
            "compare        1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }
}

fn bench_shift(rng: &mut impl Rng) {
    println!("\n--- Right Shift ---");

    // 271-bit >> 17
    {
        let bytes_a = random_bytes(rng, 34);
        let a_nb = make_nb(&bytes_a);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb >> 17u32);
        }
        let nb = t0.elapsed();

        let a_cb = make_u512(&bytes_a);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.shr(17));
        }
        let cb = t0.elapsed();

        println!(
            "shr(17)        271-bit   {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }

    // 1024-bit >> 17
    {
        let bytes_a = random_bytes(rng, 128);
        let a_nb = make_nb(&bytes_a);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(&a_nb >> 17u32);
        }
        let nb = t0.elapsed();

        let a_cb = make_u1024(&bytes_a);
        let t0 = Instant::now();
        for _ in 0..ITERS {
            black_box(a_cb.shr(17));
        }
        let cb = t0.elapsed();

        println!(
            "shr(17)        1024-bit  {:>8.1}ns/op   {:>8.1}ns/op   {:.2}x",
            ns_per_op(nb),
            ns_per_op(cb),
            ratio(nb, cb)
        );
    }
}

fn bench_modpow(rng: &mut impl Rng) {
    println!("\n--- Modular exponentiation ---");

    // 256-bit modpow (representative of Miller-Rabin witness test)
    {
        let bytes_base = random_bytes(rng, 32);
        let bytes_exp = random_bytes(rng, 32);
        let bytes_mod = random_bytes(rng, 32);

        let base_nb = make_nb(&bytes_base);
        let exp_nb = make_nb(&bytes_exp);
        let mut mod_bytes = bytes_mod.clone();
        mod_bytes[31] |= 1; // ensure odd for modpow
        let mod_nb = make_nb(&mod_bytes);
        let base_nb = base_nb % &mod_nb;

        let iters_modpow = 10_000u64;
        let t0 = Instant::now();
        for _ in 0..iters_modpow {
            black_box(base_nb.modpow(&exp_nb, &mod_nb));
        }
        let nb = t0.elapsed();
        let nb_ns = nb.as_nanos() as f64 / iters_modpow as f64;

        let mod_cb = make_u256(&mod_bytes);
        let base_cb = make_u256(&bytes_base);
        let exp_cb = make_u256(&bytes_exp);
        let params = crypto_bigint::modular::MontyParams::<{ U256::LIMBS }>::new(
            Odd::new(mod_cb).unwrap(),
        );
        let base_m = crypto_bigint::modular::MontyForm::new(&base_cb, params);
        let t0 = Instant::now();
        for _ in 0..iters_modpow {
            black_box(base_m.pow(&exp_cb));
        }
        let cb = t0.elapsed();
        let cb_ns = cb.as_nanos() as f64 / iters_modpow as f64;

        println!(
            "modpow         256-bit   {:>8.0}ns/op   {:>8.0}ns/op   {:.2}x",
            nb_ns, cb_ns, cb_ns / nb_ns
        );
    }

    // 512-bit modpow
    {
        let bytes_base = random_bytes(rng, 64);
        let bytes_exp = random_bytes(rng, 64);
        let bytes_mod = random_bytes(rng, 64);

        let base_nb = make_nb(&bytes_base);
        let exp_nb = make_nb(&bytes_exp);
        let mut mod_bytes = bytes_mod.clone();
        mod_bytes[63] |= 1;
        let mod_nb = make_nb(&mod_bytes);
        let base_nb = base_nb % &mod_nb;

        let iters_modpow = 1_000u64;
        let t0 = Instant::now();
        for _ in 0..iters_modpow {
            black_box(base_nb.modpow(&exp_nb, &mod_nb));
        }
        let nb = t0.elapsed();
        let nb_ns = nb.as_nanos() as f64 / iters_modpow as f64;

        let mod_cb = make_u512(&mod_bytes);
        let base_cb = make_u512(&bytes_base);
        let exp_cb = make_u512(&bytes_exp);
        let params = crypto_bigint::modular::MontyParams::<{ U512::LIMBS }>::new(
            Odd::new(mod_cb).unwrap(),
        );
        let base_m = crypto_bigint::modular::MontyForm::new(&base_cb, params);
        let t0 = Instant::now();
        for _ in 0..iters_modpow {
            black_box(base_m.pow(&exp_cb));
        }
        let cb = t0.elapsed();
        let cb_ns = cb.as_nanos() as f64 / iters_modpow as f64;

        println!(
            "modpow         512-bit   {:>8.0}ns/op   {:>8.0}ns/op   {:.2}x",
            nb_ns, cb_ns, cb_ns / nb_ns
        );
    }
}

fn main() {
    let mut rng = rand::thread_rng();

    println!("num-bigint vs crypto-bigint benchmark");
    println!("Iterations: {} (except where noted)", ITERS);
    println!("Ratio > 1.0 = crypto-bigint slower, < 1.0 = crypto-bigint faster");
    println!("═══════════════════════════════════════════════════════════════════");

    bench_multiply(&mut rng);
    bench_add_sub(&mut rng);
    bench_divide(&mut rng);
    bench_gcd(&mut rng);
    bench_compare(&mut rng);
    bench_shift(&mut rng);
    bench_modpow(&mut rng);

    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Notes:");
    println!("  - crypto-bigint uses fixed-width types (U128/U256/U512/U1024)");
    println!("  - num-bigint dynamically sizes to fit the value");
    println!("  - crypto-bigint compare is constant-time (subtle::CtOption)");
    println!("  - crypto-bigint gcd is constant-time (binary GCD)");
}
