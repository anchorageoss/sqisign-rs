#### Phase profile, two points with the image check (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 11.76 | 11.4% |
| Cornacchia and kernel setup | 16.40 | 15.9% |
| dimension-2 steps | 0.72 | 0.7% |
| gluing (and the gluing images of the generators) | 8.13 | 7.9% |
| plain dimension-4 steps: point pushes | 31.71 | 30.8% |
| plain dimension-4 steps: codomain computation | 25.04 | 24.3% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| image check (hd_image_check) | 9.19 | 8.9% |
| sum of the phase medians | 102.96 | 100% |
| median of the whole call (decode and verify) | 103.38 | |
| the whole minus the sum (timing calls and glue) | 0.42 | |

Inside the plain dimension-4 steps (56.75 Mcycles): point pushes 55.9%, codomain computation 44.1%.

#### Phase profile, four points (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 11.76 | 8.9% |
| Cornacchia and kernel setup | 16.39 | 12.4% |
| dimension-2 steps | 0.73 | 0.5% |
| gluing (and the gluing images of the generators) | 8.16 | 6.2% |
| plain dimension-4 steps: point pushes | 60.81 | 46.0% |
| plain dimension-4 steps: codomain computation | 25.15 | 19.0% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| image check (hd_image_check) | 9.19 | 7.0% |
| sum of the phase medians | 132.21 | 100% |
| median of the whole call (decode and verify) | 132.38 | |
| the whole minus the sum (timing calls and glue) | 0.17 | |

Inside the plain dimension-4 steps (85.97 Mcycles): point pushes 70.7%, codomain computation 29.3%.

