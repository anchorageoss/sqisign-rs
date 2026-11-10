#### Phase profile, two points with the image check (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.42 | 12.5% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 4.18 | 7.1% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.52 | 0.9% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 3.56 | 6.0% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 26.67 | 44.9% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 13.12 | 22.1% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 3.85 | 6.5% |
| sum of the phase medians | 59.35 | 100% |
| median of the whole call (decode and verify) | 59.90 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.55 | |

Inside the plain dimension-4 steps (39.79 Mcycles): point pushes 67.0%, codomain computation 33.0%.

##### Sub-slots of the setup phases, two points with the image check (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.19 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.96 |
| 2: chal: the challenge isogeny chain with its image points | 1.31 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 0.99 |
| 5: chal: P_resc and Q_resc | 0.92 |
| 6: resp: hd_torsion_basis on E_com | 0.19 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.90 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.66 |
| 9: resp: recover_response_cd and the images phi R, phi S | 0.99 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.59 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 2.51 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

#### Phase profile, four points (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.43 | 8.9% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 4.19 | 5.0% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.52 | 0.6% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 3.59 | 4.3% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 50.98 | 60.8% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 13.23 | 15.8% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 3.85 | 4.6% |
| sum of the phase medians | 83.81 | 100% |
| median of the whole call (decode and verify) | 84.40 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.59 | |

Inside the plain dimension-4 steps (64.21 Mcycles): point pushes 79.4%, codomain computation 20.6%.

##### Sub-slots of the setup phases, four points (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.19 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.96 |
| 2: chal: the challenge isogeny chain with its image points | 1.30 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 0.99 |
| 5: chal: P_resc and Q_resc | 0.92 |
| 6: resp: hd_torsion_basis on E_com | 0.19 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.90 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.66 |
| 9: resp: recover_response_cd and the images phi R, phi S | 0.99 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.59 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 2.50 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

