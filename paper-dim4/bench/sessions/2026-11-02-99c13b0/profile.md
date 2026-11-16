#### Phase profile, two points with the image check (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.48 | 12.5% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 4.23 | 7.1% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.52 | 0.9% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 3.60 | 6.0% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 26.88 | 44.9% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 13.28 | 22.2% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 3.87 | 6.5% |
| sum of the phase medians | 59.88 | 100% |
| median of the whole call (decode and verify) | 60.37 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.49 | |

Inside the plain dimension-4 steps (40.16 Mcycles): point pushes 66.9%, codomain computation 33.1%.

##### Sub-slots of the setup phases, two points with the image check (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.19 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.96 |
| 2: chal: the challenge isogeny chain with its image points | 1.31 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 1.00 |
| 5: chal: P_resc and Q_resc | 0.93 |
| 6: resp: hd_torsion_basis on E_com | 0.19 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.91 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.66 |
| 9: resp: recover_response_cd and the images phi R, phi S | 0.99 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.61 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 2.54 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

#### Phase profile, four points (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.50 | 8.9% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 4.23 | 5.0% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.52 | 0.6% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 3.62 | 4.3% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 51.47 | 60.9% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 13.35 | 15.8% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 3.86 | 4.6% |
| sum of the phase medians | 84.57 | 100% |
| median of the whole call (decode and verify) | 85.12 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.56 | |

Inside the plain dimension-4 steps (64.82 Mcycles): point pushes 79.4%, codomain computation 20.6%.

##### Sub-slots of the setup phases, four points (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.19 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.96 |
| 2: chal: the challenge isogeny chain with its image points | 1.31 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 1.00 |
| 5: chal: P_resc and Q_resc | 0.93 |
| 6: resp: hd_torsion_basis on E_com | 0.19 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.91 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.67 |
| 9: resp: recover_response_cd and the images phi R, phi S | 1.00 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.60 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 2.53 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

