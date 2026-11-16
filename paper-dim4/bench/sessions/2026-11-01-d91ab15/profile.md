#### Phase profile, two points with the image check (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.93 | 10.6% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 6.77 | 9.1% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.72 | 1.0% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 8.16 | 10.9% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 28.98 | 38.8% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 15.32 | 20.5% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 6.81 | 9.1% |
| sum of the phase medians | 74.71 | 100% |
| median of the whole call (decode and verify) | 75.22 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.51 | |

Inside the plain dimension-4 steps (44.30 Mcycles): point pushes 65.4%, codomain computation 34.6%.

##### Sub-slots of the setup phases, two points with the image check (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.25 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.97 |
| 2: chal: the challenge isogeny chain with its image points | 1.34 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 1.02 |
| 5: chal: P_resc and Q_resc | 0.94 |
| 6: resp: hd_torsion_basis on E_com | 0.25 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.94 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.69 |
| 9: resp: recover_response_cd and the images phi R, phi S | 1.20 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.82 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 4.87 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

#### Phase profile, four points (n = 100; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)

| phase | median (Mcycles) | share |
|---|---|---|
| decode, in the harness (parse_signature, parse_public_key, challenge_limbs) | 0.01 | 0.0% |
| basis generation (recover_challenge, recover_response, hd_torsion_basis) | 7.93 | 7.8% |
| Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp) | 6.77 | 6.7% |
| dimension-2 steps (IsogenyChainDim2::new) | 0.73 | 0.7% |
| gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check) | 8.13 | 8.0% |
| plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs) | 55.81 | 55.0% |
| plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute) | 15.39 | 15.2% |
| final comparisons (middle_codomain_matches) | 0.01 | 0.0% |
| the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder | 6.76 | 6.7% |
| sum of the phase medians | 101.53 | 100% |
| median of the whole call (decode and verify) | 102.10 | |
| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | 0.57 | |

Inside the plain dimension-4 steps (71.20 Mcycles): point pushes 78.4%, codomain computation 21.6%.

##### Sub-slots of the setup phases, four points (n = 100; medians in megacycles; the basis slot is sub-slots 0 to 9, the Cornacchia and kernel setup slot 10 to 15)

| sub-slot | median (Mcycles) |
|---|---|
| 0: chal: hd_torsion_basis on E_pk | 0.25 |
| 1: chal: x-only doublings of the basis and ec_ladder3pt | 0.97 |
| 2: chal: the challenge isogeny chain with its image points | 1.34 |
| 3: chal: the lifts of the image points | 0.27 |
| 4: chal: the Weil pairing(s), the lambda squarings and the sign of phi Q | 1.02 |
| 5: chal: P_resc and Q_resc | 0.94 |
| 6: resp: hd_torsion_basis on E_com | 0.25 |
| 7: resp: the (f - r) doublings to R_com and S_com | 0.93 |
| 8: resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e | 0.69 |
| 9: resp: recover_response_cd and the images phi R, phi S | 1.21 |
| 10: setup: the second hd_torsion_basis on E_com (removed in SECENG-2428, zero since) | 0.00 |
| 11: setup: norm_equation_2f_minus_q (Cornacchia) | 0.03 |
| 12: setup: the four fixed-base doubling chains, the 4-torsion points, make_canonical x2, weil4, the nulls, points_m | 1.80 |
| 13: matrices: starting_two_symplectic_matrices and the two gluing base changes | 0.02 |
| 14: b_kpp, both halves: the kernel generators from the fixed-base chains | 4.87 |
| 15: ThetaStructureDim4::new of the start, both halves | 0.00 |

