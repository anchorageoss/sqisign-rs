#### Cycles per Fp2 operation (block fp2-ops, once per round; per operation and variant the median over the rounds of the block's medians, and the minimum of its minima; dep = each operation on the previous result, lane16 = sixteen independent elements, composite = one call on one point)

| operation | variant | rounds | median (cycles) | min (cycles) |
|---|---|---|---|---|
| add | dep | 4 | 20.2 | 20.2 |
| add_out | dep | 4 | 20.2 | 20.2 |
| sub | dep | 4 | 19.3 | 19.1 |
| neg | dep | 4 | 22.9 | 22.8 |
| mul | dep | 4 | 215.5 | 215.2 |
| mul_out | dep | 4 | 220.0 | 219.4 |
| sqr | dep | 4 | 131.7 | 131.5 |
| sqr_out | dep | 4 | 132.3 | 132.1 |
| mul_small | dep | 4 | 242.7 | 242.2 |
| half | dep | 4 | 19.3 | 19.3 |
| inv | dep | 4 | 20376.1 | 20360.0 |
| sqrt | dep | 4 | 40432.7 | 40389.1 |
| legendre_plus_add | dep | 4 | 20281.6 | 20248.6 |
| add | lane16 | 4 | 21.1 | 21.1 |
| sub | lane16 | 4 | 21.1 | 21.1 |
| mul | lane16 | 4 | 214.8 | 214.5 |
| sqr | lane16 | 4 | 138.4 | 137.0 |
| hadamard16 | composite | 4 | 1527.3 | 1519.1 |
| pointwise_mul16 | composite | 4 | 3446.1 | 3441.8 |
| pointwise_sqr16 | composite | 4 | 2024.5 | 2020.5 |
| to_squared_theta16 | composite | 4 | 3554.5 | 3545.9 |
| products_of_others16 | composite | 4 | 9767.5 | 9742.7 |
| double | composite | 4 | 14487.6 | 14469.4 |
| image | composite | 4 | 8864.8 | 8851.9 |
| from_kernel_precompute | composite | 4 | 77632.6 | 77460.7 |
