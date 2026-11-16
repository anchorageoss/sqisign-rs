#### Cycles per Fp2 operation (block fp2-ops, once per round; per operation and variant the median over the rounds of the block's medians, and the minimum of its minima; dep = each operation on the previous result, lane16 = sixteen independent elements, composite = one call on one point)

| operation | variant | rounds | median (cycles) | min (cycles) |
|---|---|---|---|---|
| add | dep | 4 | 20.1 | 20.0 |
| add_out | dep | 4 | 20.1 | 20.0 |
| sub | dep | 4 | 19.1 | 19.0 |
| neg | dep | 4 | 22.8 | 22.7 |
| mul | dep | 4 | 213.9 | 213.6 |
| mul_out | dep | 4 | 218.0 | 217.6 |
| sqr | dep | 4 | 130.7 | 130.5 |
| sqr_out | dep | 4 | 132.0 | 131.2 |
| mul_small | dep | 4 | 241.0 | 240.6 |
| half | dep | 4 | 19.2 | 19.1 |
| inv | dep | 4 | 20230.9 | 20204.0 |
| sqrt | dep | 4 | 40134.0 | 40085.0 |
| legendre_plus_add | dep | 4 | 20119.4 | 20092.5 |
| add | lane16 | 4 | 21.0 | 20.9 |
| sub | lane16 | 4 | 20.9 | 20.9 |
| mul | lane16 | 4 | 213.1 | 212.8 |
| sqr | lane16 | 4 | 137.3 | 136.4 |
| hadamard16 | composite | 4 | 1513.8 | 1507.6 |
| pointwise_mul16 | composite | 4 | 3423.9 | 3421.5 |
| pointwise_sqr16 | composite | 4 | 2011.4 | 2008.9 |
| to_squared_theta16 | composite | 4 | 3528.6 | 3521.8 |
| products_of_others16 | composite | 4 | 9658.2 | 9639.3 |
| double | composite | 4 | 14472.1 | 14442.8 |
| image | composite | 4 | 8813.6 | 8800.0 |
| from_kernel_precompute | composite | 4 | 76739.2 | 76493.4 |
