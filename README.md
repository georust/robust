# robust
## Adaptive Precision Floating-Point Arithmetic and Fast Robust Predicates for Computational Geometry

See the [Interactive notebook](https://observablehq.com/@mourner/non-robust-arithmetic-as-art) for more.

[API Documentation](https://docs.rs/robust)

## Visuals

Below are visualizations comparing naive and robust predicate implementations. To learn how these images were generated and how to interpret them, see [`examples/predicate-map/`](examples/predicate-map/).

|            | Naive                   | Robust                   |
|------------|-------------------------|--------------------------|
| `incircle` | ![][incircle-naive-png] | ![][incircle-robust-png] |
| `orient2d` | ![][orient2d-naive-png] | ![][orient2d-robust-png] |

[incircle-naive-png]: https://georust.github.io/assets/incircle-naive/v1.png
[incircle-robust-png]: https://georust.github.io/assets/incircle-robust/v1.png
[orient2d-naive-png]: https://georust.github.io/assets/orient2d-naive/v1.png
[orient2d-robust-png]: https://georust.github.io/assets/orient2d-robust/v1.png

## Source

These algorithms are ported from [`predicates.c`](http://www.cs.cmu.edu/afs/cs/project/quake/public/code/predicates.c), the canonical implementation of Jonathan Richard Shewchuk's "Robust adaptive floating-point geometric predicates".

### Papers

> [Shewchuk, J.R., 1997. Adaptive precision floating-point arithmetic and fast robust geometric predicates. Discrete & Computational Geometry, 18(3), pp.305-363.](https://link.springer.com/content/pdf/10.1007/PL00009321.pdf)

> [Shewchuk, J.R., 1996, May. Robust adaptive floating-point geometric predicates. In Proceedings of the twelfth annual symposium on Computational geometry (pp. 141-150).](https://dl.acm.org/doi/abs/10.1145/237218.237337)


## License

Licensed under either of

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
