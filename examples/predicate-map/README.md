This example is ported from the demos in
https://github.com/danshapero/predicates. These provide
nice visualizations of the robustness of the predicates in
this crate, as well as the non-robustness of a naive
implementation of these predicates.

The example computes the predicates in this crate, over a
2D grid of values for one of the inputs, while keeping the
rest of the inputs at fixed values.

For example, we compute the `orient2d` predicate on `(12.0,
12.0)`, `(c[i], c[j])`, `(24.0, 24.0)`; where `i, j` varies
in `0..256`, and `c[i]` is the `i`th float after `0.5`. In
other words, `c[i]` is obtained by starting at `0.5`, and
calling
[`nextafter`](https://docs.rs/float_extras/*/float_extras/f64/fn.nextafter.html)
`i` times.

The inputs are set up so that, if the predicates are
calculated exactly, the output is a `png` with gray values on
the main diagonal, black on the lower-left, and white on
the upper-right side of it. However, the naive versions show
that the predicate is not robust: it switches values on both
sides of the main diagonals, indicating the round-off errors.

# Usage

```
cargo run --example predicate-map naive incircle naive-incircle.png
cargo run --example predicate-map naive orient2d naive-orient2d.png
cargo run --example predicate-map robust incircle robust-incircle.png
cargo run --example predicate-map robust orient2d robust-orient2d.png
```
