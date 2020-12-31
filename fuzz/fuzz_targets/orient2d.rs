#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: [(f64, f64); 3]| {
    let pa = robust::Coord { x: data[0].0, y: data[0].1 };
    let pb = robust::Coord { x: data[1].0, y: data[1].1 };
    let pc = robust::Coord { x: data[2].0, y: data[2].1 };

    let result = robust::orient2d(pa, pb, pc);

    assert!(result.is_finite());
});
