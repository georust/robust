#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: [(f64, f64); 4]| {
    let pa = robust::Coord { x: data[0].0, y: data[0].1 };
    let pb = robust::Coord { x: data[1].0, y: data[1].1 };
    let pc = robust::Coord { x: data[2].0, y: data[2].1 };
    let pd = robust::Coord { x: data[3].0, y: data[3].1 };

    let result = robust::incircle(pa, pb, pc, pd);

    assert!(result.is_finite());
});
