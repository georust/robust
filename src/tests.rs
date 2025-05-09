//! Test fixtures are from the Staged Geometric Predicates (2001) paper by Aleksandar Nanevski et al
// Fixtures copied from https://github.com/mourner/robust-predicates/tree/main/test/fixtures
// Original location: https://www.cs.cmu.edu/afs/cs/project/pscico/pscico/src/arithmetic/compiler1/test/

use super::{incircle, insphere, orient2d, orient3d, Coord, Coord3D};

#[cfg(not(feature = "no_std"))]
use std::fs::File;
#[cfg(not(feature = "no_std"))]
use std::io::{self, Read};

#[cfg(not(feature = "no_std"))]
fn filename_to_string(s: &str) -> io::Result<String> {
    let mut file = File::open(s)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;
    Ok(s)
}

#[cfg(not(feature = "no_std"))]
fn results_by_line(s: &str) -> Vec<Vec<f64>> {
    s.lines()
        .map(|line| {
            line.split_whitespace()
                .skip(1)
                .map(|foo| foo.parse::<f64>().unwrap())
                .collect()
        })
        .collect()
}

#[cfg(not(feature = "no_std"))]
#[test]
fn test_orient2d_fixtures() {
    let f = filename_to_string("fixtures/orient2d.txt").unwrap();
    let fixtures = results_by_line(&f);
    fixtures.iter().enumerate().for_each(|(idx, fixture)| {
        let ax = fixture[0];
        let ay = fixture[1];
        let bx = fixture[2];
        let by = fixture[3];
        let cx = fixture[4];
        let cy = fixture[5];
        let sign = fixture[6];
        let c1 = Coord { x: ax, y: ay };
        let c2 = Coord { x: bx, y: by };
        let c3 = Coord { x: cx, y: cy };
        let res = orient2d(c1, c2, c3);
        // result sign and fixture sign should be equal
        assert!(
            res.signum() == sign.signum(),
            "Line {line}: Result sign ({result}) and fixture sign ({sign}) should match
            \nCoord 1: {c1:?}\nCoord 2: {c2:?}\nCoord 3: {c3:?}",
            line = idx + 1,
            result = res,
            sign = sign,
            c1 = c1,
            c2 = c2,
            c3 = c3
        );
    })
}

#[cfg(not(feature = "no_std"))]
#[test]
fn test_incircle_fixtures() {
    let f = filename_to_string("fixtures/incircle.txt").unwrap();
    let fixtures = results_by_line(&f);
    fixtures.iter().enumerate().for_each(|(idx, fixture)| {
        let ax = fixture[0];
        let ay = fixture[1];
        let bx = fixture[2];
        let by = fixture[3];
        let cx = fixture[4];
        let cy = fixture[5];
        let dx = fixture[6];
        let dy = fixture[7];
        let sign = fixture[8];
        let c1 = Coord { x: ax, y: ay };
        let c2 = Coord { x: bx, y: by };
        let c3 = Coord { x: cx, y: cy };
        let c4 = Coord { x: dx, y: dy };
        let res = incircle(c1, c2, c3, c4);
        assert!(
            res.signum() == sign.signum(),
            "Line {line}: Result sign ({result}) and fixture sign ({sign}) should match
            \nCoord 1: {c1:?}\nCoord 2: {c2:?}\nCoord 3: {c3:?}\nCoord 4: {c4:?}",
            line = idx + 1,
            result = res,
            sign = sign,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            c4 = c4
        );
    })
}

#[cfg(not(feature = "no_std"))]
#[test]
fn test_orient3d_fixtures() {
    let f = filename_to_string("fixtures/orient3d.txt").unwrap();
    let fixtures = results_by_line(&f);
    fixtures.iter().enumerate().for_each(|(idx, fixture)| {
        let ax = fixture[0];
        let ay = fixture[1];
        let az = fixture[2];

        let bx = fixture[3];
        let by = fixture[4];
        let bz = fixture[5];

        let cx = fixture[6];
        let cy = fixture[7];
        let cz = fixture[8];

        let dx = fixture[9];
        let dy = fixture[10];
        let dz = fixture[11];

        let sign = fixture[12];

        let c1 = Coord3D {
            x: ax,
            y: ay,
            z: az,
        };
        let c2 = Coord3D {
            x: bx,
            y: by,
            z: bz,
        };
        let c3 = Coord3D {
            x: cx,
            y: cy,
            z: cz,
        };
        let c4 = Coord3D {
            x: dx,
            y: dy,
            z: dz,
        };
        let res = orient3d(c1, c2, c3, c4);
        assert!(
            res.signum() == sign.signum(),
            "Line {line}: Result sign ({result}) and fixture sign ({sign}) should match
            \nCoord 1: {c1:?}\nCoord 2: {c2:?}\nCoord 3: {c3:?}\nCoord 4: {c4:?}",
            line = idx + 1,
            result = res,
            sign = sign,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            c4 = c4
        );
        // symmetry
        assert!(res.signum() == orient3d(c1, c2, c3, c4).signum());
    })
}

#[cfg(not(feature = "no_std"))]
#[test]
fn test_insphere_fixtures() {
    let f = filename_to_string("fixtures/insphere.txt").unwrap();
    let fixtures = results_by_line(&f);
    fixtures.iter().enumerate().for_each(|(idx, fixture)| {
        let ax = fixture[0];
        let ay = fixture[1];
        let az = fixture[2];

        let bx = fixture[3];
        let by = fixture[4];
        let bz = fixture[5];

        let cx = fixture[6];
        let cy = fixture[7];
        let cz = fixture[8];

        let dx = fixture[9];
        let dy = fixture[10];
        let dz = fixture[11];

        let ex = fixture[12];
        let ey = fixture[13];
        let ez = fixture[14];

        let sign = fixture[15];

        let c1 = Coord3D {
            x: ax,
            y: ay,
            z: az,
        };
        let c2 = Coord3D {
            x: bx,
            y: by,
            z: bz,
        };
        let c3 = Coord3D {
            x: cx,
            y: cy,
            z: cz,
        };
        let c4 = Coord3D {
            x: dx,
            y: dy,
            z: dz,
        };
        let c5 = Coord3D {
            x: ex,
            y: ey,
            z: ez,
        };
        let res = insphere(c1, c2, c3, c4, c5);
        assert!(
            res.signum() == sign.signum(),
            "Line {line}: Result sign ({result}) and fixture sign ({sign}) should match
            \nCoord 1: {c1:?}\nCoord 2: {c2:?}\nCoord 3: {c3:?}\nCoord 4: {c4:?}\nCoord 5: {c5:?}",
            line = idx + 1,
            result = res,
            sign = sign,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            c4 = c4,
            c5 = c5
        );
    })
}

#[test]
fn test_orient2d() {
    let from = Coord { x: -1f64, y: -1.0 };
    let to = Coord { x: 1f64, y: 1.0 };
    let p1 = Coord {
        x: ::core::f64::MIN_POSITIVE,
        y: ::core::f64::MIN_POSITIVE,
    };
    let p2 = Coord {
        x: -::core::f64::MIN_POSITIVE,
        y: -::core::f64::MIN_POSITIVE,
    };
    let p3 = Coord {
        x: -::core::f64::MIN_POSITIVE,
        y: ::core::f64::MIN_POSITIVE,
    };
    let p4 = Coord {
        x: ::core::f64::MIN_POSITIVE,
        y: -::core::f64::MIN_POSITIVE,
    };

    for &(p, sign) in &[(p1, 0.0), (p2, 0.0), (p3, 1.0), (p4, -1.0)] {
        let det = orient2d(from, to, p);
        assert!(det == sign || det.signum() == sign.signum());
    }
}

#[test]
fn test_orient3d() {
    // plane
    let pa = Coord3D {
        x: 1.,
        y: 0.,
        z: 1.,
    };
    let pb = Coord3D {
        x: -1.,
        y: 0.,
        z: -1.,
    };
    let pc = Coord3D {
        x: -1.,
        y: 0.,
        z: 0.,
    };

    // above plane - negative value expected
    let p1 = Coord3D {
        x: ::core::f64::MIN_POSITIVE,
        y: ::core::f64::MIN_POSITIVE,
        z: ::core::f64::MIN_POSITIVE,
    };
    // below plane - positive value expected
    let p2 = Coord3D {
        x: -::core::f64::MIN_POSITIVE,
        y: -::core::f64::MIN_POSITIVE,
        z: -::core::f64::MIN_POSITIVE,
    };
    // collinear to plane - zero expected
    let p3 = Coord3D {
        x: 0.,
        y: 0.,
        z: 0.,
    };

    for &(p, sign) in &[(p1, -1.0), (p2, 1.0), (p3, 0.0)] {
        let det = orient3d(pa, pb, pc, p);
        assert!(det == sign || det.signum() == sign.signum());
    }
}

#[test]
fn test_incircle() {
    let from = Coord { x: -1f64, y: -1.0 };
    let to = Coord { x: 1f64, y: 1.0 };
    let p_left = Coord {
        x: -::core::f64::MIN_POSITIVE,
        y: ::core::f64::MIN_POSITIVE,
    };
    let p_right = Coord {
        x: ::core::f64::MIN_POSITIVE,
        y: -::core::f64::MIN_POSITIVE,
    };
    let p_query = Coord { x: 2.0, y: 2.0 };

    assert!(incircle(from, p_left, to, p_query) > 0.0);
    assert!(incircle(from, to, p_right, p_query) > 0.0);
}

#[test]
fn test_insphere() {
    let pa = Coord3D {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let pb = Coord3D {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };
    let pc = Coord3D {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    let pd = Coord3D {
        x: 0.0,
        y: -1.0,
        z: 0.0,
    };

    // point outside sphere
    let pe1 = Coord3D {
        x: -1.01,
        y: 0.,
        z: 0.,
    };
    // point inside sphere
    let pe2 = Coord3D {
        x: 0.,
        y: 0.,
        z: 0.99,
    };
    // cospherical point
    let pe3 = Coord3D {
        x: 0.,
        y: 0.,
        z: -1.,
    };

    assert!(insphere(pa, pb, pc, pd, pe1) < 0.0);
    assert!(insphere(pa, pb, pc, pd, pe2) > 0.0);
    assert!(insphere(pa, pb, pc, pd, pe3) == 0.0);
}

#[test]
fn test_issue48_a() {
    let pa = Coord {
        x: 2.1045541600524288e-15,
        y: -1.0000000000000016,
    };
    let pb = Coord {
        x: 1.000000000000005,
        y: -3.350874324301223e-16,
    };
    let pc = Coord {
        x: 7.553997323229233e-15,
        y: 0.9999999999999958,
    };
    let pd = Coord {
        x: -0.9999999999999922,
        y: -7.073397829693697e-15,
    };
    // the (incorrect) result from the previous version of exactpred.rs
    assert!(incircle(pa, pb, pc, pd) != 1.9217716744382023e-16f64);
    // the result predicates.c gives
    assert!(incircle(pa, pb, pc, pd) == -8.0140565430358e-30f64);
}

#[test]
fn test_issue48_b() {
    let pa = Coord {
        x: 9.128561612013288e-15,
        y: -1.0000000000000029,
    };
    let pb = Coord {
        x: 1.0000000000000044,
        y: -5.451395142523081e-15,
    };
    let pc = Coord {
        x: 3.851214418148064e-15,
        y: 0.9999999999999961,
    };
    let pd = Coord {
        x: -0.9999999999999946,
        y: -6.6797960341085084e-15,
    };
    // the (incorrect) result from the previous version of exactpred.rs
    assert!(incircle(pa, pb, pc, pd) != -1.1074731814540733e-16);
    // the result predicates.c gives
    assert!(incircle(pa, pb, pc, pd) == 7.226864249343135e-30);
}
