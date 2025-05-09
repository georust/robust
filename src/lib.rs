#![cfg_attr(feature = "no_std", no_std)]
#![doc(html_logo_url = "https://raw.githubusercontent.com/georust/meta/master/logo/logo.png")]
// Copyright 2017 The Spade Developers.
// Copyright 2020 The GeoRust Project Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
#![allow(non_snake_case)]

//! # Adaptive Precision Floating-Point Arithmetic and Fast Robust Predicates for Computational Geometry
//! This is a direct transcript of the source code and algorithms provided by
//! Jonathan Richard Shewchuk ([https://www.cs.cmu.edu/~quake/robust.html](https://www.cs.cmu.edu/~quake/robust.html))
//! See the paper and the source code for more information.
//!
//! The module offers adaptive and precise calculations for orientation queries
//! – "on which side of a line (2d) or plane (3d) does a point lie?" – and in-circle / in-sphere queries
//! – "is a given point contained in the circumference of a triangle?".  
//! The "adaptive" nature will increase performance only if a simpler calculation
//! cannot be guaranteed to be accurate enough, yielding higher performance on
//! average.
//!
//! The public API will accept both `f32` and `f64` input points for predicate checking, with input being converted to
//! `f64` values for internal use.
//! This has no effect on precision, as the [IEEE-754 standard](https://drive.google.com/file/d/0B3O3Ys97VjtxYXBCY08wanNoZ1U/view) (section 5.3)
//! guarantees that conversion from `f32` to `f64` must be exact.
//!
//! # Features
//! - `no_std`: Build without the Rust standard library

#[cfg(test)]
mod tests;

/// A two dimensional coordinate.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Coord<T: Into<f64>> {
    pub x: T,
    pub y: T,
}

/// A three dimensional coordinate.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Coord3D<T: Into<f64>> {
    pub x: T,
    pub y: T,
    pub z: T,
}

// These values are precomputed from the "exactinit" method of the c-source code. They should? be
// the same in all IEEE-754 environments, including rust f64
const SPLITTER: f64 = 134_217_729f64;
const EPSILON: f64 = 0.000_000_000_000_000_111_022_302_462_515_65;
const RESULTERRBOUND: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const CCWERRBOUND_B: f64 = (2.0 + 12.0 * EPSILON) * EPSILON;
const CCWERRBOUND_C: f64 = (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
const O3DERRBOUND_A: f64 = (7.0 + 56.0 * EPSILON) * EPSILON;
const O3DERRBOUND_B: f64 = (3.0 + 28.0 * EPSILON) * EPSILON;
const O3DERRBOUND_C: f64 = (26.0 + 288.0 * EPSILON) * EPSILON * EPSILON;
const ICCERRBOUND_A: f64 = (10.0 + 96.0 * EPSILON) * EPSILON;
const ICCERRBOUND_B: f64 = (4.0 + 48.0 * EPSILON) * EPSILON;
const ICCERRBOUND_C: f64 = (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;
const ISPERRBOUND_A: f64 = (16.0 + 224.0 * EPSILON) * EPSILON;
const ISPERRBOUND_B: f64 = (5.0 + 72.0 * EPSILON) * EPSILON;
const ISPERRBOUND_C: f64 = (71.0 + 1408.0 * EPSILON) * EPSILON * EPSILON;

#[inline(always)]
fn abs(x: f64) -> f64 {
    x.abs()
}

/// Returns a positive value if the coordinates `pa`, `pb`, and `pc` occur in counterclockwise order
/// (`pc` lies to the **left** of the directed line defined by coordinates `pa` and `pb`).  
/// Returns a negative value if they occur in clockwise order (`pc` lies to the **right** of the directed line `pa, pb`).  
/// Returns `0` if they are **collinear**.  
pub fn orient2d<T: Into<f64>>(pa: Coord<T>, pb: Coord<T>, pc: Coord<T>) -> f64 {
    let pa = Coord {
        x: pa.x.into(),
        y: pa.y.into(),
    };
    let pb = Coord {
        x: pb.x.into(),
        y: pb.y.into(),
    };
    let pc = Coord {
        x: pc.x.into(),
        y: pc.y.into(),
    };

    let detleft = (pa.x - pc.x) * (pb.y - pc.y);
    let detright = (pa.y - pc.y) * (pb.x - pc.x);
    let det = detleft - detright;

    // The errbound calculation was changed to require only one branch on the likely execution path.
    // This improves performance on modern processors as discussed by Ozaki et al in
    // https://doi.org/10.1007/s10543-015-0574-9
    // The underflow guard "+ u_N" was omitted because orient2dadapt(...) would not guarantee
    // correct results in cases of underflow, the derivation of THETA is given in the reference.
    let detsum = abs(detleft + detright);
    const THETA: f64 = 3.3306690621773722e-16;
    let errbound = THETA * detsum;
    if det >= errbound || -det >= errbound {
        det
    } else {
        orient2dadapt(pa, pb, pc, detsum)
    }
}

fn orient2dadapt(pa: Coord<f64>, pb: Coord<f64>, pc: Coord<f64>, detsum: f64) -> f64 {
    let acx = pa.x - pc.x;
    let bcx = pb.x - pc.x;
    let acy = pa.y - pc.y;
    let bcy = pb.y - pc.y;

    let (detleft, detlefttail) = two_product(acx, bcy);
    let (detright, detrighttail) = two_product(acy, bcx);

    let (B3, B2, B1, B0) = two_two_diff(detleft, detlefttail, detright, detrighttail);
    let B = [B0, B1, B2, B3];

    let mut det = estimate(&B);
    let errbound = CCWERRBOUND_B * detsum;
    if det >= errbound || (-det >= errbound) {
        return det;
    }

    let acxtail = two_diff_tail(pa.x, pc.x, acx);
    let bcxtail = two_diff_tail(pb.x, pc.x, bcx);
    let acytail = two_diff_tail(pa.y, pc.y, acy);
    let bcytail = two_diff_tail(pb.y, pc.y, bcy);

    if acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0 {
        return det;
    }

    let errbound = CCWERRBOUND_C * detsum + RESULTERRBOUND * abs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);

    if det >= errbound || -det >= errbound {
        return det;
    }

    let (s1, s0) = two_product(acxtail, bcy);
    let (t1, t0) = two_product(acytail, bcx);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];

    let mut C1 = [0.0f64; 8];
    let c1length = fast_expansion_sum_zeroelim(&B, &U, &mut C1);

    let (s1, s0) = two_product(acx, bcytail);
    let (t1, t0) = two_product(acy, bcxtail);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];

    let mut C2 = [0.0f64; 12];
    let c2length = fast_expansion_sum_zeroelim(&C1[..c1length], &U, &mut C2);

    let (s1, s0) = two_product(acxtail, bcytail);
    let (t1, t0) = two_product(acytail, bcxtail);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];
    let mut D = [0.0f64; 16];
    let dlength = fast_expansion_sum_zeroelim(&C2[..c2length], &U, &mut D);
    D[dlength - 1]
}

/// Returns a positive value if the point `pd` lies below the plane passing through `pa`, `pb`, and `pc`
/// ("below" is defined so that `pa`, `pb`, and `pc` appear in counterclockwise order when viewed from above the plane).  
/// Returns a negative value if `pd` lies above the plane.  
/// Returns `0` if they are **coplanar**.  
pub fn orient3d<T: Into<f64>>(
    pa: Coord3D<T>,
    pb: Coord3D<T>,
    pc: Coord3D<T>,
    pd: Coord3D<T>,
) -> f64 {
    let pa = Coord3D {
        x: pa.x.into(),
        y: pa.y.into(),
        z: pa.z.into(),
    };
    let pb = Coord3D {
        x: pb.x.into(),
        y: pb.y.into(),
        z: pb.z.into(),
    };
    let pc = Coord3D {
        x: pc.x.into(),
        y: pc.y.into(),
        z: pc.z.into(),
    };
    let pd = Coord3D {
        x: pd.x.into(),
        y: pd.y.into(),
        z: pd.z.into(),
    };

    let adx = pa.x - pd.x;
    let bdx = pb.x - pd.x;
    let cdx = pc.x - pd.x;
    let ady = pa.y - pd.y;
    let bdy = pb.y - pd.y;
    let cdy = pc.y - pd.y;
    let adz = pa.z - pd.z;
    let bdz = pb.z - pd.z;
    let cdz = pc.z - pd.z;

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);

    let permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz)
        + (abs(cdxady) + abs(adxcdy)) * abs(bdz)
        + (abs(adxbdy) + abs(bdxady)) * abs(cdz);

    let errbound = O3DERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return det;
    }

    orient3dadapt(pa, pb, pc, pd, permanent)
}

fn orient3dadapt(
    pa: Coord3D<f64>,
    pb: Coord3D<f64>,
    pc: Coord3D<f64>,
    pd: Coord3D<f64>,
    permanent: f64,
) -> f64 {
    let adx = pa.x - pd.x;
    let bdx = pb.x - pd.x;
    let cdx = pc.x - pd.x;
    let ady = pa.y - pd.y;
    let bdy = pb.y - pd.y;
    let cdy = pc.y - pd.y;
    let adz = pa.z - pd.z;
    let bdz = pb.z - pd.z;
    let cdz = pc.z - pd.z;

    let (bdxcdy1, bdxcdy0) = two_product(bdx, cdy);
    let (cdxbdy1, cdxbdy0) = two_product(cdx, bdy);
    let (bc3, bc2, bc1, bc0) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let bc = [bc0, bc1, bc2, bc3];
    let mut adet = [0f64; 8];
    let alen = scale_expansion_zeroelim(&bc[..4], adz, &mut adet);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let (ca3, ca2, ca1, ca0) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let ca = [ca0, ca1, ca2, ca3];
    let mut bdet = [0f64; 8];
    let blen = scale_expansion_zeroelim(&ca[..4], bdz, &mut bdet);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let (ab3, ab2, ab1, ab0) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let ab = [ab0, ab1, ab2, ab3];
    let mut cdet = [0f64; 8];
    let clen = scale_expansion_zeroelim(&ab[..4], cdz, &mut cdet);

    let mut abdet = [0f64; 16];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut fin1 = [0f64; 192];
    let finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut fin1);

    let mut det = estimate(&fin1[..finlength]);
    let mut errbound = O3DERRBOUND_B * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }

    let adxtail = two_diff_tail(pa.x, pd.x, adx);
    let bdxtail = two_diff_tail(pb.x, pd.x, bdx);
    let cdxtail = two_diff_tail(pc.x, pd.x, cdx);
    let adytail = two_diff_tail(pa.y, pd.y, ady);
    let bdytail = two_diff_tail(pb.y, pd.y, bdy);
    let cdytail = two_diff_tail(pc.y, pd.y, cdy);
    let adztail = two_diff_tail(pa.z, pd.z, adz);
    let bdztail = two_diff_tail(pb.z, pd.z, bdz);
    let cdztail = two_diff_tail(pc.z, pd.z, cdz);

    if (adxtail == 0.0)
        && (bdxtail == 0.0)
        && (cdxtail == 0.0)
        && (adytail == 0.0)
        && (bdytail == 0.0)
        && (cdytail == 0.0)
        && (adztail == 0.0)
        && (bdztail == 0.0)
        && (cdztail == 0.0)
    {
        return det;
    }

    errbound = O3DERRBOUND_C * permanent + RESULTERRBOUND * abs(det);
    det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + adztail * (bdx * cdy - bdy * cdx))
        + (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + bdztail * (cdx * ady - cdy * adx))
        + (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + cdztail * (adx * bdy - ady * bdx));
    if det >= errbound || -det >= errbound {
        return det;
    }

    let mut finnow = fin1;
    let mut finother = [0f64; 192];

    let mut at_b = [0f64; 4];
    let mut at_c = [0f64; 4];
    let mut bt_c = [0f64; 4];
    let mut bt_a = [0f64; 4];
    let mut ct_a = [0f64; 4];
    let mut ct_b = [0f64; 4];
    let at_blen: usize;
    let at_clen: usize;
    let bt_clen: usize;
    let bt_alen: usize;
    let ct_alen: usize;
    let ct_blen: usize;
    if adxtail == 0.0 {
        if adytail == 0.0 {
            at_b[0] = 0.0;
            at_blen = 1;
            at_c[0] = 0.0;
            at_clen = 1;
        } else {
            let negate = -adytail;
            (at_b[1], at_b[0]) = two_product(negate, bdx);
            at_blen = 2;
            (at_c[1], at_c[0]) = two_product(adytail, cdx);
            at_clen = 2;
        }
    } else if adytail == 0.0 {
        (at_b[1], at_b[0]) = two_product(adxtail, bdy);
        at_blen = 2;
        let negate = -adxtail;
        (at_c[1], at_c[0]) = two_product(negate, cdy);
        at_clen = 2;
    } else {
        let (adxt_bdy1, adxt_bdy0) = two_product(adxtail, bdy);
        let (adyt_bdx1, adyt_bdx0) = two_product(adytail, bdx);
        (at_b[3], at_b[2], at_b[1], at_b[0]) =
            two_two_diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0);
        at_blen = 4;
        let (adyt_cdx1, adyt_cdx0) = two_product(adytail, cdx);
        let (adxt_cdy1, adxt_cdy0) = two_product(adxtail, cdy);
        (at_c[3], at_c[2], at_c[1], at_c[0]) =
            two_two_diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0);
        at_clen = 4;
    }
    if bdxtail == 0.0 {
        if bdytail == 0.0 {
            bt_c[0] = 0.0;
            bt_clen = 1;
            bt_a[0] = 0.0;
            bt_alen = 1;
        } else {
            let negate = -bdytail;
            (bt_c[1], bt_c[0]) = two_product(negate, cdx);
            bt_clen = 2;
            (bt_a[1], bt_a[0]) = two_product(bdytail, adx);
            bt_alen = 2;
        }
    } else if bdytail == 0.0 {
        (bt_c[1], bt_c[0]) = two_product(bdxtail, cdy);
        bt_clen = 2;
        let negate = -bdxtail;
        (bt_a[1], bt_a[0]) = two_product(negate, ady);
        bt_alen = 2;
    } else {
        let (bdxt_cdy1, bdxt_cdy0) = two_product(bdxtail, cdy);
        let (bdyt_cdx1, bdyt_cdx0) = two_product(bdytail, cdx);
        (bt_c[3], bt_c[2], bt_c[1], bt_c[0]) =
            two_two_diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0);
        bt_clen = 4;
        let (bdyt_adx1, bdyt_adx0) = two_product(bdytail, adx);
        let (bdxt_ady1, bdxt_ady0) = two_product(bdxtail, ady);
        (bt_a[3], bt_a[2], bt_a[1], bt_a[0]) =
            two_two_diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0);
        bt_alen = 4;
    }
    if cdxtail == 0.0 {
        if cdytail == 0.0 {
            ct_a[0] = 0.0;
            ct_alen = 1;
            ct_b[0] = 0.0;
            ct_blen = 1;
        } else {
            let negate = -cdytail;
            (ct_a[1], ct_a[0]) = two_product(negate, adx);
            ct_alen = 2;
            (ct_b[1], ct_b[0]) = two_product(cdytail, bdx);
            ct_blen = 2;
        }
    } else if cdytail == 0.0 {
        (ct_a[1], ct_a[0]) = two_product(cdxtail, ady);
        ct_alen = 2;
        let negate = -cdxtail;
        (ct_b[1], ct_b[0]) = two_product(negate, bdy);
        ct_blen = 2;
    } else {
        let (cdxt_ady1, cdxt_ady0) = two_product(cdxtail, ady);
        let (cdyt_adx1, cdyt_adx0) = two_product(cdytail, adx);
        (ct_a[3], ct_a[2], ct_a[1], ct_a[0]) =
            two_two_diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0);
        ct_alen = 4;
        let (cdyt_bdx1, cdyt_bdx0) = two_product(cdytail, bdx);
        let (cdxt_bdy1, cdxt_bdy0) = two_product(cdxtail, bdy);
        (ct_b[3], ct_b[2], ct_b[1], ct_b[0]) =
            two_two_diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0);
        ct_blen = 4;
    }

    let mut bct = [0f64; 8];
    let mut cat = [0f64; 8];
    let mut abt = [0f64; 8];
    let mut u = [0f64; 4];
    let mut v = [0f64; 12];
    let mut w = [0f64; 16];
    let mut vlength: usize;
    let mut wlength: usize;

    let bctlen = fast_expansion_sum_zeroelim(&bt_c[..bt_clen], &ct_b[..ct_blen], &mut bct);
    wlength = scale_expansion_zeroelim(&bct[..bctlen], adz, &mut w);
    let mut finlength =
        fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
    ::core::mem::swap(&mut finnow, &mut finother);

    let catlen = fast_expansion_sum_zeroelim(&ct_a[..ct_alen], &at_c[..at_clen], &mut cat);
    wlength = scale_expansion_zeroelim(&cat[..catlen], bdz, &mut w);
    finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
    ::core::mem::swap(&mut finnow, &mut finother);

    let abtlen = fast_expansion_sum_zeroelim(&at_b[..at_blen], &bt_a[..bt_alen], &mut abt);
    wlength = scale_expansion_zeroelim(&abt[..abtlen], cdz, &mut w);
    finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
    ::core::mem::swap(&mut finnow, &mut finother);

    if adztail != 0.0 {
        vlength = scale_expansion_zeroelim(&bc[..4], adztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if bdztail != 0.0 {
        vlength = scale_expansion_zeroelim(&ca[..4], bdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if cdztail != 0.0 {
        vlength = scale_expansion_zeroelim(&ab[..4], cdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }

    if adxtail != 0.0 {
        if bdytail != 0.0 {
            let (adxt_bdyt1, adxt_bdyt0) = two_product(adxtail, bdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if cdztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if cdytail != 0.0 {
            let negate = -adxtail;
            let (adxt_cdyt1, adxt_cdyt0) = two_product(negate, cdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if bdztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
    }
    if bdxtail != 0.0 {
        if cdytail != 0.0 {
            let (bdxt_cdyt1, bdxt_cdyt0) = two_product(bdxtail, cdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if adztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if adytail != 0.0 {
            let negate = -bdxtail;
            let (bdxt_adyt1, bdxt_adyt0) = two_product(negate, adytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_adyt1, bdxt_adyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if cdztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_adyt1, bdxt_adyt0, cdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
    }
    if cdxtail != 0.0 {
        if adytail != 0.0 {
            let (cdxt_adyt1, cdxt_adyt0) = two_product(cdxtail, adytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_adyt1, cdxt_adyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if bdztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_adyt1, cdxt_adyt0, bdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if bdytail != 0.0 {
            let negate = -cdxtail;
            let (cdxt_bdyt1, cdxt_bdyt0) = two_product(negate, bdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_bdyt1, cdxt_bdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if adztail != 0.0 {
                (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_bdyt1, cdxt_bdyt0, adztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
    }

    if adztail != 0.0 {
        wlength = scale_expansion_zeroelim(&bct[..bctlen], adztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if bdztail != 0.0 {
        wlength = scale_expansion_zeroelim(&cat[..catlen], bdztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if cdztail != 0.0 {
        wlength = scale_expansion_zeroelim(&abt[..abtlen], cdztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &w[..wlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }

    finnow[finlength - 1]
}

/// Returns a positive value if the coordinate `pd` lies **inside** the circle passing through `pa`, `pb`, and `pc`.  
/// Returns a negative value if it lies **outside** the circle.  
/// Returns `0` if the four points are **cocircular**.  
/// **Note**: The points `pa`, `pb`, and `pc` must be in **counterclockwise order**, or the sign of the result will be reversed.
pub fn incircle<T: Into<f64>>(pa: Coord<T>, pb: Coord<T>, pc: Coord<T>, pd: Coord<T>) -> f64 {
    let pa = Coord {
        x: pa.x.into(),
        y: pa.y.into(),
    };
    let pb = Coord {
        x: pb.x.into(),
        y: pb.y.into(),
    };
    let pc = Coord {
        x: pc.x.into(),
        y: pc.y.into(),
    };
    let pd = Coord {
        x: pd.x.into(),
        y: pd.y.into(),
    };

    let adx = pa.x - pd.x;
    let bdx = pb.x - pd.x;
    let cdx = pc.x - pd.x;
    let ady = pa.y - pd.y;
    let bdy = pb.y - pd.y;
    let cdy = pc.y - pd.y;

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

    let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
        + (abs(cdxady) + abs(adxcdy)) * blift
        + (abs(adxbdy) + abs(bdxady)) * clift;
    let errbound = ICCERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return det;
    }
    incircleadapt(pa, pb, pc, pd, permanent)
}

fn incircleadapt(
    pa: Coord<f64>,
    pb: Coord<f64>,
    pc: Coord<f64>,
    pd: Coord<f64>,
    permanent: f64,
) -> f64 {
    let mut temp8 = [0f64; 8];
    let mut temp16a = [0f64; 16];
    let mut temp16b = [0f64; 16];
    let mut temp16c = [0f64; 16];
    let mut temp32a = [0f64; 32];
    let mut temp32b = [0f64; 32];
    let mut temp48 = [0f64; 48];
    let mut temp64 = [0f64; 64];

    let adx = pa.x - pd.x;
    let bdx = pb.x - pd.x;
    let cdx = pc.x - pd.x;
    let ady = pa.y - pd.y;
    let bdy = pb.y - pd.y;
    let cdy = pc.y - pd.y;

    let (bdxcdy1, bdxcdy0) = two_product(bdx, cdy);
    let (cdxbdy1, cdxbdy0) = two_product(cdx, bdy);
    let (bc3, bc2, bc1, bc0) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let bc = [bc0, bc1, bc2, bc3];

    let mut axbc = [0f64; 8];
    let axbclen = scale_expansion_zeroelim(&bc, adx, &mut axbc);
    let mut axxbc = [0f64; 16];
    let axxbclen = scale_expansion_zeroelim(&axbc[..axbclen], adx, &mut axxbc);
    let mut aybc = [0f64; 8];
    let aybclen = scale_expansion_zeroelim(&bc, ady, &mut aybc);
    let mut ayybc = [0f64; 16];
    let ayybclen = scale_expansion_zeroelim(&aybc[..aybclen], ady, &mut ayybc);
    let mut adet = [0f64; 32];
    let alen = fast_expansion_sum_zeroelim(&axxbc[0..axxbclen], &ayybc[0..ayybclen], &mut adet);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let (c3, c2, c1, c0) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let ca = [c0, c1, c2, c3];

    let mut bxca = [0f64; 8];
    let bxcalen = scale_expansion_zeroelim(&ca, bdx, &mut bxca);
    let mut bxxca = [0f64; 16];
    let bxxcalen = scale_expansion_zeroelim(&bxca[..bxcalen], bdx, &mut bxxca);
    let mut byca = [0f64; 8];
    let bycalen = scale_expansion_zeroelim(&ca, bdy, &mut byca);
    let mut byyca = [0f64; 16];
    let byycalen = scale_expansion_zeroelim(&byca[..bycalen], bdy, &mut byyca);
    let mut bdet = [0f64; 32];
    let blen = fast_expansion_sum_zeroelim(&bxxca[..bxxcalen], &byyca[0..byycalen], &mut bdet);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let (ab3, ab2, ab1, ab0) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let ab = [ab0, ab1, ab2, ab3];

    let mut cxab = [0f64; 8];
    let cxablen = scale_expansion_zeroelim(&ab, cdx, &mut cxab);
    let mut cxxab = [0f64; 16];
    let cxxablen = scale_expansion_zeroelim(&cxab[..cxablen], cdx, &mut cxxab);
    let mut cyab = [0f64; 8];
    let cyablen = scale_expansion_zeroelim(&ab, cdy, &mut cyab);
    let mut cyyab = [0f64; 16];
    let cyyablen = scale_expansion_zeroelim(&cyab[..cyablen], cdy, &mut cyyab);
    let mut cdet = [0f64; 32];
    let clen = fast_expansion_sum_zeroelim(&cxxab[..cxxablen], &cyyab[..cyyablen], &mut cdet);

    let mut abdet = [0f64; 64];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut fin1 = [0f64; 1152];
    let mut finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut fin1);

    let mut det = estimate(&fin1[..finlength]);
    let errbound = ICCERRBOUND_B * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }

    let adxtail = two_diff_tail(pa.x, pd.x, adx);
    let adytail = two_diff_tail(pa.y, pd.y, ady);
    let bdxtail = two_diff_tail(pb.x, pd.x, bdx);
    let bdytail = two_diff_tail(pb.y, pd.y, bdy);
    let cdxtail = two_diff_tail(pc.x, pd.x, cdx);
    let cdytail = two_diff_tail(pc.y, pd.y, cdy);
    if adxtail == 0.0
        && bdxtail == 0.0
        && cdxtail == 0.0
        && adytail == 0.0
        && bdytail == 0.0
        && cdytail == 0.0
    {
        return det;
    }

    let errbound = ICCERRBOUND_C * permanent + RESULTERRBOUND * abs(det);
    det += ((adx * adx + ady * ady)
        * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
        + ((bdx * bdx + bdy * bdy)
            * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
        + ((cdx * cdx + cdy * cdy)
            * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));

    if det >= errbound || -det >= errbound {
        return det;
    }

    let mut fin2 = [0f64; 1152];

    let aa = if bdxtail != 0.0 || bdytail != 0.0 || cdxtail != 0.0 || cdytail != 0.0 {
        let (adxadx1, adxadx0) = square(adx);
        let (adyady1, adyady0) = square(ady);
        let (aa3, aa2, aa1, aa0) = two_two_sum(adxadx1, adxadx0, adyady1, adyady0);
        [aa0, aa1, aa2, aa3]
    } else {
        [0f64; 4]
    };

    let bb = if cdxtail != 0.0 || cdytail != 0.0 || adxtail != 0.0 || adytail != 0.0 {
        let (bdxbdx1, bdxbdx0) = square(bdx);
        let (bdybdy1, bdybdy0) = square(bdy);
        let (bb3, bb2, bb1, bb0) = two_two_sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0);
        [bb0, bb1, bb2, bb3]
    } else {
        [0f64; 4]
    };

    let cc = if adxtail != 0.0 || adytail != 0.0 || bdxtail != 0.0 || bdytail != 0.0 {
        let (cdxcdx1, cdxcdx0) = square(cdx);
        let (cdycdy1, cdycdy0) = square(cdy);
        let (cc3, cc2, cc1, cc0) = two_two_sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0);
        [cc0, cc1, cc2, cc3]
    } else {
        [0f64; 4]
    };

    let mut axtbclen = 9;
    let mut axtbc = [0f64; 8];
    if adxtail != 0.0 {
        axtbclen = scale_expansion_zeroelim(&bc, adxtail, &mut axtbc);
        let mut temp16a = [0f64; 16];
        let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], 2.0 * adx, &mut temp16a);

        let mut axtcc = [0f64; 8];
        let axtcclen = scale_expansion_zeroelim(&cc, adxtail, &mut axtcc);
        let mut temp16b = [0f64; 16];
        let temp16blen = scale_expansion_zeroelim(&axtcc[..axtcclen], bdy, &mut temp16b);

        let mut axtbb = [0f64; 8];
        let axtbblen = scale_expansion_zeroelim(&bb, adxtail, &mut axtbb);
        let mut temp16c = [0f64; 16];
        let temp16clen = scale_expansion_zeroelim(&axtbb[..axtbblen], -cdy, &mut temp16c);

        let mut temp32a = [0f64; 32];
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let mut temp48 = [0f64; 48];
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2)
    }

    let mut aytbclen = 9;
    let mut aytbc = [0f64; 8];
    if adytail != 0.0 {
        aytbclen = scale_expansion_zeroelim(&bc, adytail, &mut aytbc);
        let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], 2.0 * ady, &mut temp16a);

        let mut aytbb = [0f64; 8];
        let aytbblen = scale_expansion_zeroelim(&bb, adytail, &mut aytbb);
        let temp16blen = scale_expansion_zeroelim(&aytbb[..aytbblen], cdx, &mut temp16b);

        let mut aytcc = [0f64; 8];
        let aytcclen = scale_expansion_zeroelim(&cc, adytail, &mut aytcc);
        let temp16clen = scale_expansion_zeroelim(&aytcc[..aytcclen], -bdx, &mut temp16c);

        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );

        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2)
    }

    let mut bxtcalen = 9;
    let mut bxtca = [0f64; 8];
    if bdxtail != 0.0 {
        bxtcalen = scale_expansion_zeroelim(&ca, bdxtail, &mut bxtca);
        let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], 2.0 * bdx, &mut temp16a);

        let mut bxtaa = [0f64; 8];
        let bxtaalen = scale_expansion_zeroelim(&aa, bdxtail, &mut bxtaa);
        let temp16blen = scale_expansion_zeroelim(&bxtaa[..bxtaalen], cdy, &mut temp16b);

        let mut bxtcc = [0f64; 8];
        let bxtcclen = scale_expansion_zeroelim(&cc, bdxtail, &mut bxtcc);
        let temp16clen = scale_expansion_zeroelim(&bxtcc[..bxtcclen], -ady, &mut temp16c);

        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2)
    }

    let mut bytcalen = 9;
    let mut bytca = [0f64; 8];
    if bdytail != 0.0 {
        bytcalen = scale_expansion_zeroelim(&ca, bdytail, &mut bytca);
        let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], 2.0 * bdy, &mut temp16a);

        let mut bytcc = [0f64; 8];
        let bytcclen = scale_expansion_zeroelim(&cc, bdytail, &mut bytcc);
        let temp16blen = scale_expansion_zeroelim(&bytcc[..bytcclen], adx, &mut temp16b);

        let mut bytaa = [0f64; 8];
        let bytaalen = scale_expansion_zeroelim(&aa, bdytail, &mut bytaa);
        let temp16clen = scale_expansion_zeroelim(&bytaa[..bytaalen], -cdx, &mut temp16c);

        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );

        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2)
    }

    let mut cxtab = [0f64; 8];
    let mut cxtablen = 9;
    if cdxtail != 0.0 {
        cxtablen = scale_expansion_zeroelim(&ab, cdxtail, &mut cxtab);
        let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], 2.0 * cdx, &mut temp16a);

        let mut cxtbb = [0f64; 8];
        let cxtbblen = scale_expansion_zeroelim(&bb, cdxtail, &mut cxtbb);
        let temp16blen = scale_expansion_zeroelim(&cxtbb[..cxtbblen], ady, &mut temp16b);

        let mut cxtaa = [0f64; 8];
        let cxtaalen = scale_expansion_zeroelim(&aa, cdxtail, &mut cxtaa);
        let temp16clen = scale_expansion_zeroelim(&cxtaa[..cxtaalen], -bdy, &mut temp16c);

        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2);
    }

    let mut cytab = [0f64; 8];
    let mut cytablen = 9;
    if cdytail != 0.0 {
        cytablen = scale_expansion_zeroelim(&ab, cdytail, &mut cytab);
        let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], 2.0 * cdy, &mut temp16a);

        let mut cytaa = [0f64; 8];
        let cytaalen = scale_expansion_zeroelim(&aa, cdytail, &mut cytaa);
        let temp16blen = scale_expansion_zeroelim(&cytaa[..cytaalen], bdx, &mut temp16b);

        let mut cytbb = [0f64; 8];
        let cytbblen = scale_expansion_zeroelim(&bb, cdytail, &mut cytbb);
        let temp16clen = scale_expansion_zeroelim(&cytbb[..cytbblen], -adx, &mut temp16c);

        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::core::mem::swap(&mut fin1, &mut fin2);
    }

    if adxtail != 0.0 || adytail != 0.0 {
        let mut bctt = [0f64; 4];
        let mut bct = [0f64; 8];
        let bcttlen;
        let bctlen;
        if bdxtail != 0.0 || bdytail != 0.0 || cdxtail != 0.0 || cdytail != 0.0 {
            let (ti1, ti0) = two_product(bdxtail, cdy);
            let (tj1, tj0) = two_product(bdx, cdytail);
            let (u3, u2, u1, u0) = two_two_sum(ti1, ti0, tj1, tj0);
            let u = [u0, u1, u2, u3];
            let negate = -bdy;
            let (ti1, ti0) = two_product(cdxtail, negate);
            let negate = -bdytail;
            let (tj1, tj0) = two_product(cdx, negate);
            let (v3, v2, v1, v0) = two_two_sum(ti1, ti0, tj1, tj0);
            let v = [v0, v1, v2, v3];
            bctlen = fast_expansion_sum_zeroelim(&u, &v, &mut bct);
            let (ti1, ti0) = two_product(bdxtail, cdytail);
            let (tj1, tj0) = two_product(cdxtail, bdytail);
            let (bctt3, bctt2, bctt1, bctt0) = two_two_diff(ti1, ti0, tj1, tj0);
            bctt = [bctt0, bctt1, bctt2, bctt3];
            bcttlen = 4;
        } else {
            bct[0] = 0.0;
            bctlen = 1;
            bctt[0] = 0.0;
            bcttlen = 1;
        }

        if adxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], adxtail, &mut temp16a);
            let mut axtbct = [0f64; 16];
            let axtbctlen = scale_expansion_zeroelim(&bct[..bctlen], adxtail, &mut axtbct);
            let temp32alen =
                scale_expansion_zeroelim(&axtbct[..axtbctlen], 2.0 * adx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            if bdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&cc, adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }
            if cdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&bb, -adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&axtbct[..axtbctlen], adxtail, &mut temp32a);
            let mut axtbctt = [0f64; 8];
            let axtbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adxtail, &mut axtbctt);
            let temp16alen =
                scale_expansion_zeroelim(&axtbctt[..axtbcttlen], 2.0 * adx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&axtbctt[..axtbcttlen], adxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }

        if adytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], adytail, &mut temp16a);
            let mut aytbct = [0f64; 16];
            let aytbctlen = scale_expansion_zeroelim(&bct[..bctlen], adytail, &mut aytbct);
            let temp32alen =
                scale_expansion_zeroelim(&aytbct[..aytbctlen], 2.0 * ady, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&aytbct[..aytbctlen], adytail, &mut temp32a);
            let mut aytbctt = [0f64; 8];
            let aytbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adytail, &mut aytbctt);
            let temp16alen =
                scale_expansion_zeroelim(&aytbctt[..aytbcttlen], 2.0 * ady, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&aytbctt[..aytbcttlen], adytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }
    }

    if bdxtail != 0.0 || bdytail != 0.0 {
        let mut catt = [0f64; 4];
        let mut cat = [0f64; 8];
        let cattlen;
        let catlen;

        if cdxtail != 0.0 || cdytail != 0.0 || adxtail != 0.0 || adytail != 0.0 {
            let (ti1, ti0) = two_product(cdxtail, ady);
            let (tj1, tj0) = two_product(cdx, adytail);
            let (u3, u2, u1, u0) = two_two_sum(ti1, ti0, tj1, tj0);
            let u = [u0, u1, u2, u3];
            let negate = -cdy;
            let (ti1, ti0) = two_product(adxtail, negate);
            let negate = -cdytail;
            let (tj1, tj0) = two_product(adx, negate);
            let (v3, v2, v1, v0) = two_two_sum(ti1, ti0, tj1, tj0);
            let v = [v0, v1, v2, v3];
            catlen = fast_expansion_sum_zeroelim(&u, &v, &mut cat);

            let (ti1, ti0) = two_product(cdxtail, adytail);
            let (tj1, tj0) = two_product(adxtail, cdytail);
            let (catt3, catt2, catt1, catt0) = two_two_diff(ti1, ti0, tj1, tj0);
            catt = [catt0, catt1, catt2, catt3];
            cattlen = 4;
        } else {
            cat[0] = 0.0;
            catlen = 1;
            catt[0] = 0.0;
            cattlen = 1;
        }

        if bdxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], bdxtail, &mut temp16a);
            let mut bxtcat = [0f64; 16];
            let bxtcatlen = scale_expansion_zeroelim(&cat[..catlen], bdxtail, &mut bxtcat);
            let temp32alen =
                scale_expansion_zeroelim(&bxtcat[..bxtcatlen], 2.0 * bdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            if cdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&aa, bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }
            if adytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&cc, -bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&bxtcat[..bxtcatlen], bdxtail, &mut temp32a);
            let mut bxtcatt = [0f64; 8];
            let bxtcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdxtail, &mut bxtcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bxtcatt[..bxtcattlen], 2.0 * bdx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&bxtcatt[..bxtcattlen], bdxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }
        if bdytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], bdytail, &mut temp16a);
            let mut bytcat = [0f64; 16];
            let bytcatlen = scale_expansion_zeroelim(&cat[..catlen], bdytail, &mut bytcat);
            let temp32alen =
                scale_expansion_zeroelim(&bytcat[..bytcatlen], 2.0 * bdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&bytcat[..bytcatlen], bdytail, &mut temp32a);
            let mut bytcatt = [0f64; 8];
            let bytcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdytail, &mut bytcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bytcatt[..bytcattlen], 2.0 * bdy, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&bytcatt[..bytcattlen], bdytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }
    }

    if cdxtail != 0.0 || cdytail != 0.0 {
        let mut abtt = [0f64; 4];
        let mut abt = [0f64; 8];
        let abttlen;
        let abtlen;

        if adxtail != 0.0 || adytail != 0.0 || bdxtail != 0.0 || bdytail != 0.0 {
            let (ti1, ti0) = two_product(adxtail, bdy);
            let (tj1, tj0) = two_product(adx, bdytail);
            let (u3, u2, u1, u0) = two_two_sum(ti1, ti0, tj1, tj0);
            let u = [u0, u1, u2, u3];
            let negate = -ady;
            let (ti1, ti0) = two_product(bdxtail, negate);
            let negate = -adytail;
            let (tj1, tj0) = two_product(bdx, negate);
            let (v3, v2, v1, v0) = two_two_sum(ti1, ti0, tj1, tj0);
            let v = [v0, v1, v2, v3];
            abtlen = fast_expansion_sum_zeroelim(&u, &v, &mut abt);

            let (ti1, ti0) = two_product(adxtail, bdytail);
            let (tj1, tj0) = two_product(bdxtail, adytail);
            let (abtt3, abtt2, abtt1, abtt0) = two_two_diff(ti1, ti0, tj1, tj0);
            abtt = [abtt0, abtt1, abtt2, abtt3];
            abttlen = 4;
        } else {
            abt[0] = 0.0;
            abtlen = 1;
            abtt[0] = 0.0;
            abttlen = 1;
        }

        if cdxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], cdxtail, &mut temp16a);
            let mut cxtabt = [0f64; 16];
            let cxtabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdxtail, &mut cxtabt);
            let temp32alen =
                scale_expansion_zeroelim(&cxtabt[..cxtabtlen], 2.0 * cdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            if adytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&bb, cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }
            if bdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&aa, -cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::core::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&cxtabt[..cxtabtlen], cdxtail, &mut temp32a);
            let mut cxtabtt = [0f64; 8];
            let cxtabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdxtail, &mut cxtabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cxtabtt[..cxtabttlen], 2.0 * cdx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&cxtabtt[..cxtabttlen], cdxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }
        if cdytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], cdytail, &mut temp16a);
            let mut cytabt = [0f64; 16];
            let cytabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdytail, &mut cytabt);
            let temp32alen =
                scale_expansion_zeroelim(&cytabt[..cytabtlen], 2.0 * cdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&cytabt[..cytabtlen], cdytail, &mut temp32a);
            let mut cytabtt = [0f64; 8];
            let cytabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdytail, &mut cytabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cytabtt[..cytabttlen], 2.0 * cdy, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&cytabtt[..cytabttlen], cdytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp64[..temp64len], &mut fin2);
            ::core::mem::swap(&mut fin1, &mut fin2);
        }
    }
    fin1[finlength - 1]
}

/// Returns a positive value if the point `pe` lies inside the sphere passing through `pa`, `pb`, `pc`, and `pd`.  
/// Returns a negative value if it lies outside.  
/// Returns `0` if the five points are **cospherical**.  
/// **NOTE**: The points `pa`, `pb`, `pc`, and `pd` must be ordered so that they have a positive orientation.
pub fn insphere<T: Into<f64>>(
    pa: Coord3D<T>,
    pb: Coord3D<T>,
    pc: Coord3D<T>,
    pd: Coord3D<T>,
    pe: Coord3D<T>,
) -> f64 {
    let pa = Coord3D {
        x: pa.x.into(),
        y: pa.y.into(),
        z: pa.z.into(),
    };
    let pb = Coord3D {
        x: pb.x.into(),
        y: pb.y.into(),
        z: pb.z.into(),
    };
    let pc = Coord3D {
        x: pc.x.into(),
        y: pc.y.into(),
        z: pc.z.into(),
    };
    let pd = Coord3D {
        x: pd.x.into(),
        y: pd.y.into(),
        z: pd.z.into(),
    };
    let pe = Coord3D {
        x: pe.x.into(),
        y: pe.y.into(),
        z: pe.z.into(),
    };

    let aex = pa.x - pe.x;
    let bex = pb.x - pe.x;
    let cex = pc.x - pe.x;
    let dex = pd.x - pe.x;
    let aey = pa.y - pe.y;
    let bey = pb.y - pe.y;
    let cey = pc.y - pe.y;
    let dey = pd.y - pe.y;
    let aez = pa.z - pe.z;
    let bez = pb.z - pe.z;
    let cez = pc.z - pe.z;
    let dez = pd.z - pe.z;

    let aexbey = aex * bey;
    let bexaey = bex * aey;
    let ab = aexbey - bexaey;
    let bexcey = bex * cey;
    let cexbey = cex * bey;
    let bc = bexcey - cexbey;
    let cexdey = cex * dey;
    let dexcey = dex * cey;
    let cd = cexdey - dexcey;
    let dexaey = dex * aey;
    let aexdey = aex * dey;
    let da = dexaey - aexdey;

    let aexcey = aex * cey;
    let cexaey = cex * aey;
    let ac = aexcey - cexaey;
    let bexdey = bex * dey;
    let dexbey = dex * bey;
    let bd = bexdey - dexbey;

    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;

    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;

    let det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    let aezplus = abs(aez);
    let bezplus = abs(bez);
    let cezplus = abs(cez);
    let dezplus = abs(dez);
    let aexbeyplus = abs(aexbey);
    let bexaeyplus = abs(bexaey);
    let bexceyplus = abs(bexcey);
    let cexbeyplus = abs(cexbey);
    let cexdeyplus = abs(cexdey);
    let dexceyplus = abs(dexcey);
    let dexaeyplus = abs(dexaey);
    let aexdeyplus = abs(aexdey);
    let aexceyplus = abs(aexcey);
    let cexaeyplus = abs(cexaey);
    let bexdeyplus = abs(bexdey);
    let dexbeyplus = abs(dexbey);
    let permanent = ((cexdeyplus + dexceyplus) * bezplus
        + (dexbeyplus + bexdeyplus) * cezplus
        + (bexceyplus + cexbeyplus) * dezplus)
        * alift
        + ((dexaeyplus + aexdeyplus) * cezplus
            + (aexceyplus + cexaeyplus) * dezplus
            + (cexdeyplus + dexceyplus) * aezplus)
            * blift
        + ((aexbeyplus + bexaeyplus) * dezplus
            + (bexdeyplus + dexbeyplus) * aezplus
            + (dexaeyplus + aexdeyplus) * bezplus)
            * clift
        + ((bexceyplus + cexbeyplus) * aezplus
            + (cexaeyplus + aexceyplus) * bezplus
            + (aexbeyplus + bexaeyplus) * cezplus)
            * dlift;
    let errbound = ISPERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return det;
    }

    insphereadapt(pa, pb, pc, pd, pe, permanent)
}

fn insphereadapt<T: Into<f64>>(
    pa: Coord3D<T>,
    pb: Coord3D<T>,
    pc: Coord3D<T>,
    pd: Coord3D<T>,
    pe: Coord3D<T>,
    permanent: f64,
) -> f64 {
    let pa = Coord3D {
        x: pa.x.into(),
        y: pa.y.into(),
        z: pa.z.into(),
    };
    let pb = Coord3D {
        x: pb.x.into(),
        y: pb.y.into(),
        z: pb.z.into(),
    };
    let pc = Coord3D {
        x: pc.x.into(),
        y: pc.y.into(),
        z: pc.z.into(),
    };
    let pd = Coord3D {
        x: pd.x.into(),
        y: pd.y.into(),
        z: pd.z.into(),
    };
    let pe = Coord3D {
        x: pe.x.into(),
        y: pe.y.into(),
        z: pe.z.into(),
    };

    let aex = pa.x - pe.x;
    let bex = pb.x - pe.x;
    let cex = pc.x - pe.x;
    let dex = pd.x - pe.x;
    let aey = pa.y - pe.y;
    let bey = pb.y - pe.y;
    let cey = pc.y - pe.y;
    let dey = pd.y - pe.y;
    let aez = pa.z - pe.z;
    let bez = pb.z - pe.z;
    let cez = pc.z - pe.z;
    let dez = pd.z - pe.z;

    let (aexbey1, aexbey0) = two_product(aex, bey);
    let (bexaey1, bexaey0) = two_product(bex, aey);
    let (ab3, ab2, ab1, ab0) = two_two_diff(aexbey1, aexbey0, bexaey1, bexaey0);
    let ab = [ab0, ab1, ab2, ab3];

    let (bexcey1, bexcey0) = two_product(bex, cey);
    let (cexbey1, cexbey0) = two_product(cex, bey);
    let (bc3, bc2, bc1, bc0) = two_two_diff(bexcey1, bexcey0, cexbey1, cexbey0);
    let bc = [bc0, bc1, bc2, bc3];

    let (cexdey1, cexdey0) = two_product(cex, dey);
    let (dexcey1, dexcey0) = two_product(dex, cey);
    let (cd3, cd2, cd1, cd0) = two_two_diff(cexdey1, cexdey0, dexcey1, dexcey0);
    let cd = [cd0, cd1, cd2, cd3];

    let (dexaey1, dexaey0) = two_product(dex, aey);
    let (aexdey1, aexdey0) = two_product(aex, dey);
    let (da3, da2, da1, da0) = two_two_diff(dexaey1, dexaey0, aexdey1, aexdey0);
    let da = [da0, da1, da2, da3];

    let (aexcey1, aexcey0) = two_product(aex, cey);
    let (cexaey1, cexaey0) = two_product(cex, aey);
    let (ac3, ac2, ac1, ac0) = two_two_diff(aexcey1, aexcey0, cexaey1, cexaey0);
    let ac = [ac0, ac1, ac2, ac3];

    let (bexdey1, bexdey0) = two_product(bex, dey);
    let (dexbey1, dexbey0) = two_product(dex, bey);
    let (bd3, bd2, bd1, bd0) = two_two_diff(bexdey1, bexdey0, dexbey1, dexbey0);
    let bd = [bd0, bd1, bd2, bd3];

    let mut temp8a = [0f64; 8];
    let mut temp8b = [0f64; 8];
    let mut temp8c = [0f64; 8];
    let mut temp16 = [0f64; 16];
    let mut temp24 = [0f64; 24];
    let mut temp48 = [0f64; 48];
    let mut xdet = [0f64; 96];
    let mut ydet = [0f64; 96];
    let mut zdet = [0f64; 96];
    let mut xydet = [0f64; 192];
    let mut adet = [0f64; 288];
    let mut bdet = [0f64; 288];
    let mut cdet = [0f64; 288];
    let mut ddet = [0f64; 288];
    let mut abdet = [0f64; 576];
    let mut cddet = [0f64; 576];
    let mut fin1 = [0f64; 1152];

    let mut temp8alen = scale_expansion_zeroelim(&cd[..4], bez, &mut temp8a);
    let mut temp8blen = scale_expansion_zeroelim(&bd[..4], -cez, &mut temp8b);
    let mut temp8clen = scale_expansion_zeroelim(&bc[..4], dez, &mut temp8c);
    let mut temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let mut temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aex, &mut temp48);
    let mut xlen = scale_expansion_zeroelim(&temp48[..temp48len], -aex, &mut xdet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aey, &mut temp48);
    let mut ylen = scale_expansion_zeroelim(&temp48[..temp48len], -aey, &mut ydet);
    let mut temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aez, &mut temp48);
    let mut zlen = scale_expansion_zeroelim(&temp48[..temp48len], -aez, &mut zdet);
    let mut xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let alen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut adet);

    temp8alen = scale_expansion_zeroelim(&da[..4], cez, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&ac[..4], dez, &mut temp8b);
    temp8clen = scale_expansion_zeroelim(&cd[..4], aez, &mut temp8c);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bex, &mut temp48);
    xlen = scale_expansion_zeroelim(&temp48[..temp48len], bex, &mut xdet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bey, &mut temp48);
    ylen = scale_expansion_zeroelim(&temp48[..temp48len], bey, &mut ydet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bez, &mut temp48);
    zlen = scale_expansion_zeroelim(&temp48[..temp48len], bez, &mut zdet);
    xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let blen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut bdet);

    temp8alen = scale_expansion_zeroelim(&ab[..4], dez, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&bd[..4], aez, &mut temp8b);
    temp8clen = scale_expansion_zeroelim(&da[..4], bez, &mut temp8c);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cex, &mut temp48);
    xlen = scale_expansion_zeroelim(&temp48[..temp48len], -cex, &mut xdet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cey, &mut temp48);
    ylen = scale_expansion_zeroelim(&temp48[..temp48len], -cey, &mut ydet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cez, &mut temp48);
    zlen = scale_expansion_zeroelim(&temp48[..temp48len], -cez, &mut zdet);
    xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let clen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut cdet);

    temp8alen = scale_expansion_zeroelim(&bc[..4], aez, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&ac[..4], -bez, &mut temp8b);
    temp8clen = scale_expansion_zeroelim(&ab[..4], cez, &mut temp8c);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dex, &mut temp48);
    xlen = scale_expansion_zeroelim(&temp48[..temp48len], dex, &mut xdet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dey, &mut temp48);
    ylen = scale_expansion_zeroelim(&temp48[..temp48len], dey, &mut ydet);
    temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dez, &mut temp48);
    zlen = scale_expansion_zeroelim(&temp48[..temp48len], dez, &mut zdet);
    xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let dlen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut ddet);

    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cddet[..cdlen], &mut fin1);

    let mut det = estimate(&fin1[..finlength]);
    let mut errbound = ISPERRBOUND_B * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }

    let aextail = two_diff_tail(pa.x, pe.x, aex);
    let aeytail = two_diff_tail(pa.y, pe.y, aey);
    let aeztail = two_diff_tail(pa.z, pe.z, aez);
    let bextail = two_diff_tail(pb.x, pe.x, bex);
    let beytail = two_diff_tail(pb.y, pe.y, bey);
    let beztail = two_diff_tail(pb.z, pe.z, bez);
    let cextail = two_diff_tail(pc.x, pe.x, cex);
    let ceytail = two_diff_tail(pc.y, pe.y, cey);
    let ceztail = two_diff_tail(pc.z, pe.z, cez);
    let dextail = two_diff_tail(pd.x, pe.x, dex);
    let deytail = two_diff_tail(pd.y, pe.y, dey);
    let deztail = two_diff_tail(pd.z, pe.z, dez);
    if (aextail == 0.0)
        && (aeytail == 0.0)
        && (aeztail == 0.0)
        && (bextail == 0.0)
        && (beytail == 0.0)
        && (beztail == 0.0)
        && (cextail == 0.0)
        && (ceytail == 0.0)
        && (ceztail == 0.0)
        && (dextail == 0.0)
        && (deytail == 0.0)
        && (deztail == 0.0)
    {
        return det;
    }

    errbound = ISPERRBOUND_C * permanent + RESULTERRBOUND * abs(det);
    let abeps = (aex * beytail + bey * aextail) - (aey * bextail + bex * aeytail);
    let bceps = (bex * ceytail + cey * bextail) - (bey * cextail + cex * beytail);
    let cdeps = (cex * deytail + dey * cextail) - (cey * dextail + dex * ceytail);
    let daeps = (dex * aeytail + aey * dextail) - (dey * aextail + aex * deytail);
    let aceps = (aex * ceytail + cey * aextail) - (aey * cextail + cex * aeytail);
    let bdeps = (bex * deytail + dey * bextail) - (bey * dextail + dex * beytail);
    det += (((bex * bex + bey * bey + bez * bez)
        * ((cez * daeps + dez * aceps + aez * cdeps)
            + (ceztail * da3 + deztail * ac3 + aeztail * cd3))
        + (dex * dex + dey * dey + dez * dez)
            * ((aez * bceps - bez * aceps + cez * abeps)
                + (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
        - ((aex * aex + aey * aey + aez * aez)
            * ((bez * cdeps - cez * bdeps + dez * bceps)
                + (beztail * cd3 - ceztail * bd3 + deztail * bc3))
            + (cex * cex + cey * cey + cez * cez)
                * ((dez * abeps + aez * bdeps + bez * daeps)
                    + (deztail * ab3 + aeztail * bd3 + beztail * da3))))
        + 2.0
            * (((bex * bextail + bey * beytail + bez * beztail)
                * (cez * da3 + dez * ac3 + aez * cd3)
                + (dex * dextail + dey * deytail + dez * deztail)
                    * (aez * bc3 - bez * ac3 + cez * ab3))
                - ((aex * aextail + aey * aeytail + aez * aeztail)
                    * (bez * cd3 - cez * bd3 + dez * bc3)
                    + (cex * cextail + cey * ceytail + cez * ceztail)
                        * (dez * ab3 + aez * bd3 + bez * da3)));
    if det >= errbound || -det >= errbound {
        return det;
    }

    insphereexact(pa, pb, pc, pd, pe)
}

fn insphereexact<T: Into<f64>>(
    pa: Coord3D<T>,
    pb: Coord3D<T>,
    pc: Coord3D<T>,
    pd: Coord3D<T>,
    pe: Coord3D<T>,
) -> f64 {
    let pa = Coord3D {
        x: pa.x.into(),
        y: pa.y.into(),
        z: pa.z.into(),
    };
    let pb = Coord3D {
        x: pb.x.into(),
        y: pb.y.into(),
        z: pb.z.into(),
    };
    let pc = Coord3D {
        x: pc.x.into(),
        y: pc.y.into(),
        z: pc.z.into(),
    };
    let pd = Coord3D {
        x: pd.x.into(),
        y: pd.y.into(),
        z: pd.z.into(),
    };
    let pe = Coord3D {
        x: pe.x.into(),
        y: pe.y.into(),
        z: pe.z.into(),
    };

    let (axby1, axby0) = two_product(pa.x, pb.y);
    let (bxay1, bxay0) = two_product(pb.x, pa.y);
    let (ab3, ab2, ab1, ab0) = two_two_diff(axby1, axby0, bxay1, bxay0);
    let ab = [ab0, ab1, ab2, ab3];

    let (bxcy1, bxcy0) = two_product(pb.x, pc.y);
    let (cxby1, cxby0) = two_product(pc.x, pb.y);
    let (bc3, bc2, bc1, bc0) = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);
    let bc = [bc0, bc1, bc2, bc3];

    let (cxdy1, cxdy0) = two_product(pc.x, pd.y);
    let (dxcy1, dxcy0) = two_product(pd.x, pc.y);
    let (cd3, cd2, cd1, cd0) = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);
    let cd = [cd0, cd1, cd2, cd3];

    let (dxey1, dxey0) = two_product(pd.x, pe.y);
    let (exdy1, exdy0) = two_product(pe.x, pd.y);
    let (de3, de2, de1, de0) = two_two_diff(dxey1, dxey0, exdy1, exdy0);
    let de = [de0, de1, de2, de3];

    let (exay1, exay0) = two_product(pe.x, pa.y);
    let (axey1, axey0) = two_product(pa.x, pe.y);
    let (ea3, ea2, ea1, ea0) = two_two_diff(exay1, exay0, axey1, axey0);
    let ea = [ea0, ea1, ea2, ea3];

    let (axcy1, axcy0) = two_product(pa.x, pc.y);
    let (cxay1, cxay0) = two_product(pc.x, pa.y);
    let (ac3, ac2, ac1, ac0) = two_two_diff(axcy1, axcy0, cxay1, cxay0);
    let ac = [ac0, ac1, ac2, ac3];

    let (bxdy1, bxdy0) = two_product(pb.x, pd.y);
    let (dxby1, dxby0) = two_product(pd.x, pb.y);
    let (bd3, bd2, bd1, bd0) = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);
    let bd = [bd0, bd1, bd2, bd3];

    let (cxey1, cxey0) = two_product(pc.x, pe.y);
    let (excy1, excy0) = two_product(pe.x, pc.y);
    let (ce3, ce2, ce1, ce0) = two_two_diff(cxey1, cxey0, excy1, excy0);
    let ce = [ce0, ce1, ce2, ce3];

    let (dxay1, dxay0) = two_product(pd.x, pa.y);
    let (axdy1, axdy0) = two_product(pa.x, pd.y);
    let (da3, da2, da1, da0) = two_two_diff(dxay1, dxay0, axdy1, axdy0);
    let da = [da0, da1, da2, da3];

    let (exby1, exby0) = two_product(pe.x, pb.y);
    let (bxey1, bxey0) = two_product(pb.x, pe.y);
    let (eb3, eb2, eb1, eb0) = two_two_diff(exby1, exby0, bxey1, bxey0);
    let eb = [eb0, eb1, eb2, eb3];

    let mut temp8a = [0f64; 8];
    let mut temp8b = [0f64; 8];
    let mut temp16 = [0f64; 16];
    let mut temp48a = [0f64; 48];
    let mut temp48b = [0f64; 48];

    let mut abc = [0f64; 24];
    let mut bcd = [0f64; 24];
    let mut cde = [0f64; 24];
    let mut dea = [0f64; 24];
    let mut eab = [0f64; 24];
    let mut abd = [0f64; 24];
    let mut bce = [0f64; 24];
    let mut cda = [0f64; 24];
    let mut deb = [0f64; 24];
    let mut eac = [0f64; 24];

    let mut abcd = [0f64; 96];
    let mut bcde = [0f64; 96];
    let mut cdea = [0f64; 96];
    let mut deab = [0f64; 96];
    let mut eabc = [0f64; 96];

    let mut temp192 = [0f64; 192];
    let mut det384x = [0f64; 384];
    let mut det384y = [0f64; 384];
    let mut det384z = [0f64; 384];

    let mut detxy = [0f64; 768];

    let mut adet = [0f64; 1152];
    let mut bdet = [0f64; 1152];
    let mut cdet = [0f64; 1152];
    let mut ddet = [0f64; 1152];
    let mut edet = [0f64; 1152];

    let mut abdet = [0f64; 2304];
    let mut cddet = [0f64; 2304];
    let mut cdedet = [0f64; 3456];

    let mut deter = [0f64; 5760];

    let mut temp8alen = scale_expansion_zeroelim(&bc[..4], pa.z, &mut temp8a);
    let mut temp8blen = scale_expansion_zeroelim(&ac[..4], -pb.z, &mut temp8b);
    let mut temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&ab[..4], pc.z, &mut temp8a);
    let abclen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut abc);

    temp8alen = scale_expansion_zeroelim(&cd[..4], pb.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&bd[..4], -pc.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&bc[..4], pd.z, &mut temp8a);
    let bcdlen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut bcd);

    temp8alen = scale_expansion_zeroelim(&de[..4], pc.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&ce[..4], -pd.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&cd[..4], pe.z, &mut temp8a);
    let cdelen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut cde);

    temp8alen = scale_expansion_zeroelim(&ea[..4], pd.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&da[..4], -pe.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&de[..4], pa.z, &mut temp8a);
    let dealen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut dea);

    temp8alen = scale_expansion_zeroelim(&ab[..4], pe.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&eb[..4], -pa.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&ea[..4], pb.z, &mut temp8a);
    let eablen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut eab);

    temp8alen = scale_expansion_zeroelim(&bd[..4], pa.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&da[..4], pb.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&ab[..4], pd.z, &mut temp8a);
    let abdlen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut abd);

    temp8alen = scale_expansion_zeroelim(&ce[..4], pb.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&eb[..4], pc.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&bc[..4], pe.z, &mut temp8a);
    let bcelen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut bce);

    temp8alen = scale_expansion_zeroelim(&da[..4], pc.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&ac[..4], pd.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&cd[..4], pa.z, &mut temp8a);
    let cdalen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut cda);

    temp8alen = scale_expansion_zeroelim(&eb[..4], pd.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&bd[..4], pe.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    temp8alen = scale_expansion_zeroelim(&de[..4], pb.z, &mut temp8a);
    let deblen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut deb);

    temp8alen = scale_expansion_zeroelim(&ac[..4], pe.z, &mut temp8a);
    temp8blen = scale_expansion_zeroelim(&ce[..4], pa.z, &mut temp8b);
    temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&ea[..4], pc.z, &mut temp8a);
    let eaclen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut eac);

    let mut temp48alen = fast_expansion_sum_zeroelim(&cde[..cdelen], &bce[..bcelen], &mut temp48a);
    let mut temp48blen = fast_expansion_sum_zeroelim(&deb[..deblen], &bcd[..bcdlen], &mut temp48b);
    for i in 0..temp48blen {
        temp48b[i] = -temp48b[i];
    }

    let bcdelen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut bcde);
    let mut xlen = scale_expansion_zeroelim(&bcde[..bcdelen], pa.x, &mut temp192);
    xlen = scale_expansion_zeroelim(&temp192[..xlen], pa.x, &mut det384x);
    let mut ylen = scale_expansion_zeroelim(&bcde[..bcdelen], pa.y, &mut temp192);
    ylen = scale_expansion_zeroelim(&temp192[..ylen], pa.y, &mut det384y);
    let mut zlen = scale_expansion_zeroelim(&bcde[..bcdelen], pa.z, &mut temp192);
    zlen = scale_expansion_zeroelim(&temp192[..zlen], pa.z, &mut det384z);
    let mut xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let alen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut adet);

    temp48alen = fast_expansion_sum_zeroelim(&dea[..dealen], &cda[..cdalen], &mut temp48a);
    temp48blen = fast_expansion_sum_zeroelim(&eac[..eaclen], &cde[..cdelen], &mut temp48b);
    for i in 0..temp48blen {
        temp48b[i] = -temp48b[i];
    }
    let cdealen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut cdea);
    xlen = scale_expansion_zeroelim(&cdea[..cdealen], pb.x, &mut temp192);
    xlen = scale_expansion_zeroelim(&temp192[..xlen], pb.x, &mut det384x);
    ylen = scale_expansion_zeroelim(&cdea[..cdealen], pb.y, &mut temp192);
    ylen = scale_expansion_zeroelim(&temp192[..ylen], pb.y, &mut det384y);
    zlen = scale_expansion_zeroelim(&cdea[..cdealen], pb.z, &mut temp192);
    zlen = scale_expansion_zeroelim(&temp192[..zlen], pb.z, &mut det384z);
    xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let blen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut bdet);

    temp48alen = fast_expansion_sum_zeroelim(&eab[..eablen], &deb[..deblen], &mut temp48a);
    temp48blen = fast_expansion_sum_zeroelim(&abd[..abdlen], &dea[..dealen], &mut temp48b);
    for i in 0..temp48blen {
        temp48b[i] = -temp48b[i];
    }
    let deablen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut deab);
    xlen = scale_expansion_zeroelim(&deab[..deablen], pc.x, &mut temp192);
    xlen = scale_expansion_zeroelim(&temp192[..xlen], pc.x, &mut det384x);
    ylen = scale_expansion_zeroelim(&deab[..deablen], pc.y, &mut temp192);
    ylen = scale_expansion_zeroelim(&temp192[..ylen], pc.y, &mut det384y);
    zlen = scale_expansion_zeroelim(&deab[..deablen], pc.z, &mut temp192);
    zlen = scale_expansion_zeroelim(&temp192[..zlen], pc.z, &mut det384z);
    xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let clen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut cdet);

    temp48alen = fast_expansion_sum_zeroelim(&abc[..abclen], &eac[..eaclen], &mut temp48a);
    temp48blen = fast_expansion_sum_zeroelim(&bce[..bcelen], &eab[..eablen], &mut temp48b);
    for i in 0..temp48blen {
        temp48b[i] = -temp48b[i];
    }
    let eabclen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut eabc);
    xlen = scale_expansion_zeroelim(&eabc[..eabclen], pd.x, &mut temp192);
    xlen = scale_expansion_zeroelim(&temp192[..xlen], pd.x, &mut det384x);
    ylen = scale_expansion_zeroelim(&eabc[..eabclen], pd.y, &mut temp192);
    ylen = scale_expansion_zeroelim(&temp192[..ylen], pd.y, &mut det384y);
    zlen = scale_expansion_zeroelim(&eabc[..eabclen], pd.z, &mut temp192);
    zlen = scale_expansion_zeroelim(&temp192[..zlen], pd.z, &mut det384z);
    xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let dlen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut ddet);

    temp48alen = fast_expansion_sum_zeroelim(&bcd[..bcdlen], &abd[..abdlen], &mut temp48a);
    temp48blen = fast_expansion_sum_zeroelim(&cda[..cdalen], &abc[..abclen], &mut temp48b);
    for i in 0..temp48blen {
        temp48b[i] = -temp48b[i];
    }
    let abcdlen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut abcd);
    xlen = scale_expansion_zeroelim(&abcd[..abcdlen], pe.x, &mut temp192);
    xlen = scale_expansion_zeroelim(&temp192[..xlen], pe.x, &mut det384x);
    ylen = scale_expansion_zeroelim(&abcd[..abcdlen], pe.y, &mut temp192);
    ylen = scale_expansion_zeroelim(&temp192[..ylen], pe.y, &mut det384y);
    zlen = scale_expansion_zeroelim(&abcd[..abcdlen], pe.z, &mut temp192);
    zlen = scale_expansion_zeroelim(&temp192[..zlen], pe.z, &mut det384z);
    xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let elen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut edet);

    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let cdelen = fast_expansion_sum_zeroelim(&cddet[..cdlen], &edet[..elen], &mut cdedet);
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdedet[..cdelen], &mut deter);

    deter[deterlen - 1]
}

fn scale_expansion_zeroelim(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    let (bhi, blo) = split(b);
    let (mut Q, hh) = two_product_presplit(e[0], b, bhi, blo);
    let mut hindex = 0;
    if hh != 0.0 {
        h[hindex] = hh;
        hindex += 1;
    }
    for eindex in 1..e.len() {
        let enow = e[eindex];
        let (product1, product0) = two_product_presplit(enow, b, bhi, blo);
        let (sum, hh) = two_sum(Q, product0);
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1;
        }
        let (new_q, hh) = fast_two_sum(product1, sum);
        Q = new_q;
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1;
        }
    }
    if Q != 0.0 || hindex == 0 {
        h[hindex] = Q;
        hindex += 1;
    }
    hindex
}

#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    (x, two_product_tail(a, b, x))
}

#[inline]
fn two_product_tail(a: f64, b: f64, x: f64) -> f64 {
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    (alo * blo) - err3
}

#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = SPLITTER * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;
    (ahi, alo)
}

#[inline]
fn two_product_presplit(a: f64, b: f64, bhi: f64, blo: f64) -> (f64, f64) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let err1 = x - ahi * bhi;
    let err2 = err1 - alo * bhi;
    let err3 = err2 - ahi * blo;
    let y = alo * blo - err3;
    (x, y)
}

#[inline]
fn two_two_diff(a1: f64, a0: f64, b1: f64, b0: f64) -> (f64, f64, f64, f64) {
    let (j, _r0, x0) = two_one_diff(a1, a0, b0);
    let (x3, x2, x1) = two_one_diff(j, _r0, b1);
    (x3, x2, x1, x0)
}

#[inline]
fn two_one_diff(a1: f64, a0: f64, b: f64) -> (f64, f64, f64) {
    let (i, x0) = two_diff(a0, b);
    let (x2, x1) = two_sum(a1, i);
    (x2, x1, x0)
}

#[inline]
fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    (x, two_diff_tail(a, b, x))
}

#[inline]
fn two_diff_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt = a - x;
    let avirt = x + bvirt;
    let bround = bvirt - b;
    let around = a - avirt;
    around + bround
}

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    (x, two_sum_tail(a, b, x))
}

#[inline]
fn two_sum_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt = x - a;
    let avirt = x - bvirt;
    let bround = b - bvirt;
    let around = a - avirt;
    around + bround
}

fn estimate(e: &[f64]) -> f64 {
    let mut q = e[0];
    for cur in &e[1..] {
        q += *cur;
    }
    q
}

fn fast_expansion_sum_zeroelim(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let mut enow = e[0];
    let mut fnow = f[0];
    let mut eindex = 0;
    let mut findex = 0;
    let mut Qnew;
    let mut hh;
    let mut Q;
    if (fnow > enow) == (fnow > -enow) {
        Q = enow;
        eindex += 1;
    } else {
        Q = fnow;
        findex += 1;
    }

    let mut hindex = 0;
    if eindex < e.len() && findex < f.len() {
        enow = e[eindex];
        fnow = f[findex];
        if (fnow > enow) == (fnow > -enow) {
            let r = fast_two_sum(enow, Q);
            Qnew = r.0;
            hh = r.1;
            eindex += 1;
        } else {
            let r = fast_two_sum(fnow, Q);
            Qnew = r.0;
            hh = r.1;
            findex += 1;
        }
        Q = Qnew;
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1;
        }

        while eindex < e.len() && findex < f.len() {
            enow = e[eindex];
            fnow = f[findex];
            if (fnow > enow) == (fnow > -enow) {
                let r = two_sum(Q, enow);
                Qnew = r.0;
                hh = r.1;
                eindex += 1;
            } else {
                let r = two_sum(Q, fnow);
                Qnew = r.0;
                hh = r.1;
                findex += 1;
            };
            Q = Qnew;
            if hh != 0.0 {
                h[hindex] = hh;
                hindex += 1;
            }
        }
    }

    while eindex < e.len() {
        enow = e[eindex];
        let r = two_sum(Q, enow);
        Qnew = r.0;
        hh = r.1;
        Q = Qnew;
        eindex += 1;
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1
        }
    }

    while findex < f.len() {
        fnow = f[findex];
        let r = two_sum(Q, fnow);
        Qnew = r.0;
        hh = r.1;
        Q = Qnew;
        findex += 1;
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1
        }
    }

    if Q != 0.0 || hindex == 0 {
        h[hindex] = Q;
        hindex += 1;
    }
    hindex
}

#[inline]
fn fast_two_sum_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt = x - a;
    b - bvirt
}

#[inline]
fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    (x, fast_two_sum_tail(a, b, x))
}

#[inline]
fn square_tail(a: f64, x: f64) -> f64 {
    let (ahi, alo) = split(a);
    let err1 = x - ahi * ahi;
    let err3 = err1 - (ahi + ahi) * alo;
    alo * alo - err3
}

#[inline]
fn square(a: f64) -> (f64, f64) {
    let x = a * a;
    (x, square_tail(a, x))
}

#[inline]
fn two_one_sum(a1: f64, a0: f64, b: f64) -> (f64, f64, f64) {
    let (_i, x0) = two_sum(a0, b);
    let (x2, x1) = two_sum(a1, _i);
    (x2, x1, x0)
}

#[inline]
fn two_two_sum(a1: f64, a0: f64, b1: f64, b0: f64) -> (f64, f64, f64, f64) {
    let (_j, _r0, x0) = two_one_sum(a1, a0, b0);
    let (x3, x2, x1) = two_one_sum(_j, _r0, b1);
    (x3, x2, x1, x0)
}

#[inline]
fn two_one_product(a1: f64, a0: f64, b: f64) -> (f64, f64, f64, f64) {
    let (bhi, blo) = split(b);
    let (mut _i, x0) = two_product_presplit(a0, b, bhi, blo);
    let (mut _j, _0) = two_product_presplit(a1, b, bhi, blo);
    let (_k, x1) = two_sum(_i, _0);
    let (x3, x2) = fast_two_sum(_j, _k);
    (x3, x2, x1, x0)
}
