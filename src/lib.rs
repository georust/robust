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

//! This is a direct transcript of the sourcecode and algorithms provided by
//! Jonathan Richard Shewchuk ([https://www.cs.cmu.edu/~quake/robust.html](https://www.cs.cmu.edu/~quake/robust.html))
//! See the paper and the source code for more information.
//!
//! The module offers adaptive and precise calculations for orientation queries
//! (on which side of a line does a point lie?) and in-circle queries
//! (is a given point contained in the circumference of a triangle?)
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
//!
//! - `no_std`: Build without the Rust standard library

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
const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;
const CCWERRBOUND_B: f64 = (2.0 + 12.0 * EPSILON) * EPSILON;
const CCWERRBOUND_C: f64 = (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
const O3DERRBOUND_A: f64 = (7.0 + 56.0 * EPSILON) * EPSILON;
const O3DERRBOUND_B: f64 = (3.0 + 28.0 * EPSILON) * EPSILON;
const O3DERRBOUND_C: f64 = (26.0 + 288.0 * EPSILON) * EPSILON * EPSILON;
const ICCERRBOUND_A: f64 = (10.0 + 96.0 * EPSILON) * EPSILON;
const ICCERRBOUND_B: f64 = (4.0 + 48.0 * EPSILON) * EPSILON;
const ICCERRBOUND_C: f64 = (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;

/// Returns a positive value if the coordinates `pa`, `pb`, and `pc` occur in counterclockwise order
/// (pc lies to the **left** of the directed line defined by coordinates pa and pb).  
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

    let detsum = if detleft > 0.0 {
        if detright <= 0.0 {
            return det;
        } else {
            detleft + detright
        }
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return det;
        } else {
            -detleft - detright
        }
    } else {
        return det;
    };
    let errbound = CCWERRBOUND_A * detsum;
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
/// (pc lies to the **left** of the directed line defined by coordinates pa and pb).  
/// Returns a negative value if `pd` lies above the plane
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
    if ((det > errbound) || (-det > errbound)) {
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
    if ((det >= errbound) || (-det >= errbound)) {
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
    if ((det >= errbound) || (-det >= errbound)) {
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
        if (adytail == 0.0) {
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
    } else {
        if adytail == 0.0 {
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
    } else {
        if bdytail == 0.0 {
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
    } else {
        if cdytail == 0.0 {
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
    }

    let mut bct = [0f64; 8];
    let mut cat = [0f64; 8];
    let mut abt = [0f64; 8];
    let mut u = [0f64; 4];
    let mut v = [0f64; 12];
    let mut w = [0f64; 16];
    let mut vlength: usize;
    let wlength: usize;

    let bctlen = fast_expansion_sum_zeroelim(&bt_c[..bt_clen], &ct_b[..ct_blen], &mut bct);
    let mut wlength = scale_expansion_zeroelim(&bct[..bctlen], adz, &mut w);
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

    if (adztail != 0.0) {
        vlength = scale_expansion_zeroelim(&bc[..4], adztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if (bdztail != 0.0) {
        vlength = scale_expansion_zeroelim(&ca[..4], bdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }
    if (cdztail != 0.0) {
        vlength = scale_expansion_zeroelim(&ab[..4], cdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &v[..vlength], &mut finother);
        ::core::mem::swap(&mut finnow, &mut finother);
    }

    if (adxtail != 0.0) {
        if (bdytail != 0.0) {
            let (adxt_bdyt1, adxt_bdyt0) = two_product(adxtail, bdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (cdztail != 0.0) {
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if (cdytail != 0.0) {
            let negate = -adxtail;
            let (adxt_cdyt1, adxt_cdyt0) = two_product(negate, cdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (bdztail != 0.0) {
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
    }
    if (bdxtail != 0.0) {
        if (cdytail != 0.0) {
            let (bdxt_cdyt1, bdxt_cdyt0) = two_product(bdxtail, cdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (adztail != 0.0) {
                (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if (adytail != 0.0) {
            let negate = -bdxtail;
            let (bdxt_adyt1, bdxt_adyt0) = two_product(negate, adytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_adyt1, bdxt_adyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (cdztail != 0.0) {
                (u[3], u[2], u[1], u[0]) = two_one_product(bdxt_adyt1, bdxt_adyt0, cdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
    }
    if (cdxtail != 0.0) {
        if (adytail != 0.0) {
            let (cdxt_adyt1, cdxt_adyt0) = two_product(cdxtail, adytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_adyt1, cdxt_adyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (bdztail != 0.0) {
                (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_adyt1, cdxt_adyt0, bdztail);
                finlength =
                    fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
                ::core::mem::swap(&mut finnow, &mut finother);
            }
        }
        if (bdytail != 0.0) {
            let negate = -cdxtail;
            let (cdxt_bdyt1, cdxt_bdyt0) = two_product(negate, bdytail);
            (u[3], u[2], u[1], u[0]) = two_one_product(cdxt_bdyt1, cdxt_bdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&finnow[..finlength], &u[..4], &mut finother);
            ::core::mem::swap(&mut finnow, &mut finother);
            if (adztail != 0.0) {
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

    return finnow[finlength - 1];
}

/// Returns a positive value if the coordinate `pd` lies **outside** the circle passing through `pa`, `pb`, and `pc`.  
/// Returns a negative value if it lies **inside** the circle.  
/// Returns `0` if the four points are **cocircular**.
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

// f64::abs is not yet available in libcore, so we’ll need to use a separate
// crate for this functionality.
//
// https://github.com/rust-lang/rust/issues/50145
#[cfg(feature = "no_std")]
#[inline]
fn abs(x: f64) -> f64 {
    ieee754::Ieee754::abs(x)
}
#[cfg(not(feature = "no_std"))]
#[inline]
fn abs(x: f64) -> f64 {
    x.abs()
}

#[cfg(test)]
mod test {
    use super::{incircle, orient2d, orient3d, Coord, Coord3D};

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
        let from = Coord3D { x: -1f64, y: -1.0 };
        let to = Coord3D { x: 1f64, y: 1.0 };
        let p1 = Coord3D {
            x: ::core::f64::MIN_POSITIVE,
            y: ::core::f64::MIN_POSITIVE,
        };
        let p2 = Coord3D {
            x: -::core::f64::MIN_POSITIVE,
            y: -::core::f64::MIN_POSITIVE,
        };
        let p3 = Coord3D {
            x: -::core::f64::MIN_POSITIVE,
            y: ::core::f64::MIN_POSITIVE,
        };
        let p4 = Coord3D {
            x: ::core::f64::MIN_POSITIVE,
            y: -::core::f64::MIN_POSITIVE,
        };

        // TODO: make this test work
        for &(p, sign) in &[(p1, 0.0), (p2, 0.0), (p3, 1.0), (p4, -1.0)] {
            let det = orient2d(from, to, p);
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
}
