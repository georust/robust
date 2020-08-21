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

extern crate num_traits;
use num_traits::Float;

pub trait CoordinateType: Float {
    fn two() -> Self {
        Self::one() + Self::one()
    }
}
impl<T: Float> CoordinateType for T {}

/// A two dimensional coordinate.
pub trait Coord<T: CoordinateType> {
    fn x(&self) -> T;
    fn y(&self) -> T;
}

// These values are precomputed from the "exactinit" method of the c-source code. They should? be
// the same in all IEEE-754 environments, including rust f64
const SPLITTER: f64 = 134_217_729f64;
const EPSILON: f64 = 0.000_000_000_000_000_111_022_302_462_515_65;
const RESULTERRBOUND: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;
const CCWERRBOUND_B: f64 = (2.0 + 12.0 * EPSILON) * EPSILON;
const CCWERRBOUND_C: f64 = (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
const ICCERRBOUND_A: f64 = (10.0 + 96.0 * EPSILON) * EPSILON;
const ICCERRBOUND_B: f64 = (4.0 + 48.0 * EPSILON) * EPSILON;
const ICCERRBOUND_C: f64 = (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;

/// Returns a positive value if the coordinates `pa`, `pb`, and `pc` occur in counterclockwise order
/// (pc lies to the **left** of the directed line defined by coordinates pa and pb).  
/// Returns a negative value if they occur in clockwise order (`pc` lies to the **right** of the directed line `pa, pb`).  
/// Returns `0` if they are **collinear**.  
pub fn orient2d<T: CoordinateType, C: Coord<T>>(pa: C, pb: C, pc: C) -> T {
    let detleft = (pa.x() - pc.x()) * (pb.y() - pc.y());
    let detright = (pa.y() - pc.y()) * (pb.x() - pc.x());
    let det = detleft - detright;

    let detsum = if detleft > T::zero() {
        if detright <= T::zero() {
            return det;
        } else {
            detleft + detright
        }
    } else if detleft < T::zero() {
        if detright >= T::zero() {
            return det;
        } else {
            -detleft - detright
        }
    } else {
        return det;
    };
    let errbound = T::from(CCWERRBOUND_A).unwrap() * detsum;
    if det >= errbound || -det >= errbound {
        det
    } else {
        orient2dadapt(pa, pb, pc, detsum)
    }
}

fn orient2dadapt<T: CoordinateType, C: Coord<T>>(pa: C, pb: C, pc: C, detsum: T) -> T {
    let acx = pa.x() - pc.x();
    let bcx = pb.x() - pc.x();
    let acy = pa.y() - pc.y();
    let bcy = pb.y() - pc.y();

    let (detleft, detlefttail) = two_product(acx, bcy);
    let (detright, detrighttail) = two_product(acy, bcx);

    let (B3, B2, B1, B0) = two_two_diff(detleft, detlefttail, detright, detrighttail);
    let B = [B0, B1, B2, B3];

    let mut det = estimate(&B);
    let errbound = T::from(CCWERRBOUND_B).unwrap() * detsum;
    if det >= errbound || (-det >= errbound) {
        return det;
    }

    let acxtail = two_diff_tail(pa.x(), pc.x(), acx);
    let bcxtail = two_diff_tail(pb.x(), pc.x(), bcx);
    let acytail = two_diff_tail(pa.y(), pc.y(), acy);
    let bcytail = two_diff_tail(pb.y(), pc.y(), bcy);

    if acxtail == T::zero() && acytail == T::zero() && bcxtail == T::zero() && bcytail == T::zero()
    {
        return det;
    }

    let errbound =
        T::from(CCWERRBOUND_C).unwrap() * detsum + T::from(RESULTERRBOUND).unwrap() * det.abs();
    det = det + (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);

    if det >= errbound || -det >= errbound {
        return det;
    }

    let (s1, s0) = two_product(acxtail, bcy);
    let (t1, t0) = two_product(acytail, bcx);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];

    let mut C1 = [T::zero(); 8];
    let c1length = fast_expansion_sum_zeroelim(&B, &U, &mut C1);

    let (s1, s0) = two_product(acx, bcytail);
    let (t1, t0) = two_product(acy, bcxtail);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];

    let mut C2 = [T::zero(); 12];
    let c2length = fast_expansion_sum_zeroelim(&C1[..c1length], &U, &mut C2);

    let (s1, s0) = two_product(acxtail, bcytail);
    let (t1, t0) = two_product(acytail, bcxtail);
    let (u3, u2, u1, u0) = two_two_diff(s1, s0, t1, t0);
    let U = [u0, u1, u2, u3];
    let mut D = [T::zero(); 16];
    let dlength = fast_expansion_sum_zeroelim(&C2[..c2length], &U, &mut D);
    D[dlength - 1]
}

/// Returns a positive value if the coordinate `pd` lies **outside** the circle passing through `pa`, `pb`, and `pc`.  
/// Returns a negative value if it lies **inside** the circle.  
/// Returns `0` if the four points are **cocircular**.
pub fn incircle<T: CoordinateType, C: Coord<T>>(pa: C, pb: C, pc: C, pd: C) -> T {
    let adx = pa.x() - pd.x();
    let bdx = pb.x() - pd.x();
    let cdx = pc.x() - pd.x();
    let ady = pa.y() - pd.y();
    let bdy = pb.y() - pd.y();
    let cdy = pc.y() - pd.y();

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

    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;
    let errbound = T::from(ICCERRBOUND_A).unwrap() * permanent;
    if det > errbound || -det > errbound {
        return det;
    }
    incircleadapt(pa, pb, pc, pd, permanent)
}

fn incircleadapt<T: CoordinateType, C: Coord<T>>(pa: C, pb: C, pc: C, pd: C, permanent: T) -> T {
    let mut temp8 = [T::zero(); 8];
    let mut temp16a = [T::zero(); 16];
    let mut temp16b = [T::zero(); 16];
    let mut temp16c = [T::zero(); 16];
    let mut temp32a = [T::zero(); 32];
    let mut temp32b = [T::zero(); 32];
    let mut temp48 = [T::zero(); 48];
    let mut temp64 = [T::zero(); 64];

    let adx = pa.x() - pd.x();
    let bdx = pb.x() - pd.x();
    let cdx = pc.x() - pd.x();
    let ady = pa.y() - pd.y();
    let bdy = pb.y() - pd.y();
    let cdy = pc.y() - pd.y();

    let (bdxcdy1, bdxcdy0) = two_product(bdx, cdy);
    let (cdxbdy1, cdxbdy0) = two_product(cdx, bdy);
    let (bc3, bc2, bc1, bc0) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let bc = [bc0, bc1, bc2, bc3];

    let mut axbc = [T::zero(); 8];
    let axbclen = scale_expansion_zeroelim(&bc, adx, &mut axbc);
    let mut axxbc = [T::zero(); 16];
    let axxbclen = scale_expansion_zeroelim(&axbc[..axbclen], adx, &mut axxbc);
    let mut aybc = [T::zero(); 8];
    let aybclen = scale_expansion_zeroelim(&bc, ady, &mut aybc);
    let mut ayybc = [T::zero(); 16];
    let ayybclen = scale_expansion_zeroelim(&aybc[..aybclen], ady, &mut ayybc);
    let mut adet = [T::zero(); 32];
    let alen = fast_expansion_sum_zeroelim(&axxbc[0..axxbclen], &ayybc[0..ayybclen], &mut adet);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let (c3, c2, c1, c0) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let ca = [c0, c1, c2, c3];

    let mut bxca = [T::zero(); 8];
    let bxcalen = scale_expansion_zeroelim(&ca, bdx, &mut bxca);
    let mut bxxca = [T::zero(); 16];
    let bxxcalen = scale_expansion_zeroelim(&bxca[..bxcalen], bdx, &mut bxxca);
    let mut byca = [T::zero(); 8];
    let bycalen = scale_expansion_zeroelim(&ca, bdy, &mut byca);
    let mut byyca = [T::zero(); 16];
    let byycalen = scale_expansion_zeroelim(&byca[..bycalen], bdy, &mut byyca);
    let mut bdet = [T::zero(); 32];
    let blen = fast_expansion_sum_zeroelim(&bxxca[..bxxcalen], &byyca[0..byycalen], &mut bdet);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let (ab3, ab2, ab1, ab0) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let ab = [ab0, ab1, ab2, ab3];

    let mut cxab = [T::zero(); 8];
    let cxablen = scale_expansion_zeroelim(&ab, cdx, &mut cxab);
    let mut cxxab = [T::zero(); 16];
    let cxxablen = scale_expansion_zeroelim(&cxab[..cxablen], cdx, &mut cxxab);
    let mut cyab = [T::zero(); 8];
    let cyablen = scale_expansion_zeroelim(&ab, cdy, &mut cyab);
    let mut cyyab = [T::zero(); 16];
    let cyyablen = scale_expansion_zeroelim(&cyab[..cyablen], cdy, &mut cyyab);
    let mut cdet = [T::zero(); 32];
    let clen = fast_expansion_sum_zeroelim(&cxxab[..cxxablen], &cyyab[..cyyablen], &mut cdet);

    let mut abdet = [T::zero(); 64];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut fin1 = [T::zero(); 1152];
    let mut finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut fin1);

    let mut det = estimate(&fin1[..finlength]);
    let errbound = T::from(ICCERRBOUND_B).unwrap() * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }

    let adxtail = two_diff_tail(pa.x(), pd.x(), adx);
    let adytail = two_diff_tail(pa.y(), pd.y(), ady);
    let bdxtail = two_diff_tail(pb.x(), pd.x(), bdx);
    let bdytail = two_diff_tail(pb.y(), pd.y(), bdy);
    let cdxtail = two_diff_tail(pc.x(), pd.x(), cdx);
    let cdytail = two_diff_tail(pc.y(), pd.y(), cdy);
    if adxtail == T::zero()
        && bdxtail == T::zero()
        && cdxtail == T::zero()
        && adytail == T::zero()
        && bdytail == T::zero()
        && cdytail == T::zero()
    {
        return det;
    }

    let errbound =
        T::from(ICCERRBOUND_C).unwrap() * permanent + T::from(RESULTERRBOUND).unwrap() * det.abs();
    det = det
        + ((adx * adx + ady * ady)
            * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
            + T::two() * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
        + ((bdx * bdx + bdy * bdy)
            * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + T::two() * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
        + ((cdx * cdx + cdy * cdy)
            * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + T::two() * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));

    if det >= errbound || -det >= errbound {
        return det;
    }

    let mut fin2 = [T::zero(); 1152];

    let aa = if bdxtail != T::zero()
        || bdytail != T::zero()
        || cdxtail != T::zero()
        || cdytail != T::zero()
    {
        let (adxadx1, adxadx0) = square(adx);
        let (adyady1, adyady0) = square(ady);
        let (aa3, aa2, aa1, aa0) = two_two_sum(adxadx1, adxadx0, adyady1, adyady0);
        [aa0, aa1, aa2, aa3]
    } else {
        [T::zero(); 4]
    };

    let bb = if cdxtail != T::zero()
        || cdytail != T::zero()
        || adxtail != T::zero()
        || adytail != T::zero()
    {
        let (bdxbdx1, bdxbdx0) = square(bdx);
        let (bdybdy1, bdybdy0) = square(bdy);
        let (bb3, bb2, bb1, bb0) = two_two_sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0);
        [bb0, bb1, bb2, bb3]
    } else {
        [T::zero(); 4]
    };

    let cc = if adxtail != T::zero()
        || adytail != T::zero()
        || bdxtail != T::zero()
        || bdytail != T::zero()
    {
        let (cdxcdx1, cdxcdx0) = square(cdx);
        let (cdycdy1, cdycdy0) = square(cdy);
        let (cc3, cc2, cc1, cc0) = two_two_sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0);
        [cc0, cc1, cc2, cc3]
    } else {
        [T::zero(); 4]
    };

    let mut axtbclen = 9;
    let mut axtbc = [T::zero(); 8];
    if adxtail != T::zero() {
        axtbclen = scale_expansion_zeroelim(&bc, adxtail, &mut axtbc);
        let mut temp16a = [T::zero(); 16];
        let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], T::two() * adx, &mut temp16a);

        let mut axtcc = [T::zero(); 8];
        let axtcclen = scale_expansion_zeroelim(&cc, adxtail, &mut axtcc);
        let mut temp16b = [T::zero(); 16];
        let temp16blen = scale_expansion_zeroelim(&axtcc[..axtcclen], bdy, &mut temp16b);

        let mut axtbb = [T::zero(); 8];
        let axtbblen = scale_expansion_zeroelim(&bb, adxtail, &mut axtbb);
        let mut temp16c = [T::zero(); 16];
        let temp16clen = scale_expansion_zeroelim(&axtbb[..axtbblen], -cdy, &mut temp16c);

        let mut temp32a = [T::zero(); 32];
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let mut temp48 = [T::zero(); 48];
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
        ::std::mem::swap(&mut fin1, &mut fin2)
    }

    let mut aytbclen = 9;
    let mut aytbc = [T::zero(); 8];
    if adytail != T::zero() {
        aytbclen = scale_expansion_zeroelim(&bc, adytail, &mut aytbc);
        let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], T::two() * ady, &mut temp16a);

        let mut aytbb = [T::zero(); 8];
        let aytbblen = scale_expansion_zeroelim(&bb, adytail, &mut aytbb);
        let temp16blen = scale_expansion_zeroelim(&aytbb[..aytbblen], cdx, &mut temp16b);

        let mut aytcc = [T::zero(); 8];
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
        ::std::mem::swap(&mut fin1, &mut fin2)
    }

    let mut bxtcalen = 9;
    let mut bxtca = [T::zero(); 8];
    if bdxtail != T::zero() {
        bxtcalen = scale_expansion_zeroelim(&ca, bdxtail, &mut bxtca);
        let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], T::two() * bdx, &mut temp16a);

        let mut bxtaa = [T::zero(); 8];
        let bxtaalen = scale_expansion_zeroelim(&aa, bdxtail, &mut bxtaa);
        let temp16blen = scale_expansion_zeroelim(&bxtaa[..bxtaalen], cdy, &mut temp16b);

        let mut bxtcc = [T::zero(); 8];
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
        ::std::mem::swap(&mut fin1, &mut fin2)
    }

    let mut bytcalen = 9;
    let mut bytca = [T::zero(); 8];
    if bdytail != T::zero() {
        bytcalen = scale_expansion_zeroelim(&ca, bdytail, &mut bytca);
        let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], T::two() * bdy, &mut temp16a);

        let mut bytcc = [T::zero(); 8];
        let bytcclen = scale_expansion_zeroelim(&cc, bdytail, &mut bytcc);
        let temp16blen = scale_expansion_zeroelim(&bytcc[..bytcclen], adx, &mut temp16b);

        let mut bytaa = [T::zero(); 8];
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
        ::std::mem::swap(&mut fin1, &mut fin2)
    }

    let mut cxtab = [T::zero(); 8];
    let mut cxtablen = 9;
    if cdxtail != T::zero() {
        cxtablen = scale_expansion_zeroelim(&ab, cdxtail, &mut cxtab);
        let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], T::two() * cdx, &mut temp16a);

        let mut cxtbb = [T::zero(); 8];
        let cxtbblen = scale_expansion_zeroelim(&bb, cdxtail, &mut cxtbb);
        let temp16blen = scale_expansion_zeroelim(&cxtbb[..cxtbblen], ady, &mut temp16b);

        let mut cxtaa = [T::zero(); 8];
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
        ::std::mem::swap(&mut fin1, &mut fin2);
    }

    let mut cytab = [T::zero(); 8];
    let mut cytablen = 9;
    if cdytail != T::zero() {
        cytablen = scale_expansion_zeroelim(&ab, cdytail, &mut cytab);
        let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], T::two() * cdy, &mut temp16a);

        let mut cytaa = [T::zero(); 8];
        let cytaalen = scale_expansion_zeroelim(&aa, cdytail, &mut cytaa);
        let temp16blen = scale_expansion_zeroelim(&cytaa[..cytaalen], bdx, &mut temp16b);

        let mut cytbb = [T::zero(); 8];
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
        ::std::mem::swap(&mut fin1, &mut fin2);
    }

    if adxtail != T::zero() || adytail != T::zero() {
        let mut bctt = [T::zero(); 4];
        let mut bct = [T::zero(); 8];
        let bcttlen;
        let bctlen;
        if bdxtail != T::zero()
            || bdytail != T::zero()
            || cdxtail != T::zero()
            || cdytail != T::zero()
        {
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
            bct[0] = T::zero();
            bctlen = 1;
            bctt[0] = T::zero();
            bcttlen = 1;
        }

        if adxtail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], adxtail, &mut temp16a);
            let mut axtbct = [T::zero(); 16];
            let axtbctlen = scale_expansion_zeroelim(&bct[..bctlen], adxtail, &mut axtbct);
            let temp32alen =
                scale_expansion_zeroelim(&axtbct[..axtbctlen], T::two() * adx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            if bdytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&cc, adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }
            if cdytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&bb, -adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&axtbct[..axtbctlen], adxtail, &mut temp32a);
            let mut axtbctt = [T::zero(); 8];
            let axtbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adxtail, &mut axtbctt);
            let temp16alen =
                scale_expansion_zeroelim(&axtbctt[..axtbcttlen], T::two() * adx, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }

        if adytail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], adytail, &mut temp16a);
            let mut aytbct = [T::zero(); 16];
            let aytbctlen = scale_expansion_zeroelim(&bct[..bctlen], adytail, &mut aytbct);
            let temp32alen =
                scale_expansion_zeroelim(&aytbct[..aytbctlen], T::two() * ady, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&aytbct[..aytbctlen], adytail, &mut temp32a);
            let mut aytbctt = [T::zero(); 8];
            let aytbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adytail, &mut aytbctt);
            let temp16alen =
                scale_expansion_zeroelim(&aytbctt[..aytbcttlen], T::two() * ady, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }
    }

    if bdxtail != T::zero() || bdytail != T::zero() {
        let mut catt = [T::zero(); 4];
        let mut cat = [T::zero(); 8];
        let cattlen;
        let catlen;

        if cdxtail != T::zero()
            || cdytail != T::zero()
            || adxtail != T::zero()
            || adytail != T::zero()
        {
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
            cat[0] = T::zero();
            catlen = 1;
            catt[0] = T::zero();
            cattlen = 1;
        }

        if bdxtail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], bdxtail, &mut temp16a);
            let mut bxtcat = [T::zero(); 16];
            let bxtcatlen = scale_expansion_zeroelim(&cat[..catlen], bdxtail, &mut bxtcat);
            let temp32alen =
                scale_expansion_zeroelim(&bxtcat[..bxtcatlen], T::two() * bdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            if cdytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&aa, bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }
            if adytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&cc, -bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&bxtcat[..bxtcatlen], bdxtail, &mut temp32a);
            let mut bxtcatt = [T::zero(); 8];
            let bxtcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdxtail, &mut bxtcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bxtcatt[..bxtcattlen], T::two() * bdx, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }
        if bdytail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], bdytail, &mut temp16a);
            let mut bytcat = [T::zero(); 16];
            let bytcatlen = scale_expansion_zeroelim(&cat[..catlen], bdytail, &mut bytcat);
            let temp32alen =
                scale_expansion_zeroelim(&bytcat[..bytcatlen], T::two() * bdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&bytcat[..bytcatlen], bdytail, &mut temp32a);
            let mut bytcatt = [T::zero(); 8];
            let bytcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdytail, &mut bytcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bytcatt[..bytcattlen], T::two() * bdy, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }
    }

    if cdxtail != T::zero() || cdytail != T::zero() {
        let mut abtt = [T::zero(); 4];
        let mut abt = [T::zero(); 8];
        let abttlen;
        let abtlen;

        if adxtail != T::zero()
            || adytail != T::zero()
            || bdxtail != T::zero()
            || bdytail != T::zero()
        {
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
            abt[0] = T::zero();
            abtlen = 1;
            abtt[0] = T::zero();
            abttlen = 1;
        }

        if cdxtail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], cdxtail, &mut temp16a);
            let mut cxtabt = [T::zero(); 16];
            let cxtabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdxtail, &mut cxtabt);
            let temp32alen =
                scale_expansion_zeroelim(&cxtabt[..cxtabtlen], T::two() * cdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            if adytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&bb, cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }
            if bdytail != T::zero() {
                let temp8len = scale_expansion_zeroelim(&aa, -cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &fin1[..finlength],
                    &temp16a[..temp16alen],
                    &mut fin2,
                );
                ::std::mem::swap(&mut fin1, &mut fin2);
            }

            let temp32alen = scale_expansion_zeroelim(&cxtabt[..cxtabtlen], cdxtail, &mut temp32a);
            let mut cxtabtt = [T::zero(); 8];
            let cxtabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdxtail, &mut cxtabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cxtabtt[..cxtabttlen], T::two() * cdx, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }
        if cdytail != T::zero() {
            let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], cdytail, &mut temp16a);
            let mut cytabt = [T::zero(); 16];
            let cytabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdytail, &mut cytabt);
            let temp32alen =
                scale_expansion_zeroelim(&cytabt[..cytabtlen], T::two() * cdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&fin1[..finlength], &temp48[..temp48len], &mut fin2);
            ::std::mem::swap(&mut fin1, &mut fin2);

            let temp32alen = scale_expansion_zeroelim(&cytabt[..cytabtlen], cdytail, &mut temp32a);
            let mut cytabtt = [T::zero(); 8];
            let cytabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdytail, &mut cytabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cytabtt[..cytabttlen], T::two() * cdy, &mut temp16a);
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
            ::std::mem::swap(&mut fin1, &mut fin2);
        }
    }
    fin1[finlength - 1]
}

fn scale_expansion_zeroelim<T: CoordinateType>(e: &[T], b: T, h: &mut [T]) -> usize {
    let (bhi, blo) = split(b);
    let (mut Q, hh) = two_product_presplit(e[0], b, bhi, blo);
    let mut hindex = 0;
    if hh != T::zero() {
        h[hindex] = hh;
        hindex += 1;
    }
    for eindex in 1..e.len() {
        let enow = e[eindex];
        let (product1, product0) = two_product_presplit(enow, b, bhi, blo);
        let (sum, hh) = two_sum(Q, product0);
        if hh != T::zero() {
            h[hindex] = hh;
            hindex += 1;
        }
        let (new_q, hh) = fast_two_sum(product1, sum);
        Q = new_q;
        if hh != T::zero() {
            h[hindex] = hh;
            hindex += 1;
        }
    }
    if Q != T::zero() || hindex == 0 {
        h[hindex] = Q;
        hindex += 1;
    }
    hindex
}

#[inline]
fn two_product<T: CoordinateType>(a: T, b: T) -> (T, T) {
    let x = a * b;
    (x, two_product_tail(a, b, x))
}

#[inline]
fn two_product_tail<T: CoordinateType>(a: T, b: T, x: T) -> T {
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    (alo * blo) - err3
}

#[inline]
fn split<T: CoordinateType>(a: T) -> (T, T) {
    let c = T::from(SPLITTER).unwrap() * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;
    (ahi, alo)
}

#[inline]
fn two_product_presplit<T: CoordinateType>(a: T, b: T, bhi: T, blo: T) -> (T, T) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let err1 = x - ahi * bhi;
    let err2 = err1 - alo * bhi;
    let err3 = err2 - ahi * blo;
    let y = alo * blo - err3;
    (x, y)
}

#[inline]
fn two_two_diff<T: CoordinateType>(a1: T, a0: T, b1: T, b0: T) -> (T, T, T, T) {
    let (j, _r0, x0) = two_one_diff(a1, a0, b0);
    let (x3, x2, x1) = two_one_diff(j, _r0, b1);
    (x3, x2, x1, x0)
}

#[inline]
fn two_one_diff<T: CoordinateType>(a1: T, a0: T, b: T) -> (T, T, T) {
    let (i, x0) = two_diff(a0, b);
    let (x2, x1) = two_sum(a1, i);
    (x2, x1, x0)
}

#[inline]
fn two_diff<T: CoordinateType>(a: T, b: T) -> (T, T) {
    let x = a - b;
    (x, two_diff_tail(a, b, x))
}

#[inline]
fn two_diff_tail<T: CoordinateType>(a: T, b: T, x: T) -> T {
    let bvirt = a - x;
    let avirt = x + bvirt;
    let bround = bvirt - b;
    let around = a - avirt;
    around + bround
}

#[inline]
fn two_sum<T: CoordinateType>(a: T, b: T) -> (T, T) {
    let x = a + b;
    (x, two_sum_tail(a, b, x))
}

#[inline]
fn two_sum_tail<T: CoordinateType>(a: T, b: T, x: T) -> T {
    let bvirt = x - a;
    let avirt = x - bvirt;
    let bround = b - bvirt;
    let around = a - avirt;
    around + bround
}

fn estimate<T: CoordinateType>(e: &[T]) -> T {
    let mut q = e[0];
    for cur in &e[1..] {
        q = q + *cur;
    }
    q
}

fn fast_expansion_sum_zeroelim<T: CoordinateType>(e: &[T], f: &[T], h: &mut [T]) -> usize {
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
        if hh != T::zero() {
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
            if hh != T::zero() {
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
        if hh != T::zero() {
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
        if hh != T::zero() {
            h[hindex] = hh;
            hindex += 1
        }
    }

    if Q != T::zero() || hindex == 0 {
        h[hindex] = Q;
        hindex += 1;
    }
    hindex
}

#[inline]
fn fast_two_sum_tail<T: CoordinateType>(a: T, b: T, x: T) -> T {
    let bvirt = x - a;
    b - bvirt
}

#[inline]
fn fast_two_sum<T: CoordinateType>(a: T, b: T) -> (T, T) {
    let x = a + b;
    (x, fast_two_sum_tail(a, b, x))
}

#[inline]
fn square_tail<T: CoordinateType>(a: T, x: T) -> T {
    let (ahi, alo) = split(a);
    let err1 = x - ahi * ahi;
    let err3 = err1 - (ahi + ahi) * alo;
    alo * alo - err3
}

#[inline]
fn square<T: CoordinateType>(a: T) -> (T, T) {
    let x = a * a;
    (x, square_tail(a, x))
}

#[inline]
fn two_one_sum<T: CoordinateType>(a1: T, a0: T, b: T) -> (T, T, T) {
    let (_i, x0) = two_sum(a0, b);
    let (x2, x1) = two_sum(a1, _i);
    (x2, x1, x0)
}

#[inline]
fn two_two_sum<T: CoordinateType>(a1: T, a0: T, b1: T, b0: T) -> (T, T, T, T) {
    let (_j, _r0, x0) = two_one_sum(a1, a0, b0);
    let (x3, x2, x1) = two_one_sum(_j, _r0, b1);
    (x3, x2, x1, x0)
}

#[cfg(test)]
mod test {
    use super::{incircle, orient2d, Coord};

    #[derive(Clone, Copy)]
    struct MyCoord {
        x: f64,
        y: f64,
    }

    impl Coord<f64> for MyCoord {
        fn x(&self) -> f64 {
            self.x
        }

        fn y(&self) -> f64 {
            self.y
        }
    }

    #[test]
    fn test_orient2d() {
        let from = MyCoord { x: -1f64, y: -1.0 };
        let to = MyCoord { x: 1f64, y: 1.0 };
        let p1 = MyCoord {
            x: ::std::f64::MIN_POSITIVE,
            y: ::std::f64::MIN_POSITIVE,
        };
        let p2 = MyCoord {
            x: -::std::f64::MIN_POSITIVE,
            y: -::std::f64::MIN_POSITIVE,
        };
        let p3 = MyCoord {
            x: -::std::f64::MIN_POSITIVE,
            y: ::std::f64::MIN_POSITIVE,
        };
        let p4 = MyCoord {
            x: ::std::f64::MIN_POSITIVE,
            y: -::std::f64::MIN_POSITIVE,
        };

        for &(p, sign) in &[(p1, 0.0), (p2, 0.0), (p3, 1.0), (p4, -1.0)] {
            let det = orient2d(from, to, p);
            assert!(det == sign || det.signum() == sign.signum());
        }
    }

    #[test]
    fn test_incircle() {
        let from = MyCoord { x: -1f64, y: -1.0 };
        let to = MyCoord { x: 1f64, y: 1.0 };
        let p_left = MyCoord {
            x: -::std::f64::MIN_POSITIVE,
            y: ::std::f64::MIN_POSITIVE,
        };
        let p_right = MyCoord {
            x: ::std::f64::MIN_POSITIVE,
            y: -::std::f64::MIN_POSITIVE,
        };
        let p_query = MyCoord { x: 2.0, y: 2.0 };

        assert!(incircle(from, p_left, to, p_query) > 0.0);
        assert!(incircle(from, to, p_right, p_query) > 0.0);
    }

    #[test]
    fn test_issue48_a() {
        let pa = MyCoord {
            x: 2.1045541600524288e-15,
            y: -1.0000000000000016,
        };
        let pb = MyCoord {
            x: 1.000000000000005,
            y: -3.350874324301223e-16,
        };
        let pc = MyCoord {
            x: 7.553997323229233e-15,
            y: 0.9999999999999958,
        };
        let pd = MyCoord {
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
        let pa = MyCoord {
            x: 9.128561612013288e-15,
            y: -1.0000000000000029,
        };
        let pb = MyCoord {
            x: 1.0000000000000044,
            y: -5.451395142523081e-15,
        };
        let pc = MyCoord {
            x: 3.851214418148064e-15,
            y: 0.9999999999999961,
        };
        let pd = MyCoord {
            x: -0.9999999999999946,
            y: -6.6797960341085084e-15,
        };
        // the (incorrect) result from the previous version of exactpred.rs
        assert!(incircle(pa, pb, pc, pd) != -1.1074731814540733e-16);
        // the result predicates.c gives
        assert!(incircle(pa, pb, pc, pd) == 7.226864249343135e-30);
    }
}
