use ::libc;

use crate::bezctx;
unsafe extern "C" {
    fn atan2(_: libc::c_double, _: libc::c_double) -> libc::c_double;
    fn cos(_: libc::c_double) -> libc::c_double;
    fn sin(_: libc::c_double) -> libc::c_double;
    fn fabs(_: libc::c_double) -> libc::c_double;
    fn hypot(_: libc::c_double, _: libc::c_double) -> libc::c_double;
    fn floor(_: libc::c_double) -> libc::c_double;
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn memcpy(_: *mut libc::c_void, _: *const libc::c_void, _: libc::c_ulong) -> *mut libc::c_void;
    fn bezctx_moveto(
        bc: *mut bezctx,
        x: libc::c_double,
        y: libc::c_double,
        is_open: libc::c_int,
        si: libc::c_int,
    );
    fn bezctx_lineto(bc: *mut bezctx, x: libc::c_double, y: libc::c_double, si: libc::c_int);
    fn bezctx_quadto(
        bc: *mut bezctx,
        x1: libc::c_double,
        y1: libc::c_double,
        x2: libc::c_double,
        y2: libc::c_double,
        si: libc::c_int,
    );
    fn bezctx_curveto(
        bc: *mut bezctx,
        x1: libc::c_double,
        y1: libc::c_double,
        x2: libc::c_double,
        y2: libc::c_double,
        x3: libc::c_double,
        y3: libc::c_double,
        si: libc::c_int,
    );
    fn bezctx_mark_knot(bc: *mut bezctx, knot_idx: libc::c_int, i: libc::c_int);
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct spiro_cp {
    pub x: libc::c_double,
    pub y: libc::c_double,
    pub ty: libc::c_char,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct spiro_seg_s {
    pub x: libc::c_double,
    pub y: libc::c_double,
    pub ty: libc::c_char,
    pub bend_th: libc::c_double,
    pub ks: [libc::c_double; 4],
    pub seg_ch: libc::c_double,
    pub seg_th: libc::c_double,
    pub l: libc::c_double,
}
pub type spiro_seg = spiro_seg_s;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct bandmat {
    pub a: [libc::c_double; 11],
    pub al: [libc::c_double; 5],
}
#[inline(always)]
unsafe extern "C" fn __inline_isfinitef(mut __x: libc::c_float) -> libc::c_int {
    return (__x == __x && __x.abs() != ::core::f32::INFINITY) as libc::c_int;
}
#[inline(always)]
unsafe extern "C" fn __inline_isfinited(mut __x: libc::c_double) -> libc::c_int {
    return (__x == __x && __x.abs() != ::core::f64::INFINITY) as libc::c_int;
}
unsafe extern "C" fn integrate_spiro(
    mut ks: *const libc::c_double,
    mut xy: *mut libc::c_double,
    mut n: libc::c_int,
) {
    unsafe {
        let mut th1: libc::c_double = *ks.offset(0 as libc::c_int as isize);
        let mut th2: libc::c_double = 0.5f64 * *ks.offset(1 as libc::c_int as isize);
        let mut th3: libc::c_double =
            1.0f64 / 6 as libc::c_int as libc::c_double * *ks.offset(2 as libc::c_int as isize);
        let mut th4: libc::c_double =
            1.0f64 / 24 as libc::c_int as libc::c_double * *ks.offset(3 as libc::c_int as isize);
        let mut x: libc::c_double = 0.;
        let mut y: libc::c_double = 0.;
        let mut ds: libc::c_double = 1.0f64 / n as libc::c_double;
        let mut ds2: libc::c_double = ds * ds;
        let mut ds3: libc::c_double = ds2 * ds;
        let mut k0: libc::c_double = *ks.offset(0 as libc::c_int as isize) * ds;
        let mut k1: libc::c_double = *ks.offset(1 as libc::c_int as isize) * ds;
        let mut k2: libc::c_double = *ks.offset(2 as libc::c_int as isize) * ds;
        let mut k3: libc::c_double = *ks.offset(3 as libc::c_int as isize) * ds;
        let mut i: libc::c_int = 0;
        let mut s: libc::c_double = 0.5f64 * ds - 0.5f64;
        x = 0 as libc::c_int as libc::c_double;
        y = 0 as libc::c_int as libc::c_double;
        i = 0 as libc::c_int;
        while i < n {
            let mut u: libc::c_double = 0.;
            let mut v: libc::c_double = 0.;
            let mut km0: libc::c_double = 0.;
            let mut km1: libc::c_double = 0.;
            let mut km2: libc::c_double = 0.;
            let mut km3: libc::c_double = 0.;
            if n == 1 as libc::c_int {
                km0 = k0;
                km1 = k1 * ds;
                km2 = k2 * ds2;
            } else {
                km0 = ((1.0f64 / 6 as libc::c_int as libc::c_double * k3 * s + 0.5f64 * k2) * s
                    + k1)
                    * s
                    + k0;
                km1 = ((0.5f64 * k3 * s + k2) * s + k1) * ds;
                km2 = (k3 * s + k2) * ds2;
            }
            km3 = k3 * ds3;
            let mut t1_1: libc::c_double = km0;
            let mut t1_2: libc::c_double = 0.5f64 * km1;
            let mut t1_3: libc::c_double = 1.0f64 / 6 as libc::c_int as libc::c_double * km2;
            let mut t1_4: libc::c_double = 1.0f64 / 24 as libc::c_int as libc::c_double * km3;
            let mut t2_2: libc::c_double = t1_1 * t1_1;
            let mut t2_3: libc::c_double = 2 as libc::c_int as libc::c_double * (t1_1 * t1_2);
            let mut t2_4: libc::c_double =
                2 as libc::c_int as libc::c_double * (t1_1 * t1_3) + t1_2 * t1_2;
            let mut t2_5: libc::c_double =
                2 as libc::c_int as libc::c_double * (t1_1 * t1_4 + t1_2 * t1_3);
            let mut t2_6: libc::c_double =
                2 as libc::c_int as libc::c_double * (t1_2 * t1_4) + t1_3 * t1_3;
            let mut t2_7: libc::c_double = 2 as libc::c_int as libc::c_double * (t1_3 * t1_4);
            let mut t2_8: libc::c_double = t1_4 * t1_4;
            let mut t3_4: libc::c_double = t2_2 * t1_2 + t2_3 * t1_1;
            let mut t3_6: libc::c_double = t2_2 * t1_4 + t2_3 * t1_3 + t2_4 * t1_2 + t2_5 * t1_1;
            let mut t3_8: libc::c_double = t2_4 * t1_4 + t2_5 * t1_3 + t2_6 * t1_2 + t2_7 * t1_1;
            let mut t3_10: libc::c_double = t2_6 * t1_4 + t2_7 * t1_3 + t2_8 * t1_2;
            let mut t4_4: libc::c_double = t2_2 * t2_2;
            let mut t4_5: libc::c_double = 2 as libc::c_int as libc::c_double * (t2_2 * t2_3);
            let mut t4_6: libc::c_double =
                2 as libc::c_int as libc::c_double * (t2_2 * t2_4) + t2_3 * t2_3;
            let mut t4_7: libc::c_double =
                2 as libc::c_int as libc::c_double * (t2_2 * t2_5 + t2_3 * t2_4);
            let mut t4_8: libc::c_double =
                2 as libc::c_int as libc::c_double * (t2_2 * t2_6 + t2_3 * t2_5) + t2_4 * t2_4;
            let mut t4_9: libc::c_double =
                2 as libc::c_int as libc::c_double * (t2_2 * t2_7 + t2_3 * t2_6 + t2_4 * t2_5);
            let mut t4_10: libc::c_double = 2 as libc::c_int as libc::c_double
                * (t2_2 * t2_8 + t2_3 * t2_7 + t2_4 * t2_6)
                + t2_5 * t2_5;
            let mut t5_6: libc::c_double = t4_4 * t1_2 + t4_5 * t1_1;
            let mut t5_8: libc::c_double = t4_4 * t1_4 + t4_5 * t1_3 + t4_6 * t1_2 + t4_7 * t1_1;
            let mut t5_10: libc::c_double = t4_6 * t1_4 + t4_7 * t1_3 + t4_8 * t1_2 + t4_9 * t1_1;
            let mut t6_6: libc::c_double = t4_4 * t2_2;
            let mut t6_7: libc::c_double = t4_4 * t2_3 + t4_5 * t2_2;
            let mut t6_8: libc::c_double = t4_4 * t2_4 + t4_5 * t2_3 + t4_6 * t2_2;
            let mut t6_9: libc::c_double = t4_4 * t2_5 + t4_5 * t2_4 + t4_6 * t2_3 + t4_7 * t2_2;
            let mut t6_10: libc::c_double =
                t4_4 * t2_6 + t4_5 * t2_5 + t4_6 * t2_4 + t4_7 * t2_3 + t4_8 * t2_2;
            let mut t7_8: libc::c_double = t6_6 * t1_2 + t6_7 * t1_1;
            let mut t7_10: libc::c_double = t6_6 * t1_4 + t6_7 * t1_3 + t6_8 * t1_2 + t6_9 * t1_1;
            let mut t8_8: libc::c_double = t6_6 * t2_2;
            let mut t8_9: libc::c_double = t6_6 * t2_3 + t6_7 * t2_2;
            let mut t8_10: libc::c_double = t6_6 * t2_4 + t6_7 * t2_3 + t6_8 * t2_2;
            let mut t9_10: libc::c_double = t8_8 * t1_2 + t8_9 * t1_1;
            let mut t10_10: libc::c_double = t8_8 * t2_2;
            u = 1 as libc::c_int as libc::c_double;
            v = 0 as libc::c_int as libc::c_double;
            v += 1.0f64 / 12 as libc::c_int as libc::c_double * t1_2
                + 1.0f64 / 80 as libc::c_int as libc::c_double * t1_4;
            u -= 1.0f64 / 24 as libc::c_int as libc::c_double * t2_2
                + 1.0f64 / 160 as libc::c_int as libc::c_double * t2_4
                + 1.0f64 / 896 as libc::c_int as libc::c_double * t2_6
                + 1.0f64 / 4608 as libc::c_int as libc::c_double * t2_8;
            v -= 1.0f64 / 480 as libc::c_int as libc::c_double * t3_4
                + 1.0f64 / 2688 as libc::c_int as libc::c_double * t3_6
                + 1.0f64 / 13824 as libc::c_int as libc::c_double * t3_8
                + 1.0f64 / 67584 as libc::c_int as libc::c_double * t3_10;
            u += 1.0f64 / 1920 as libc::c_int as libc::c_double * t4_4
                + 1.0f64 / 10752 as libc::c_int as libc::c_double * t4_6
                + 1.0f64 / 55296 as libc::c_int as libc::c_double * t4_8
                + 1.0f64 / 270336 as libc::c_int as libc::c_double * t4_10;
            v += 1.0f64 / 53760 as libc::c_int as libc::c_double * t5_6
                + 1.0f64 / 276480 as libc::c_int as libc::c_double * t5_8
                + 1.0f64 / 1.35168e+06f64 * t5_10;
            u -= 1.0f64 / 322560 as libc::c_int as libc::c_double * t6_6
                + 1.0f64 / 1.65888e+06f64 * t6_8
                + 1.0f64 / 8.11008e+06f64 * t6_10;
            v -= 1.0f64 / 1.16122e+07f64 * t7_8 + 1.0f64 / 5.67706e+07f64 * t7_10;
            u += 1.0f64 / 9.28973e+07f64 * t8_8 + 1.0f64 / 4.54164e+08f64 * t8_10;
            v += 1.0f64 / 4.08748e+09f64 * t9_10;
            u -= 1.0f64 / 4.08748e+10f64 * t10_10;
            if n == 1 as libc::c_int {
                x = u;
                y = v;
            } else {
                let mut th: libc::c_double = (((th4 * s + th3) * s + th2) * s + th1) * s;
                let mut cth: libc::c_double = cos(th);
                let mut sth: libc::c_double = sin(th);
                x += cth * u - sth * v;
                y += cth * v + sth * u;
                s += ds;
            }
            i += 1;
            i;
        }
        *xy.offset(0 as libc::c_int as isize) = x * ds;
        *xy.offset(1 as libc::c_int as isize) = y * ds;
    }
}
unsafe extern "C" fn set_dm_to_1(mut dm: *mut libc::c_double) {
    unsafe {
        *dm.offset(0 as libc::c_int as isize) = 1.0f64;
        let ref mut fresh0 = *dm.offset(2 as libc::c_int as isize);
        *fresh0 = 0.0f64;
        *dm.offset(1 as libc::c_int as isize) = *fresh0;
    }
}
unsafe extern "C" fn set_di_to_x1y1(
    mut di: *mut libc::c_double,
    mut dm: *mut libc::c_double,
    mut x1: libc::c_double,
    mut y1: libc::c_double,
) {
    unsafe {
        let ref mut fresh1 = *di.offset(4 as libc::c_int as isize);
        *fresh1 = x1;
        *di.offset(3 as libc::c_int as isize) = *fresh1;
        let ref mut fresh2 = *di.offset(7 as libc::c_int as isize);
        *fresh2 = y1;
        *di.offset(6 as libc::c_int as isize) = *fresh2;
        *di.offset(3 as libc::c_int as isize) -= *di.offset(1 as libc::c_int as isize);
        *di.offset(4 as libc::c_int as isize) += *di.offset(1 as libc::c_int as isize);
        *di.offset(6 as libc::c_int as isize) -= *di.offset(1 as libc::c_int as isize);
        *di.offset(7 as libc::c_int as isize) += *di.offset(1 as libc::c_int as isize);
        *di.offset(2 as libc::c_int as isize) =
            x1 * *dm.offset(0 as libc::c_int as isize) + *dm.offset(1 as libc::c_int as isize);
        *di.offset(5 as libc::c_int as isize) =
            y1 * *dm.offset(0 as libc::c_int as isize) + *dm.offset(2 as libc::c_int as isize);
    }
}
unsafe extern "C" fn compute_ends(
    mut ks: *const libc::c_double,
    mut ends: *mut [libc::c_double; 4],
    mut seg_ch: libc::c_double,
) -> libc::c_double {
    unsafe {
        let mut xy: [libc::c_double; 2] = [0.; 2];
        let mut ch: libc::c_double = 0.;
        let mut th: libc::c_double = 0.;
        let mut l: libc::c_double = 0.;
        let mut l2: libc::c_double = 0.;
        let mut l3: libc::c_double = 0.;
        let mut th_even: libc::c_double = 0.;
        let mut th_odd: libc::c_double = 0.;
        let mut k0_even: libc::c_double = 0.;
        let mut k0_odd: libc::c_double = 0.;
        let mut k1_even: libc::c_double = 0.;
        let mut k1_odd: libc::c_double = 0.;
        let mut k2_even: libc::c_double = 0.;
        let mut k2_odd: libc::c_double = 0.;
        integrate_spiro(ks, xy.as_mut_ptr(), 4 as libc::c_int);
        ch = hypot(xy[0 as libc::c_int as usize], xy[1 as libc::c_int as usize]);
        th = atan2(xy[1 as libc::c_int as usize], xy[0 as libc::c_int as usize]);
        l = ch / seg_ch;
        th_even = 0.5f64 * *ks.offset(0 as libc::c_int as isize)
            + 1.0f64 / 48 as libc::c_int as libc::c_double * *ks.offset(2 as libc::c_int as isize);
        th_odd = 0.125f64 * *ks.offset(1 as libc::c_int as isize)
            + 1.0f64 / 384 as libc::c_int as libc::c_double * *ks.offset(3 as libc::c_int as isize)
            - th;
        (*ends.offset(0 as libc::c_int as isize))[0 as libc::c_int as usize] = th_even - th_odd;
        (*ends.offset(1 as libc::c_int as isize))[0 as libc::c_int as usize] = th_even + th_odd;
        k0_even = l
            * (*ks.offset(0 as libc::c_int as isize)
                + 0.125f64 * *ks.offset(2 as libc::c_int as isize));
        k0_odd = l
            * (0.5f64 * *ks.offset(1 as libc::c_int as isize)
                + 1.0f64 / 48 as libc::c_int as libc::c_double
                    * *ks.offset(3 as libc::c_int as isize));
        (*ends.offset(0 as libc::c_int as isize))[1 as libc::c_int as usize] = k0_even - k0_odd;
        (*ends.offset(1 as libc::c_int as isize))[1 as libc::c_int as usize] = k0_even + k0_odd;
        l2 = l * l;
        k1_even = l2
            * (*ks.offset(1 as libc::c_int as isize)
                + 0.125f64 * *ks.offset(3 as libc::c_int as isize));
        k1_odd = l2 * 0.5f64 * *ks.offset(2 as libc::c_int as isize);
        (*ends.offset(0 as libc::c_int as isize))[2 as libc::c_int as usize] = k1_even - k1_odd;
        (*ends.offset(1 as libc::c_int as isize))[2 as libc::c_int as usize] = k1_even + k1_odd;
        l3 = l2 * l;
        k2_even = l3 * *ks.offset(2 as libc::c_int as isize);
        k2_odd = l3 * 0.5f64 * *ks.offset(3 as libc::c_int as isize);
        (*ends.offset(0 as libc::c_int as isize))[3 as libc::c_int as usize] = k2_even - k2_odd;
        (*ends.offset(1 as libc::c_int as isize))[3 as libc::c_int as usize] = k2_even + k2_odd;
        return l;
    }
}
unsafe extern "C" fn compute_pderivs(
    mut s: *const spiro_seg,
    mut ends: *mut [libc::c_double; 4],
    mut derivs: *mut [[libc::c_double; 4]; 2],
    mut jinc: libc::c_int,
) {
    unsafe {
        let mut recip_d: libc::c_double = 2e6f64;
        let mut delta: libc::c_double = 1.0f64 / recip_d;
        let mut try_ks: [libc::c_double; 4] = [0.; 4];
        let mut try_ends: [[libc::c_double; 4]; 2] = [[0.; 4]; 2];
        let mut i: libc::c_int = 0;
        let mut j: libc::c_int = 0;
        let mut k: libc::c_int = 0;
        compute_ends(((*s).ks).as_ptr(), ends, (*s).seg_ch);
        i = 0 as libc::c_int;
        while i < jinc {
            j = 0 as libc::c_int;
            while j < 4 as libc::c_int {
                try_ks[j as usize] = (*s).ks[j as usize];
                j += 1;
                j;
            }
            try_ks[i as usize] += delta;
            compute_ends(
                try_ks.as_mut_ptr() as *const libc::c_double,
                try_ends.as_mut_ptr(),
                (*s).seg_ch,
            );
            k = 0 as libc::c_int;
            while k < 2 as libc::c_int {
                j = 0 as libc::c_int;
                while j < 4 as libc::c_int {
                    (*derivs.offset(j as isize))[k as usize][i as usize] = recip_d
                        * (try_ends[k as usize][j as usize]
                            - (*ends.offset(k as isize))[j as usize]);
                    j += 1;
                    j;
                }
                k += 1;
                k;
            }
            i += 1;
            i;
        }
    }
}
unsafe extern "C" fn mod_2pi(mut th: libc::c_double) -> libc::c_double {
    unsafe {
        let mut u: libc::c_double =
            th / (2 as libc::c_int as libc::c_double * 3.14159265358979323846264338327950288f64);
        return 2 as libc::c_int as libc::c_double
            * 3.14159265358979323846264338327950288f64
            * (u - floor(u + 0.5f64));
    }
}
unsafe extern "C" fn setup_path0(
    mut src: *const spiro_cp,
    mut dm: *mut libc::c_double,
    mut n: libc::c_int,
) -> *mut spiro_seg {
    unsafe {
        let mut current_block: u64;
        let mut i: libc::c_int = 0;
        let mut ilast: libc::c_int = 0;
        let mut n_seg: libc::c_int = 0;
        let mut z: libc::c_int = 0;
        let mut dx: libc::c_double = 0.;
        let mut dy: libc::c_double = 0.;
        let mut xmin: libc::c_double = 0.;
        let mut xmax: libc::c_double = 0.;
        let mut ymin: libc::c_double = 0.;
        let mut ymax: libc::c_double = 0.;
        let mut r: *mut spiro_seg = 0 as *mut spiro_seg;
        z = -(1 as libc::c_int);
        if (*src.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == 'z' as i32 {
            n -= 1;
            z = n;
        }
        if (*src.offset(0 as libc::c_int as isize)).ty as libc::c_int == ']' as i32
            || (*src.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == '[' as i32
        {
            return 0 as *mut spiro_seg;
        }
        if (*src.offset(0 as libc::c_int as isize)).ty as libc::c_int == 'h' as i32
            || (*src.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == 'a' as i32
        {
            return 0 as *mut spiro_seg;
        }
        n_seg = if (*src.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32 {
            n - 1 as libc::c_int
        } else {
            n
        };
        i = ((n_seg + 1 as libc::c_int) as libc::c_uint as libc::c_ulong)
            .wrapping_mul(::core::mem::size_of::<spiro_seg>() as libc::c_ulong)
            as libc::c_int;
        if i <= 0 as libc::c_int || {
            r = malloc(i as libc::c_uint as libc::c_ulong) as *mut spiro_seg;
            r.is_null()
        } {
            return 0 as *mut spiro_seg;
        }
        if *dm.offset(0 as libc::c_int as isize) < 0.9f64 {
            xmax = (*src.offset(0 as libc::c_int as isize)).x;
            xmin = xmax;
            ymax = (*src.offset(0 as libc::c_int as isize)).y;
            ymin = ymax;
            i = 0 as libc::c_int;
            while i < n_seg {
                if (*src.offset(i as isize)).ty as libc::c_int != 'z' as i32
                    && (*src.offset(i as isize)).ty as libc::c_int != 'h' as i32
                {
                    if (*src.offset(i as isize)).x < xmin {
                        xmin = (*src.offset(i as isize)).x;
                    } else if (*src.offset(i as isize)).x > xmax {
                        xmax = (*src.offset(i as isize)).x;
                    }
                    if (*src.offset(i as isize)).y < ymin {
                        ymin = (*src.offset(i as isize)).y;
                    } else if (*src.offset(i as isize)).y > ymax {
                        ymax = (*src.offset(i as isize)).y;
                    }
                }
                i += 1;
                i;
            }
            *dm.offset(1 as libc::c_int as isize) =
                (xmin + xmax) / 2 as libc::c_int as libc::c_double;
            xmax -= xmin;
            *dm.offset(2 as libc::c_int as isize) =
                (ymin + ymax) / 2 as libc::c_int as libc::c_double;
            ymax -= ymin;
            *dm.offset(0 as libc::c_int as isize) =
                fabs(if fabs(xmax) >= fabs(ymax) { xmax } else { ymax });
            *dm.offset(0 as libc::c_int as isize) /= 500.0f64;
        }
        i = 0 as libc::c_int;
        loop {
            if !(i < n_seg) {
                current_block = 721385680381463314;
                break;
            }
            if (*src.offset(i as isize)).ty as libc::c_int == 'a' as i32 {
                if (*src.offset((i + 1 as libc::c_int) as isize)).ty as libc::c_int == 'h' as i32
                    || i == n_seg - 1 as libc::c_int
                        && (*src.offset((i + 1 as libc::c_int) as isize)).ty as libc::c_int
                            == '}' as i32
                {
                    if (*src.offset(i as isize)).x
                        == (*src.offset((i + 1 as libc::c_int) as isize)).x
                        && (*src.offset(i as isize)).y
                            == (*src.offset((i + 1 as libc::c_int) as isize)).y
                    {
                        current_block = 8622554240757135458;
                        break;
                    }
                } else {
                    current_block = 8622554240757135458;
                    break;
                }
            } else if (*src.offset(i as isize)).ty as libc::c_int == 'h' as i32 {
                if !((*src.offset((i - 1 as libc::c_int) as isize)).ty as libc::c_int == 'a' as i32
                    || (i == 1 as libc::c_int
                        || (*src.offset(0 as libc::c_int as isize)).ty as libc::c_int
                            == '{' as i32))
                {
                    current_block = 8622554240757135458;
                    break;
                }
                if (*src.offset((i - 1 as libc::c_int) as isize)).x == (*src.offset(i as isize)).x
                    && (*src.offset((i - 1 as libc::c_int) as isize)).y
                        == (*src.offset(i as isize)).y
                {
                    current_block = 8622554240757135458;
                    break;
                }
            }
            (*r.offset(i as isize)).ty = (*src.offset(i as isize)).ty;
            (*r.offset(i as isize)).x = ((*src.offset(i as isize)).x
                - *dm.offset(1 as libc::c_int as isize))
                / *dm.offset(0 as libc::c_int as isize);
            (*r.offset(i as isize)).y = ((*src.offset(i as isize)).y
                - *dm.offset(2 as libc::c_int as isize))
                / *dm.offset(0 as libc::c_int as isize);
            (*r.offset(i as isize)).ks[0 as libc::c_int as usize] = 0.0f64;
            (*r.offset(i as isize)).ks[1 as libc::c_int as usize] = 0.0f64;
            (*r.offset(i as isize)).ks[2 as libc::c_int as usize] = 0.0f64;
            (*r.offset(i as isize)).ks[3 as libc::c_int as usize] = 0.0f64;
            i += 1;
            i;
        }
        match current_block {
            8622554240757135458 => {
                free(r as *mut libc::c_void);
                return 0 as *mut spiro_seg;
            }
            _ => {
                (*r.offset(n_seg as isize)).x = ((*src.offset((n_seg % n) as isize)).x
                    - *dm.offset(1 as libc::c_int as isize))
                    / *dm.offset(0 as libc::c_int as isize);
                (*r.offset(n_seg as isize)).y = ((*src.offset((n_seg % n) as isize)).y
                    - *dm.offset(2 as libc::c_int as isize))
                    / *dm.offset(0 as libc::c_int as isize);
                (*r.offset(n_seg as isize)).ty = (*src.offset((n_seg % n) as isize)).ty;
                i = 0 as libc::c_int;
                while i < n_seg {
                    if (*r.offset(i as isize)).ty as libc::c_int == 'h' as i32
                        || i == n_seg - 1 as libc::c_int
                            && i > 0 as libc::c_int
                            && (*r.offset(i as isize)).ty as libc::c_int == '}' as i32
                            && (*r.offset((i - 1 as libc::c_int) as isize)).ty as libc::c_int
                                == 'a' as i32
                    {
                        (*r.offset(i as isize)).x = (*r.offset((i - 1 as libc::c_int) as isize)).x;
                        (*r.offset(i as isize)).y = (*r.offset((i - 1 as libc::c_int) as isize)).y;
                    }
                    dx = (*r.offset((i + 1 as libc::c_int) as isize)).x - (*r.offset(i as isize)).x;
                    dy = (*r.offset((i + 1 as libc::c_int) as isize)).y - (*r.offset(i as isize)).y;
                    (*r.offset(i as isize)).seg_ch = hypot(dx, dy);
                    (*r.offset(i as isize)).seg_th = atan2(dy, dx);
                    i += 1;
                    i;
                }
                ilast = n_seg - 1 as libc::c_int;
                i = 0 as libc::c_int;
                while i < n_seg {
                    if (*r.offset(i as isize)).ty as libc::c_int == '{' as i32
                        || (*r.offset(i as isize)).ty as libc::c_int == '}' as i32
                        || (*r.offset(i as isize)).ty as libc::c_int == 'v' as i32
                    {
                        (*r.offset(i as isize)).bend_th = 0.0f64;
                    } else {
                        (*r.offset(i as isize)).bend_th = mod_2pi(
                            (*r.offset(i as isize)).seg_th - (*r.offset(ilast as isize)).seg_th,
                        );
                    }
                    ilast = i;
                    i += 1;
                    i;
                }
                if z >= 0 as libc::c_int {
                    (*r.offset(z as isize)).ty = 'z' as i32 as libc::c_char;
                }
                return r;
            }
        };
    }
}
unsafe extern "C" fn bandec11(mut m: *mut bandmat, mut perm: *mut libc::c_int, mut n: libc::c_int) {
    unsafe {
        let mut i: libc::c_int = 0;
        let mut j: libc::c_int = 0;
        let mut k: libc::c_int = 0;
        let mut l: libc::c_int = 0;
        let mut pivot: libc::c_int = 0;
        let mut pivot_val: libc::c_double = 0.;
        let mut pivot_scale: libc::c_double = 0.;
        let mut tmp: libc::c_double = 0.;
        let mut x: libc::c_double = 0.;
        i = 0 as libc::c_int;
        while i < 5 as libc::c_int {
            j = 0 as libc::c_int;
            while j < i + 6 as libc::c_int {
                (*m.offset(i as isize)).a[j as usize] =
                    (*m.offset(i as isize)).a[(j + 5 as libc::c_int - i) as usize];
                j += 1;
                j;
            }
            while j < 11 as libc::c_int {
                (*m.offset(i as isize)).a[j as usize] = 0.0f64;
                j += 1;
                j;
            }
            i += 1;
            i;
        }
        l = 5 as libc::c_int;
        k = 0 as libc::c_int;
        while k < n {
            pivot = k;
            pivot_val = (*m.offset(k as isize)).a[0 as libc::c_int as usize];
            l = if l < n { l + 1 as libc::c_int } else { n };
            j = k + 1 as libc::c_int;
            while j < l {
                if fabs((*m.offset(j as isize)).a[0 as libc::c_int as usize]) > fabs(pivot_val) {
                    pivot_val = (*m.offset(j as isize)).a[0 as libc::c_int as usize];
                    pivot = j;
                }
                j += 1;
                j;
            }
            *perm.offset(k as isize) = pivot;
            if pivot != k {
                j = 0 as libc::c_int;
                while j < 11 as libc::c_int {
                    tmp = (*m.offset(k as isize)).a[j as usize];
                    (*m.offset(k as isize)).a[j as usize] =
                        (*m.offset(pivot as isize)).a[j as usize];
                    (*m.offset(pivot as isize)).a[j as usize] = tmp;
                    j += 1;
                    j;
                }
            }
            if fabs(pivot_val) < 1e-12f64 {
                pivot_val = 1e-12f64;
            }
            pivot_scale = 1.0f64 / pivot_val;
            i = k + 1 as libc::c_int;
            while i < l {
                x = (*m.offset(i as isize)).a[0 as libc::c_int as usize] * pivot_scale;
                (*m.offset(k as isize)).al[(i - k - 1 as libc::c_int) as usize] = x;
                j = 1 as libc::c_int;
                while j < 11 as libc::c_int {
                    (*m.offset(i as isize)).a[(j - 1 as libc::c_int) as usize] =
                        (*m.offset(i as isize)).a[j as usize]
                            - x * (*m.offset(k as isize)).a[j as usize];
                    j += 1;
                    j;
                }
                (*m.offset(i as isize)).a[10 as libc::c_int as usize] = 0.0f64;
                i += 1;
                i;
            }
            k += 1;
            k;
        }
    }
}
unsafe extern "C" fn banbks11(
    mut m: *const bandmat,
    mut perm: *const libc::c_int,
    mut v: *mut libc::c_double,
    mut n: libc::c_int,
) {
    unsafe {
        let mut i: libc::c_int = 0;
        let mut k: libc::c_int = 0;
        let mut l: libc::c_int = 0;
        let mut tmp: libc::c_double = 0.;
        let mut x: libc::c_double = 0.;
        l = 5 as libc::c_int;
        k = 0 as libc::c_int;
        while k < n {
            i = *perm.offset(k as isize);
            if i != k {
                tmp = *v.offset(k as isize);
                *v.offset(k as isize) = *v.offset(i as isize);
                *v.offset(i as isize) = tmp;
            }
            if l < n {
                l += 1;
                l;
            }
            i = k + 1 as libc::c_int;
            while i < l {
                *v.offset(i as isize) -= (*m.offset(k as isize)).al
                    [(i - k - 1 as libc::c_int) as usize]
                    * *v.offset(k as isize);
                i += 1;
                i;
            }
            k += 1;
            k;
        }
        l = 1 as libc::c_int;
        i = n - 1 as libc::c_int;
        while i >= 0 as libc::c_int {
            x = *v.offset(i as isize);
            k = 1 as libc::c_int;
            while k < l {
                x -= (*m.offset(i as isize)).a[k as usize] * *v.offset((k + i) as isize);
                k += 1;
                k;
            }
            *v.offset(i as isize) = x / (*m.offset(i as isize)).a[0 as libc::c_int as usize];
            if l < 11 as libc::c_int {
                l += 1;
                l;
            }
            i -= 1;
            i;
        }
    }
}
unsafe extern "C" fn compute_jinc(mut ty0: libc::c_char, mut ty1: libc::c_char) -> libc::c_int {
    if ty0 as libc::c_int == 'o' as i32
        || ty1 as libc::c_int == 'o' as i32
        || ty0 as libc::c_int == ']' as i32
        || ty1 as libc::c_int == '[' as i32
        || ty0 as libc::c_int == 'h' as i32
        || ty1 as libc::c_int == 'a' as i32
    {
        return 4 as libc::c_int;
    } else if ty0 as libc::c_int == 'c' as i32 && ty1 as libc::c_int == 'c' as i32 {
        return 2 as libc::c_int;
    } else if ty1 as libc::c_int == 'c' as i32
        && (ty0 as libc::c_int == '{' as i32
            || ty0 as libc::c_int == 'v' as i32
            || ty0 as libc::c_int == '[' as i32
            || ty0 as libc::c_int == 'a' as i32)
        || ty0 as libc::c_int == 'c' as i32
            && (ty1 as libc::c_int == '}' as i32
                || ty1 as libc::c_int == 'v' as i32
                || ty1 as libc::c_int == ']' as i32
                || ty1 as libc::c_int == 'h' as i32)
    {
        return 1 as libc::c_int;
    } else {
        return 0 as libc::c_int;
    };
}
unsafe extern "C" fn count_vec(
    mut s: *const spiro_seg,
    mut jinca: *mut libc::c_int,
    mut nseg: libc::c_int,
) -> libc::c_int {
    unsafe {
        let mut i: libc::c_int = 0;
        let mut n: libc::c_int = 0;
        n = 0 as libc::c_int;
        i = 0 as libc::c_int;
        while i < nseg {
            let ref mut fresh3 = *jinca.offset(i as isize);
            *fresh3 = compute_jinc(
                (*s.offset(i as isize)).ty,
                (*s.offset((i + 1 as libc::c_int) as isize)).ty,
            );
            n += *fresh3;
            i += 1;
            i;
        }
        return n;
    }
}
unsafe extern "C" fn add_mat_line(
    mut m: *mut bandmat,
    mut v: *mut libc::c_double,
    mut derivs: *mut libc::c_double,
    mut x: libc::c_double,
    mut y: libc::c_double,
    mut j: libc::c_int,
    mut jj: libc::c_int,
    mut jinc: libc::c_int,
    mut nmat: libc::c_int,
) {
    unsafe {
        let mut joff: libc::c_int = 0;
        let mut k: libc::c_int = 0;
        if jj >= 0 as libc::c_int {
            jj %= nmat;
            joff = (j + 5 as libc::c_int - jj + nmat) % nmat;
            if nmat < 6 as libc::c_int {
                joff = j + 5 as libc::c_int - jj;
            } else if nmat == 6 as libc::c_int {
                joff = 2 as libc::c_int + (j + 3 as libc::c_int - jj + nmat) % nmat;
            }
            *v.offset(jj as isize) += x;
            k = 0 as libc::c_int;
            while k < jinc {
                (*m.offset(jj as isize)).a[(joff + k) as usize] += y * *derivs.offset(k as isize);
                k += 1;
                k;
            }
        }
    }
}
unsafe extern "C" fn spiro_iter(
    mut s: *mut spiro_seg,
    mut m: *mut bandmat,
    mut perm: *mut libc::c_int,
    mut v: *mut libc::c_double,
    mut jinca: *mut libc::c_int,
    mut n: libc::c_int,
    mut cyclic: libc::c_int,
    mut nmat: libc::c_int,
) -> libc::c_double {
    unsafe {
        let mut l: libc::c_uint = 0;
        let mut i: libc::c_int = 0;
        let mut j: libc::c_int = 0;
        let mut jthl: libc::c_int = 0;
        let mut jthr: libc::c_int = 0;
        let mut jk0l: libc::c_int = 0;
        let mut jk0r: libc::c_int = 0;
        let mut jk1l: libc::c_int = 0;
        let mut jk1r: libc::c_int = 0;
        let mut jk2l: libc::c_int = 0;
        let mut jk2r: libc::c_int = 0;
        let mut jinc: libc::c_int = 0;
        let mut jj: libc::c_int = 0;
        let mut k: libc::c_int = 0;
        let mut n_invert: libc::c_int = 0;
        let mut ty0: libc::c_char = 0;
        let mut ty1: libc::c_char = 0;
        let mut dk: libc::c_double = 0.;
        let mut norm: libc::c_double = 0.;
        let mut th: libc::c_double = 0.;
        let mut ends: [[libc::c_double; 4]; 2] = [[0.; 4]; 2];
        let mut derivs: [[[libc::c_double; 4]; 2]; 4] = [[[0.; 4]; 2]; 4];
        i = 0 as libc::c_int;
        while i < nmat {
            *v.offset(i as isize) = 0.0f64;
            j = 0 as libc::c_int;
            while j < 11 as libc::c_int {
                (*m.offset(i as isize)).a[j as usize] = 0.0f64;
                j += 1;
                j;
            }
            j = 0 as libc::c_int;
            while j < 5 as libc::c_int {
                (*m.offset(i as isize)).al[j as usize] = 0.0f64;
                j += 1;
                j;
            }
            i += 1;
            i;
        }
        jj = 0 as libc::c_int;
        j = jj;
        i = j;
        if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == 'o' as i32 {
            jj = nmat - 2 as libc::c_int;
        } else if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == 'c' as i32 {
            jj = nmat - 1 as libc::c_int;
        } else if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == '[' as i32
            || (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == 'a' as i32
        {
            if cyclic != 0 {
                let mut current_block_21: u64;
                i = 0 as libc::c_int;
                loop {
                    if !(i < n) {
                        current_block_21 = 7172762164747879670;
                        break;
                    }
                    match (*s.offset(i as isize)).ty as libc::c_int {
                        111 => {
                            jj -= 1;
                            jj;
                            current_block_21 = 17656848751191360210;
                        }
                        99 => {
                            current_block_21 = 17656848751191360210;
                        }
                        118 => {
                            current_block_21 = 17445073145183142282;
                        }
                        _ => {
                            jj += *jinca.offset(i as isize);
                            j += *jinca.offset(i as isize);
                            i += 1;
                            i;
                            continue;
                        }
                    }
                    match current_block_21 {
                        17656848751191360210 => {
                            jj -= 1;
                            jj;
                        }
                        _ => {}
                    }
                    jj = (jj + nmat) % nmat;
                    j %= nmat;
                    current_block_21 = 14763689060501151050;
                    break;
                }
                match current_block_21 {
                    7172762164747879670 => {
                        jj = 0 as libc::c_int;
                        j = jj;
                        i = j;
                    }
                    _ => {}
                }
            }
        }
        k = 0 as libc::c_int;
        while k < n {
            i %= n;
            ty0 = (*s.offset(i as isize)).ty;
            ty1 = (*s.offset((i + 1 as libc::c_int) as isize)).ty;
            jinc = *jinca.offset(i as isize);
            th = (*s.offset(i as isize)).bend_th;
            jk2l = -(1 as libc::c_int);
            jk1l = jk2l;
            jk0l = jk1l;
            jthl = jk0l;
            jk2r = -(1 as libc::c_int);
            jk1r = jk2r;
            jk0r = jk1r;
            jthr = jk0r;
            compute_pderivs(
                &mut *s.offset(i as isize),
                ends.as_mut_ptr(),
                derivs.as_mut_ptr(),
                jinc,
            );
            if ty0 as libc::c_int == 'o' as i32
                || ty0 as libc::c_int == 'c' as i32
                || ty0 as libc::c_int == '[' as i32
                || ty0 as libc::c_int == ']' as i32
                || ty0 as libc::c_int == 'a' as i32
                || ty0 as libc::c_int == 'h' as i32
            {
                let fresh4 = jj;
                jj = jj + 1;
                jthl = fresh4;
                jthl %= nmat;
                jj %= nmat;
                let fresh5 = jj;
                jj = jj + 1;
                jk0l = fresh5;
                if ty0 as libc::c_int == 'o' as i32 {
                    jj %= nmat;
                    let fresh6 = jj;
                    jj = jj + 1;
                    jk1l = fresh6;
                    let fresh7 = jj;
                    jj = jj + 1;
                    jk2l = fresh7;
                }
            }
            if jinc == 4 as libc::c_int {
                if ty0 as libc::c_int == 'c' as i32
                    || ty0 as libc::c_int == 'v' as i32
                    || ty0 as libc::c_int == '[' as i32
                    || ty0 as libc::c_int == 'a' as i32
                    || ty0 as libc::c_int == '{' as i32
                {
                    if ty0 as libc::c_int != 'c' as i32 {
                        let fresh8 = jj;
                        jj = jj + 1;
                        jk1l = fresh8;
                    }
                    let fresh9 = jj;
                    jj = jj + 1;
                    jk2l = fresh9;
                }
                if ty1 as libc::c_int == 'c' as i32
                    || ty1 as libc::c_int == 'v' as i32
                    || ty1 as libc::c_int == ']' as i32
                    || ty1 as libc::c_int == 'h' as i32
                    || ty1 as libc::c_int == '}' as i32
                {
                    if ty1 as libc::c_int != 'c' as i32 {
                        let fresh10 = jj;
                        jj = jj + 1;
                        jk1r = fresh10;
                    }
                    let fresh11 = jj;
                    jj = jj + 1;
                    jk2r = fresh11;
                }
            }
            if ty1 as libc::c_int == 'o' as i32
                || ty1 as libc::c_int == 'c' as i32
                || ty1 as libc::c_int == '[' as i32
                || ty1 as libc::c_int == ']' as i32
                || ty1 as libc::c_int == 'a' as i32
                || ty1 as libc::c_int == 'h' as i32
            {
                jj %= nmat;
                jthr = jj;
                jk0r = (jj + 1 as libc::c_int) % nmat;
                if ty1 as libc::c_int == 'o' as i32 {
                    jk1r = (jj + 2 as libc::c_int) % nmat;
                    jk2r = (jj + 3 as libc::c_int) % nmat;
                }
            }
            add_mat_line(
                m,
                v,
                (derivs[0 as libc::c_int as usize][0 as libc::c_int as usize]).as_mut_ptr(),
                th - ends[0 as libc::c_int as usize][0 as libc::c_int as usize],
                1 as libc::c_int as libc::c_double,
                j,
                jthl,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[1 as libc::c_int as usize][0 as libc::c_int as usize]).as_mut_ptr(),
                ends[0 as libc::c_int as usize][1 as libc::c_int as usize],
                -(1 as libc::c_int) as libc::c_double,
                j,
                jk0l,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[2 as libc::c_int as usize][0 as libc::c_int as usize]).as_mut_ptr(),
                ends[0 as libc::c_int as usize][2 as libc::c_int as usize],
                -(1 as libc::c_int) as libc::c_double,
                j,
                jk1l,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[3 as libc::c_int as usize][0 as libc::c_int as usize]).as_mut_ptr(),
                ends[0 as libc::c_int as usize][3 as libc::c_int as usize],
                -(1 as libc::c_int) as libc::c_double,
                j,
                jk2l,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[0 as libc::c_int as usize][1 as libc::c_int as usize]).as_mut_ptr(),
                -ends[1 as libc::c_int as usize][0 as libc::c_int as usize],
                1 as libc::c_int as libc::c_double,
                j,
                jthr,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[1 as libc::c_int as usize][1 as libc::c_int as usize]).as_mut_ptr(),
                -ends[1 as libc::c_int as usize][1 as libc::c_int as usize],
                1 as libc::c_int as libc::c_double,
                j,
                jk0r,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[2 as libc::c_int as usize][1 as libc::c_int as usize]).as_mut_ptr(),
                -ends[1 as libc::c_int as usize][2 as libc::c_int as usize],
                1 as libc::c_int as libc::c_double,
                j,
                jk1r,
                jinc,
                nmat,
            );
            add_mat_line(
                m,
                v,
                (derivs[3 as libc::c_int as usize][1 as libc::c_int as usize]).as_mut_ptr(),
                -ends[1 as libc::c_int as usize][3 as libc::c_int as usize],
                1 as libc::c_int as libc::c_double,
                j,
                jk2r,
                jinc,
                nmat,
            );
            if jthl >= 0 as libc::c_int {
                *v.offset(jthl as isize) = mod_2pi(*v.offset(jthl as isize));
            }
            if jthr >= 0 as libc::c_int {
                *v.offset(jthr as isize) = mod_2pi(*v.offset(jthr as isize));
            }
            j += jinc;
            j %= nmat;
            i += 1;
            i;
            k += 1;
            k;
        }
        if cyclic != 0 {
            l = (::core::mem::size_of::<bandmat>() as libc::c_ulong)
                .wrapping_mul(nmat as libc::c_uint as libc::c_ulong)
                as libc::c_uint;
            memcpy(
                m.offset(nmat as isize) as *mut libc::c_void,
                m as *const libc::c_void,
                l as libc::c_ulong,
            );
            memcpy(
                m.offset((2 as libc::c_int * nmat) as isize) as *mut libc::c_void,
                m as *const libc::c_void,
                l as libc::c_ulong,
            );
            l = (::core::mem::size_of::<libc::c_double>() as libc::c_ulong)
                .wrapping_mul(nmat as libc::c_uint as libc::c_ulong)
                as libc::c_uint;
            memcpy(
                v.offset(nmat as isize) as *mut libc::c_void,
                v as *const libc::c_void,
                l as libc::c_ulong,
            );
            memcpy(
                v.offset((2 as libc::c_int * nmat) as isize) as *mut libc::c_void,
                v as *const libc::c_void,
                l as libc::c_ulong,
            );
            n_invert = 3 as libc::c_int * nmat;
            j = nmat;
        } else {
            n_invert = nmat;
            j = 0 as libc::c_int;
        }
        bandec11(m, perm, n_invert);
        banbks11(m, perm, v, n_invert);
        norm = 0.0f64;
        i = 0 as libc::c_int;
        while i < n {
            jinc = *jinca.offset(i as isize);
            k = 0 as libc::c_int;
            while k < jinc {
                let fresh12 = j;
                j = j + 1;
                dk = *v.offset(fresh12 as isize);
                (*s.offset(i as isize)).ks[k as usize] += dk;
                norm += dk * dk;
                if __inline_isfinited((*s.offset(i as isize)).ks[k as usize]) == 0 as libc::c_int {
                    return (*s.offset(i as isize)).ks[k as usize];
                }
                k += 1;
                k;
            }
            (*s.offset(i as isize)).ks[0 as libc::c_int as usize] =
                2.0f64 * mod_2pi((*s.offset(i as isize)).ks[0 as libc::c_int as usize] / 2.0f64);
            i += 1;
            i;
        }
        return norm;
    }
}
unsafe extern "C" fn solve_spiro(mut s: *mut spiro_seg, mut n: libc::c_int) -> libc::c_int {
    unsafe {
        let mut i: libc::c_int = 0;
        let mut converged: libc::c_int = 0;
        let mut cyclic: libc::c_int = 0;
        let mut nmat: libc::c_int = 0;
        let mut n_alloc: libc::c_int = 0;
        let mut nseg: libc::c_int = 0;
        let mut z: libc::c_int = 0;
        let mut m: *mut bandmat = 0 as *mut bandmat;
        let mut v: *mut libc::c_double = 0 as *mut libc::c_double;
        let mut perm: *mut libc::c_int = 0 as *mut libc::c_int;
        let mut jinca: *mut libc::c_int = 0 as *mut libc::c_int;
        let mut norm: libc::c_double = 0.;
        converged = 0 as libc::c_int;
        i = converged;
        z = -(1 as libc::c_int);
        if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32 {
            nseg = n - 1 as libc::c_int;
        } else {
            if (*s.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == 'z' as i32 {
                n -= 1;
                z = n;
                (*s.offset(z as isize)).ty = (*s.offset(0 as libc::c_int as isize)).ty;
            }
            nseg = n;
        }
        if nseg <= 1 as libc::c_int {
            converged = 1 as libc::c_int;
        } else {
            jinca = malloc(
                (::core::mem::size_of::<libc::c_int>() as libc::c_ulong)
                    .wrapping_mul(nseg as libc::c_ulong) as libc::c_uint
                    as libc::c_ulong,
            ) as *mut libc::c_int;
            if !jinca.is_null() {
                nmat = count_vec(s, jinca, nseg);
                if nmat == 0 as libc::c_int {
                    converged = 1 as libc::c_int;
                } else {
                    n_alloc = nmat;
                    cyclic = 0 as libc::c_int;
                    if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int != '{' as i32
                        && (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int != 'v' as i32
                    {
                        n_alloc *= 3 as libc::c_int;
                        cyclic += 1;
                        cyclic;
                    }
                    if n_alloc < 5 as libc::c_int {
                        n_alloc = 5 as libc::c_int;
                    }
                    m = malloc(
                        (::core::mem::size_of::<bandmat>() as libc::c_ulong)
                            .wrapping_mul(n_alloc as libc::c_uint as libc::c_ulong),
                    ) as *mut bandmat;
                    v = malloc(
                        (::core::mem::size_of::<libc::c_double>() as libc::c_ulong)
                            .wrapping_mul(n_alloc as libc::c_uint as libc::c_ulong),
                    ) as *mut libc::c_double;
                    perm = malloc(
                        (::core::mem::size_of::<libc::c_int>() as libc::c_ulong)
                            .wrapping_mul(n_alloc as libc::c_uint as libc::c_ulong),
                    ) as *mut libc::c_int;
                    if !m.is_null() && !v.is_null() && !perm.is_null() {
                        loop {
                            let fresh13 = i;
                            i = i + 1;
                            if !(fresh13 < 60 as libc::c_int) {
                                break;
                            }
                            norm = spiro_iter(s, m, perm, v, jinca, nseg, cyclic, nmat);
                            if (__inline_isfinited(norm)) == 0 as libc::c_int {
                                break;
                            }
                            if !(norm < 1e-12f64) {
                                continue;
                            }
                            converged = 1 as libc::c_int;
                            break;
                        }
                    }
                    free(perm as *mut libc::c_void);
                    free(v as *mut libc::c_void);
                    free(m as *mut libc::c_void);
                }
                free(jinca as *mut libc::c_void);
            }
        }
        if z >= 0 as libc::c_int {
            (*s.offset(z as isize)).ty = 'z' as i32 as libc::c_char;
        }
        return converged;
    }
}
unsafe extern "C" fn spiro_seg_to_bpath1(
    mut ks: *const libc::c_double,
    mut dm: *mut libc::c_double,
    mut di: *mut libc::c_double,
    mut x0: libc::c_double,
    mut y0: libc::c_double,
    mut x1: libc::c_double,
    mut y1: libc::c_double,
    mut bc: *mut bezctx,
    mut ncq: libc::c_int,
    mut si: libc::c_int,
    mut depth: libc::c_int,
) {
    unsafe {
        let mut bend: libc::c_double = 0.;
        let mut seg_ch: libc::c_double = 0.;
        let mut seg_th: libc::c_double = 0.;
        let mut ch: libc::c_double = 0.;
        let mut th: libc::c_double = 0.;
        let mut scale: libc::c_double = 0.;
        let mut rot: libc::c_double = 0.;
        let mut th_even: libc::c_double = 0.;
        let mut th_odd: libc::c_double = 0.;
        let mut ul: libc::c_double = 0.;
        let mut vl: libc::c_double = 0.;
        let mut ur: libc::c_double = 0.;
        let mut vr: libc::c_double = 0.;
        let mut thsub: libc::c_double = 0.;
        let mut xmid: libc::c_double = 0.;
        let mut ymid: libc::c_double = 0.;
        let mut cth: libc::c_double = 0.;
        let mut sth: libc::c_double = 0.;
        let mut ksub: [libc::c_double; 4] = [0.; 4];
        let mut xysub: [libc::c_double; 2] = [0.; 2];
        let mut xy: [libc::c_double; 2] = [0.; 2];
        bend = fabs(*ks.offset(0 as libc::c_int as isize))
            + fabs(0.5f64 * *ks.offset(1 as libc::c_int as isize))
            + fabs(0.125f64 * *ks.offset(2 as libc::c_int as isize))
            + fabs(
                1.0f64 / 48 as libc::c_int as libc::c_double
                    * *ks.offset(3 as libc::c_int as isize),
            );
        if bend <= 1e-8f64 {
            if *di.offset(3 as libc::c_int as isize) < x1
                && x1 < *di.offset(4 as libc::c_int as isize)
                && *di.offset(6 as libc::c_int as isize) < y1
                && y1 < *di.offset(7 as libc::c_int as isize)
            {
                bezctx_lineto(
                    bc,
                    *di.offset(2 as libc::c_int as isize),
                    *di.offset(5 as libc::c_int as isize),
                    si,
                );
            } else {
                bezctx_lineto(
                    bc,
                    x1 * *dm.offset(0 as libc::c_int as isize)
                        + *dm.offset(1 as libc::c_int as isize),
                    y1 * *dm.offset(0 as libc::c_int as isize)
                        + *dm.offset(2 as libc::c_int as isize),
                    si,
                );
            }
        } else {
            seg_ch = hypot(x1 - x0, y1 - y0);
            seg_th = atan2(y1 - y0, x1 - x0);
            integrate_spiro(ks, xy.as_mut_ptr(), 4 as libc::c_int);
            ch = hypot(xy[0 as libc::c_int as usize], xy[1 as libc::c_int as usize]);
            th = atan2(xy[1 as libc::c_int as usize], xy[0 as libc::c_int as usize]);
            scale = seg_ch / ch;
            rot = seg_th - th;
            if ncq == 0 as libc::c_int
                && (depth > 5 as libc::c_int || bend < *di.offset(0 as libc::c_int as isize))
            {
                th_even = 1.0f64 / 384 as libc::c_int as libc::c_double
                    * *ks.offset(3 as libc::c_int as isize)
                    + 1.0f64 / 8 as libc::c_int as libc::c_double
                        * *ks.offset(1 as libc::c_int as isize)
                    + rot;
                th_odd = 1.0f64 / 48 as libc::c_int as libc::c_double
                    * *ks.offset(2 as libc::c_int as isize)
                    + 0.5f64 * *ks.offset(0 as libc::c_int as isize);
                ul = scale * (1.0f64 / 3 as libc::c_int as libc::c_double) * cos(th_even - th_odd);
                vl = scale * (1.0f64 / 3 as libc::c_int as libc::c_double) * sin(th_even - th_odd);
                ur = scale * (1.0f64 / 3 as libc::c_int as libc::c_double) * cos(th_even + th_odd);
                vr = scale * (1.0f64 / 3 as libc::c_int as libc::c_double) * sin(th_even + th_odd);
                if *di.offset(3 as libc::c_int as isize) < x1
                    && x1 < *di.offset(4 as libc::c_int as isize)
                    && *di.offset(6 as libc::c_int as isize) < y1
                    && y1 < *di.offset(7 as libc::c_int as isize)
                {
                    bezctx_curveto(
                        bc,
                        (x0 + ul) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        (y0 + vl) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        (x1 - ur) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        (y1 - vr) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        *di.offset(2 as libc::c_int as isize),
                        *di.offset(5 as libc::c_int as isize),
                        si,
                    );
                } else {
                    bezctx_curveto(
                        bc,
                        (x0 + ul) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        (y0 + vl) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        (x1 - ur) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        (y1 - vr) * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        x1 * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        y1 * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        si,
                    );
                }
            } else {
                ksub[0 as libc::c_int as usize] = 0.5f64 * *ks.offset(0 as libc::c_int as isize)
                    - 0.125f64 * *ks.offset(1 as libc::c_int as isize)
                    + 1.0f64 / 64 as libc::c_int as libc::c_double
                        * *ks.offset(2 as libc::c_int as isize)
                    - 1.0f64 / 768 as libc::c_int as libc::c_double
                        * *ks.offset(3 as libc::c_int as isize);
                ksub[1 as libc::c_int as usize] = 0.25f64 * *ks.offset(1 as libc::c_int as isize)
                    - 1.0f64 / 16 as libc::c_int as libc::c_double
                        * *ks.offset(2 as libc::c_int as isize)
                    + 1.0f64 / 128 as libc::c_int as libc::c_double
                        * *ks.offset(3 as libc::c_int as isize);
                ksub[2 as libc::c_int as usize] = 0.125f64 * *ks.offset(2 as libc::c_int as isize)
                    - 1.0f64 / 32 as libc::c_int as libc::c_double
                        * *ks.offset(3 as libc::c_int as isize);
                ksub[3 as libc::c_int as usize] = 1.0f64 / 16 as libc::c_int as libc::c_double
                    * *ks.offset(3 as libc::c_int as isize);
                thsub = rot - 0.25f64 * *ks.offset(0 as libc::c_int as isize)
                    + 1.0f64 / 32 as libc::c_int as libc::c_double
                        * *ks.offset(1 as libc::c_int as isize)
                    - 1.0f64 / 384 as libc::c_int as libc::c_double
                        * *ks.offset(2 as libc::c_int as isize)
                    + 1.0f64 / 6144 as libc::c_int as libc::c_double
                        * *ks.offset(3 as libc::c_int as isize);
                cth = 0.5f64 * scale * cos(thsub);
                sth = 0.5f64 * scale * sin(thsub);
                integrate_spiro(
                    ksub.as_mut_ptr() as *const libc::c_double,
                    xysub.as_mut_ptr(),
                    4 as libc::c_int,
                );
                xmid = x0 + cth * xysub[0 as libc::c_int as usize]
                    - sth * xysub[1 as libc::c_int as usize];
                ymid = y0
                    + cth * xysub[1 as libc::c_int as usize]
                    + sth * xysub[0 as libc::c_int as usize];
                if ncq != 0 as libc::c_int
                    && (depth > 5 as libc::c_int || bend < *di.offset(0 as libc::c_int as isize))
                {
                    if ncq < 0 as libc::c_int {
                        if *di.offset(3 as libc::c_int as isize) < x1
                            && x1 < *di.offset(4 as libc::c_int as isize)
                            && *di.offset(6 as libc::c_int as isize) < y1
                            && y1 < *di.offset(7 as libc::c_int as isize)
                        {
                            bezctx_quadto(
                                bc,
                                xmid * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                ymid * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                *di.offset(2 as libc::c_int as isize),
                                *di.offset(5 as libc::c_int as isize),
                                si,
                            );
                        } else {
                            bezctx_quadto(
                                bc,
                                xmid * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                ymid * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                x1 * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                y1 * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                si,
                            );
                        }
                    } else {
                        th_even = 1.0f64 / 384 as libc::c_int as libc::c_double
                            * *ks.offset(3 as libc::c_int as isize)
                            + 1.0f64 / 8 as libc::c_int as libc::c_double
                                * *ks.offset(1 as libc::c_int as isize)
                            + rot;
                        th_odd = 1.0f64 / 48 as libc::c_int as libc::c_double
                            * *ks.offset(2 as libc::c_int as isize)
                            + 0.5f64 * *ks.offset(0 as libc::c_int as isize);
                        ul = scale
                            * (1.0f64 / 6 as libc::c_int as libc::c_double)
                            * cos(th_even - th_odd);
                        vl = scale
                            * (1.0f64 / 6 as libc::c_int as libc::c_double)
                            * sin(th_even - th_odd);
                        ur = scale
                            * (1.0f64 / 6 as libc::c_int as libc::c_double)
                            * cos(th_even + th_odd);
                        vr = scale
                            * (1.0f64 / 6 as libc::c_int as libc::c_double)
                            * sin(th_even + th_odd);
                        bezctx_quadto(
                            bc,
                            (x0 + ul) * *dm.offset(0 as libc::c_int as isize)
                                + *dm.offset(1 as libc::c_int as isize),
                            (y0 + vl) * *dm.offset(0 as libc::c_int as isize)
                                + *dm.offset(2 as libc::c_int as isize),
                            xmid * *dm.offset(0 as libc::c_int as isize)
                                + *dm.offset(1 as libc::c_int as isize),
                            ymid * *dm.offset(0 as libc::c_int as isize)
                                + *dm.offset(2 as libc::c_int as isize),
                            si,
                        );
                        if *di.offset(3 as libc::c_int as isize) < x1
                            && x1 < *di.offset(4 as libc::c_int as isize)
                            && *di.offset(6 as libc::c_int as isize) < y1
                            && y1 < *di.offset(7 as libc::c_int as isize)
                        {
                            bezctx_quadto(
                                bc,
                                (x1 - ur) * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                (y1 - vr) * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                *di.offset(2 as libc::c_int as isize),
                                *di.offset(5 as libc::c_int as isize),
                                si,
                            );
                        } else {
                            bezctx_quadto(
                                bc,
                                (x1 - ur) * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                (y1 - vr) * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                x1 * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(1 as libc::c_int as isize),
                                y1 * *dm.offset(0 as libc::c_int as isize)
                                    + *dm.offset(2 as libc::c_int as isize),
                                si,
                            );
                        }
                    }
                } else {
                    spiro_seg_to_bpath1(
                        ksub.as_mut_ptr() as *const libc::c_double,
                        dm,
                        di,
                        x0,
                        y0,
                        xmid,
                        ymid,
                        bc,
                        ncq,
                        si,
                        depth + 1 as libc::c_int,
                    );
                    ksub[0 as libc::c_int as usize] += 0.25f64
                        * *ks.offset(1 as libc::c_int as isize)
                        + 1.0f64 / 384 as libc::c_int as libc::c_double
                            * *ks.offset(3 as libc::c_int as isize);
                    ksub[1 as libc::c_int as usize] +=
                        0.125f64 * *ks.offset(2 as libc::c_int as isize);
                    ksub[2 as libc::c_int as usize] += 1.0f64 / 16 as libc::c_int as libc::c_double
                        * *ks.offset(3 as libc::c_int as isize);
                    spiro_seg_to_bpath1(
                        ksub.as_mut_ptr() as *const libc::c_double,
                        dm,
                        di,
                        xmid,
                        ymid,
                        x1,
                        y1,
                        bc,
                        ncq,
                        si,
                        depth + 1 as libc::c_int,
                    );
                }
            }
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn spiroreverse(mut src: *mut spiro_cp, mut n: libc::c_int) -> libc::c_int {
    unsafe {
        let mut current_block: u64;
        let mut c: libc::c_char = 0;
        let mut i: libc::c_int = 0;
        let mut j: libc::c_int = 0;
        let mut x: libc::c_double = 0.;
        let mut y: libc::c_double = 0.;
        let mut tmp: *mut spiro_cp = 0 as *mut spiro_cp;
        if n > 2 as libc::c_int
            && (*src.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32
            && ((*src.offset(1 as libc::c_int as isize)).ty as libc::c_int == 'h' as i32
                || (*src.offset((n - 2 as libc::c_int) as isize)).ty as libc::c_int == 'a' as i32)
        {
            return -(1 as libc::c_int);
        }
        if (*src.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == 'z' as i32 {
            n -= 1;
            n;
        }
        i = (n as libc::c_uint as libc::c_ulong)
            .wrapping_mul(::core::mem::size_of::<spiro_cp>() as libc::c_ulong)
            as libc::c_int;
        if i <= 0 as libc::c_int || {
            tmp = malloc(i as libc::c_uint as libc::c_ulong) as *mut spiro_cp;
            tmp.is_null()
        } {
            return -(1 as libc::c_int);
        }
        i = 0 as libc::c_int;
        n -= 1;
        j = n;
        while i <= j {
            (*tmp.offset(j as isize)).ty = (*src.offset(i as isize)).ty;
            (*tmp.offset(j as isize)).x = (*src.offset(i as isize)).x;
            (*tmp.offset(j as isize)).y = (*src.offset(i as isize)).y;
            if i == j {
                break;
            }
            (*tmp.offset(i as isize)).ty = (*src.offset(j as isize)).ty;
            (*tmp.offset(i as isize)).x = (*src.offset(j as isize)).x;
            (*tmp.offset(i as isize)).y = (*src.offset(j as isize)).y;
            i += 1;
            i;
            j -= 1;
            j;
        }
        i = 0 as libc::c_int;
        loop {
            if !(i <= n) {
                current_block = 14136749492126903395;
                break;
            }
            c = (*tmp.offset(i as isize)).ty;
            if c as libc::c_int == '[' as i32 {
                (*tmp.offset(i as isize)).ty = ']' as i32 as libc::c_char;
            } else if c as libc::c_int == ']' as i32 {
                (*tmp.offset(i as isize)).ty = '[' as i32 as libc::c_char;
            } else if c as libc::c_int == '{' as i32 {
                (*tmp.offset(i as isize)).ty = '}' as i32 as libc::c_char;
            } else if c as libc::c_int == '}' as i32 {
                (*tmp.offset(i as isize)).ty = '{' as i32 as libc::c_char;
            } else if c as libc::c_int == 'h' as i32 {
                (*tmp.offset(i as isize)).ty = 'a' as i32 as libc::c_char;
                x = (*tmp.offset(i as isize)).x;
                (*tmp.offset(i as isize)).x = (*tmp.offset((i + 1 as libc::c_int) as isize)).x;
                x -= (*tmp.offset(i as isize)).x;
                y = (*tmp.offset(i as isize)).y;
                (*tmp.offset(i as isize)).y = (*tmp.offset((i + 1 as libc::c_int) as isize)).y;
                y -= (*tmp.offset(i as isize)).y;
                i += 1;
                if (*tmp.offset(i as isize)).ty as libc::c_int != 'a' as i32 {
                    current_block = 3640159645872305465;
                    break;
                }
                (*tmp.offset(i as isize)).ty = 'h' as i32 as libc::c_char;
                (*tmp.offset(i as isize)).x -= x;
                (*tmp.offset(i as isize)).y -= y;
            } else if c as libc::c_int == 'a' as i32 {
                current_block = 3640159645872305465;
                break;
            }
            i += 1;
            i;
        }
        match current_block {
            3640159645872305465 => {
                free(tmp as *mut libc::c_void);
                return -(1 as libc::c_int);
            }
            _ => {
                i = 0 as libc::c_int;
                while i <= n {
                    (*src.offset(i as isize)).ty = (*tmp.offset(i as isize)).ty;
                    (*src.offset(i as isize)).x = (*tmp.offset(i as isize)).x;
                    (*src.offset(i as isize)).y = (*tmp.offset(i as isize)).y;
                    i += 1;
                    i;
                }
                free(tmp as *mut libc::c_void);
                return 0 as libc::c_int;
            }
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn run_spiro0(
    mut src: *const spiro_cp,
    mut dm: *mut libc::c_double,
    mut ncq: libc::c_int,
    mut n: libc::c_int,
) -> *mut spiro_seg {
    unsafe {
        let mut s: *mut spiro_seg = 0 as *mut spiro_seg;
        if src.is_null() || n <= 0 as libc::c_int || ncq < 0 as libc::c_int {
            return 0 as *mut spiro_seg;
        }
        if ncq & 0x400 as libc::c_int != 0 {
            set_dm_to_1(dm);
        } else {
            *dm.offset(0 as libc::c_int as isize) = -1.0f64;
        }
        s = setup_path0(src, dm, n);
        if !s.is_null() {
            if solve_spiro(s, n) != 0 {
                return s;
            }
            free(s as *mut libc::c_void);
        }
        return 0 as *mut spiro_seg;
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn run_spiro(mut src: *const spiro_cp, mut n: libc::c_int) -> *mut spiro_seg {
    unsafe {
        let mut dm: [libc::c_double; 6] = [0.; 6];
        return run_spiro0(src, dm.as_mut_ptr(), 0x400 as libc::c_int, n);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_spiro(mut s: *mut spiro_seg) {
    unsafe {
        if !s.is_null() {
            free(s as *mut libc::c_void);
        }
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn spiro_to_bpath0(
    mut src: *const spiro_cp,
    mut s: *const spiro_seg,
    mut dm: *mut libc::c_double,
    mut ncq: libc::c_int,
    mut n: libc::c_int,
    mut bc: *mut bezctx,
) -> libc::c_int {
    unsafe {
        let mut i: libc::c_int = 0;
        let mut j: libc::c_int = 0;
        let mut lk: libc::c_int = 0;
        let mut nsegs: libc::c_int = 0;
        let mut si: libc::c_int = 0;
        let mut z: libc::c_int = 0;
        let mut di: [libc::c_double; 8] = [0.; 8];
        let mut x0: libc::c_double = 0.;
        let mut y0: libc::c_double = 0.;
        let mut x1: libc::c_double = 0.;
        let mut y1: libc::c_double = 0.;
        if s.is_null() || n <= 0 as libc::c_int || ncq < 0 as libc::c_int || bc.is_null() {
            return 0 as libc::c_int;
        }
        nsegs = n;
        if (*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32 {
            if n >= 2 as libc::c_int
                && (*s.offset((n - 2 as libc::c_int) as isize)).ty as libc::c_int == 'a' as i32
            {
                nsegs -= 1;
                nsegs;
            }
            nsegs -= 1;
            nsegs;
            z = -(1 as libc::c_int);
        } else {
            if (*s.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == 'z' as i32 {
                nsegs -= 1;
                nsegs;
            }
            z = nsegs - 1 as libc::c_int;
        }
        x1 = (*s.offset(0 as libc::c_int as isize)).x;
        x0 = x1;
        y1 = (*s.offset(0 as libc::c_int as isize)).y;
        y0 = y1;
        i = 1 as libc::c_int;
        while i < nsegs {
            if (*s.offset(i as isize)).ty as libc::c_int != 'z' as i32
                && (*s.offset(i as isize)).ty as libc::c_int != 'h' as i32
            {
                if (*s.offset(i as isize)).x < x0 {
                    x0 = (*s.offset(i as isize)).x;
                } else if (*s.offset(i as isize)).x > x1 {
                    x1 = (*s.offset(i as isize)).x;
                }
                if (*s.offset(i as isize)).y < y0 {
                    y0 = (*s.offset(i as isize)).y;
                } else if (*s.offset(i as isize)).y > y1 {
                    y1 = (*s.offset(i as isize)).y;
                }
            }
            i += 1;
            i;
        }
        x1 -= x0;
        y1 -= y0;
        di[1 as libc::c_int as usize] = if x1 >= y1 { x1 } else { y1 };
        di[1 as libc::c_int as usize] *= 0.0005f64;
        di[0 as libc::c_int as usize] = 1.0f64;
        lk = if ncq & 0x100 as libc::c_int != 0
            && (*s.offset((n - 1 as libc::c_int) as isize)).ty as libc::c_int == '}' as i32
        {
            1 as libc::c_int
        } else {
            0 as libc::c_int
        };
        si = if ncq & 0x200 as libc::c_int != 0 {
            1 as libc::c_int
        } else {
            0 as libc::c_int
        };
        ncq &= 0x7000 as libc::c_int;
        if !(ncq == 0 as libc::c_int) {
            if ncq == 0x1000 as libc::c_int {
                ncq = 0 as libc::c_int;
                di[0 as libc::c_int as usize] = 3.14159265358979323846264338327950288f64
                    / 2 as libc::c_int as libc::c_double
                    + 0.000001f64;
            } else if ncq == 0x2000 as libc::c_int {
                ncq = -(1 as libc::c_int);
            } else if ncq == 0x3000 as libc::c_int {
                ncq = -(1 as libc::c_int);
                di[0 as libc::c_int as usize] = 3.14159265358979323846264338327950288f64
                    / 2 as libc::c_int as libc::c_double
                    + 0.000001f64;
            } else if ncq == 0x4000 as libc::c_int {
                ncq = 0x10 as libc::c_int;
            }
        }
        j = 0 as libc::c_int;
        i = j;
        while i < nsegs {
            x0 = (*s.offset(i as isize)).x;
            y0 = (*s.offset(i as isize)).y;
            if i == 0 as libc::c_int {
                if !src.is_null() {
                    bezctx_moveto(
                        bc,
                        (*src.offset(0 as libc::c_int as isize)).x,
                        (*src.offset(0 as libc::c_int as isize)).y,
                        ((*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32)
                            as libc::c_int,
                        si,
                    );
                } else {
                    bezctx_moveto(
                        bc,
                        x0 * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(1 as libc::c_int as isize),
                        y0 * *dm.offset(0 as libc::c_int as isize)
                            + *dm.offset(2 as libc::c_int as isize),
                        ((*s.offset(0 as libc::c_int as isize)).ty as libc::c_int == '{' as i32)
                            as libc::c_int,
                        si,
                    );
                }
                if nsegs > 1 as libc::c_int
                    && (*s.offset(1 as libc::c_int as isize)).ty as libc::c_int == 'h' as i32
                {
                    i += 1;
                    i;
                }
            } else if (*s.offset(i as isize)).ty as libc::c_int == 'a' as i32 {
                i += 1;
                i;
            }
            if i == z {
                x1 = (*s.offset(0 as libc::c_int as isize)).x;
                y1 = (*s.offset(0 as libc::c_int as isize)).y;
            } else {
                x1 = (*s.offset((i + 1 as libc::c_int) as isize)).x;
                y1 = (*s.offset((i + 1 as libc::c_int) as isize)).y;
            }
            set_di_to_x1y1(di.as_mut_ptr(), dm, x1, y1);
            if !src.is_null() {
                if i == z {
                    di[2 as libc::c_int as usize] = (*src.offset(0 as libc::c_int as isize)).x;
                    di[5 as libc::c_int as usize] = (*src.offset(0 as libc::c_int as isize)).y;
                } else {
                    di[2 as libc::c_int as usize] =
                        (*src.offset((i + 1 as libc::c_int) as isize)).x;
                    di[5 as libc::c_int as usize] =
                        (*src.offset((i + 1 as libc::c_int) as isize)).y;
                }
            }
            bezctx_mark_knot(bc, j, si);
            spiro_seg_to_bpath1(
                ((*s.offset(i as isize)).ks).as_ptr(),
                dm,
                di.as_mut_ptr(),
                x0,
                y0,
                x1,
                y1,
                bc,
                ncq,
                si,
                0 as libc::c_int,
            );
            i += 1;
            i;
            j += 1;
            j;
        }
        if lk != 0 {
            bezctx_mark_knot(bc, j, si);
        }
        return 1 as libc::c_int;
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn spiro_to_bpath(
    mut s: *const spiro_seg,
    mut n: libc::c_int,
    mut bc: *mut bezctx,
) {
    unsafe {
        let mut dm: [libc::c_double; 6] = [0.; 6];
        set_dm_to_1(dm.as_mut_ptr());
        spiro_to_bpath0(
            0 as *const spiro_cp,
            s,
            dm.as_mut_ptr(),
            0x400 as libc::c_int,
            n,
            bc,
        );
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn LibSpiroVersion() -> *const libc::c_char {
    return b"1.5\0" as *const u8 as *const libc::c_char;
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_knot_th(
    mut s: *const spiro_seg,
    mut i: libc::c_int,
) -> libc::c_double {
    unsafe {
        let mut ends: [[libc::c_double; 4]; 2] = [[0.; 4]; 2];
        if i == 0 as libc::c_int {
            compute_ends(
                ((*s.offset(i as isize)).ks).as_ptr(),
                ends.as_mut_ptr(),
                (*s.offset(i as isize)).seg_ch,
            );
            return (*s.offset(i as isize)).seg_th
                - ends[0 as libc::c_int as usize][0 as libc::c_int as usize];
        } else {
            compute_ends(
                ((*s.offset((i - 1 as libc::c_int) as isize)).ks).as_ptr(),
                ends.as_mut_ptr(),
                (*s.offset((i - 1 as libc::c_int) as isize)).seg_ch,
            );
            return (*s.offset((i - 1 as libc::c_int) as isize)).seg_th
                + ends[1 as libc::c_int as usize][0 as libc::c_int as usize];
        };
    }
}
