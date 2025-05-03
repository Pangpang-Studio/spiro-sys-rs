use ::libc;
unsafe extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _bezctx {
    pub moveto: Option<
        unsafe extern "C" fn(*mut bezctx, libc::c_double, libc::c_double, libc::c_int) -> (),
    >,
    pub lineto: Option<unsafe extern "C" fn(*mut bezctx, libc::c_double, libc::c_double) -> ()>,
    pub quadto: Option<
        unsafe extern "C" fn(
            *mut bezctx,
            libc::c_double,
            libc::c_double,
            libc::c_double,
            libc::c_double,
        ) -> (),
    >,
    pub curveto: Option<
        unsafe extern "C" fn(
            *mut bezctx,
            libc::c_double,
            libc::c_double,
            libc::c_double,
            libc::c_double,
            libc::c_double,
            libc::c_double,
        ) -> (),
    >,
    pub mark_knot: Option<unsafe extern "C" fn(*mut bezctx, libc::c_int) -> ()>,
}
pub type bezctx = _bezctx;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct curve_data {
    pub ty: libc::c_char,
    pub x0: libc::c_double,
    pub y0: libc::c_double,
    pub x1: libc::c_double,
    pub y1: libc::c_double,
    pub x2: libc::c_double,
    pub y2: libc::c_double,
    pub x3: libc::c_double,
    pub y3: libc::c_double,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct ls_bezctx {
    pub base: bezctx,
    pub cd: *mut curve_data,
    pub l: libc::c_int,
    pub max: libc::c_int,
    pub is_open: libc::c_int,
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn new_ls_bezctx(
    mut max: libc::c_int,
    mut ncq: libc::c_int,
) -> *mut ls_bezctx {
    unsafe {
        let mut r: *mut ls_bezctx = 0 as *mut ls_bezctx;
        if ncq & 0x200 as libc::c_int == 0 as libc::c_int || max < 1 as libc::c_int || {
            r = calloc(
                1 as libc::c_int as libc::c_ulong,
                ::core::mem::size_of::<ls_bezctx>() as libc::c_ulong,
            ) as *mut ls_bezctx;
            r.is_null()
        } {
            return 0 as *mut ls_bezctx;
        }
        (*r).cd = malloc(
            (max as libc::c_ulong)
                .wrapping_mul(::core::mem::size_of::<curve_data>() as libc::c_ulong),
        ) as *mut curve_data;
        if ((*r).cd).is_null() {
            free(r as *mut libc::c_void);
        } else {
            (*r).l = 0 as libc::c_int;
            (*r).max = max;
        }
        return r;
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_ls_bezctx(mut bd: *mut ls_bezctx) {
    unsafe {
        free((*bd).cd as *mut libc::c_void);
        free(bd as *mut libc::c_void);
    }
}
unsafe extern "C" fn prep_row_bc(mut bd: *mut ls_bezctx) -> libc::c_int {
    unsafe {
        let mut x: libc::c_double = 0.;
        let mut y: libc::c_double = 0.;
        let mut i: libc::c_int = 0;
        let mut t: libc::c_char = 0;
        i = (*bd).l;
        if i < 0 as libc::c_int {
            return 0 as libc::c_int;
        }
        if i >= (*bd).max || (*bd).max < 1 as libc::c_int {
            (*bd).l = -(1 as libc::c_int);
            return 0 as libc::c_int;
        }
        if i != 0 {
            t = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).ty;
            match t as libc::c_int {
                107 => {
                    x = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).x0;
                    y = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).y0;
                }
                108 | 109 => {
                    x = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).x1;
                    y = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).y1;
                }
                113 => {
                    x = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).x2;
                    y = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).y2;
                }
                99 => {
                    x = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).x3;
                    y = (*((*bd).cd).offset((i - 1 as libc::c_int) as isize)).y3;
                }
                _ => {
                    y = 0.0f64;
                    x = y;
                }
            }
            (*((*bd).cd).offset(i as isize)).x0 = x;
            (*((*bd).cd).offset(i as isize)).y0 = y;
        } else {
            let ref mut fresh0 = (*((*bd).cd).offset(0 as libc::c_int as isize)).y0;
            *fresh0 = 0.0f64;
            (*((*bd).cd).offset(0 as libc::c_int as isize)).x0 = *fresh0;
        }
        return 1 as libc::c_int;
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bezctx_moveto(
    mut bc: *mut bezctx,
    mut x: libc::c_double,
    mut y: libc::c_double,
    mut is_open: libc::c_int,
    mut si: libc::c_int,
) {
    unsafe {
        if si != 0 {
            let mut bd: *mut ls_bezctx = bc as *mut ls_bezctx;
            if prep_row_bc(bd) != 0 {
                let mut i: libc::c_int = (*bd).l;
                (*((*bd).cd).offset(i as isize)).x1 = x;
                (*((*bd).cd).offset(i as isize)).y1 = y;
                let ref mut fresh1 = (*((*bd).cd).offset(i as isize)).y3;
                *fresh1 = 0.0f64;
                let ref mut fresh2 = (*((*bd).cd).offset(i as isize)).x3;
                *fresh2 = *fresh1;
                let ref mut fresh3 = (*((*bd).cd).offset(i as isize)).y2;
                *fresh3 = *fresh2;
                (*((*bd).cd).offset(i as isize)).x2 = *fresh3;
                (*((*bd).cd).offset(i as isize)).ty = 'm' as i32 as libc::c_char;
                i += 1;
                i;
                (*bd).l = i;
                (*bd).is_open = is_open;
            }
        } else {
            ((*bc).moveto).expect("non-null function pointer")(bc, x, y, is_open);
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bezctx_lineto(
    mut bc: *mut bezctx,
    mut x: libc::c_double,
    mut y: libc::c_double,
    mut si: libc::c_int,
) {
    unsafe {
        if si != 0 {
            let mut bd: *mut ls_bezctx = bc as *mut ls_bezctx;
            if prep_row_bc(bd) != 0 {
                let mut i: libc::c_int = (*bd).l;
                (*((*bd).cd).offset(i as isize)).x1 = x;
                (*((*bd).cd).offset(i as isize)).y1 = y;
                let ref mut fresh4 = (*((*bd).cd).offset(i as isize)).y3;
                *fresh4 = 0.0f64;
                let ref mut fresh5 = (*((*bd).cd).offset(i as isize)).x3;
                *fresh5 = *fresh4;
                let ref mut fresh6 = (*((*bd).cd).offset(i as isize)).y2;
                *fresh6 = *fresh5;
                (*((*bd).cd).offset(i as isize)).x2 = *fresh6;
                (*((*bd).cd).offset(i as isize)).ty = 'l' as i32 as libc::c_char;
                i += 1;
                i;
                (*bd).l = i;
            }
        } else {
            ((*bc).lineto).expect("non-null function pointer")(bc, x, y);
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bezctx_quadto(
    mut bc: *mut bezctx,
    mut x1: libc::c_double,
    mut y1: libc::c_double,
    mut x2: libc::c_double,
    mut y2: libc::c_double,
    mut si: libc::c_int,
) {
    unsafe {
        if si != 0 {
            let mut bd: *mut ls_bezctx = bc as *mut ls_bezctx;
            if prep_row_bc(bd) != 0 {
                let mut i: libc::c_int = (*bd).l;
                (*((*bd).cd).offset(i as isize)).x1 = x1;
                (*((*bd).cd).offset(i as isize)).y1 = y1;
                (*((*bd).cd).offset(i as isize)).x2 = x2;
                (*((*bd).cd).offset(i as isize)).y2 = y2;
                let ref mut fresh7 = (*((*bd).cd).offset(i as isize)).y3;
                *fresh7 = 0.0f64;
                (*((*bd).cd).offset(i as isize)).x3 = *fresh7;
                (*((*bd).cd).offset(i as isize)).ty = 'q' as i32 as libc::c_char;
                i += 1;
                i;
                (*bd).l = i;
            }
        } else {
            ((*bc).quadto).expect("non-null function pointer")(bc, x1, y1, x2, y2);
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bezctx_curveto(
    mut bc: *mut bezctx,
    mut x1: libc::c_double,
    mut y1: libc::c_double,
    mut x2: libc::c_double,
    mut y2: libc::c_double,
    mut x3: libc::c_double,
    mut y3: libc::c_double,
    mut si: libc::c_int,
) {
    unsafe {
        if si != 0 {
            let mut bd: *mut ls_bezctx = bc as *mut ls_bezctx;
            if prep_row_bc(bd) != 0 {
                let mut i: libc::c_int = (*bd).l;
                (*((*bd).cd).offset(i as isize)).x1 = x1;
                (*((*bd).cd).offset(i as isize)).y1 = y1;
                (*((*bd).cd).offset(i as isize)).x2 = x2;
                (*((*bd).cd).offset(i as isize)).y2 = y2;
                (*((*bd).cd).offset(i as isize)).x3 = x3;
                (*((*bd).cd).offset(i as isize)).y3 = y3;
                (*((*bd).cd).offset(i as isize)).ty = 'c' as i32 as libc::c_char;
                i += 1;
                i;
                (*bd).l = i;
            }
        } else {
            ((*bc).curveto).expect("non-null function pointer")(bc, x1, y1, x2, y2, x3, y3);
        };
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bezctx_mark_knot(
    mut bc: *mut bezctx,
    mut knot_idx: libc::c_int,
    mut si: libc::c_int,
) {
    unsafe {
        if si != 0 {
            let mut bd: *mut ls_bezctx = bc as *mut ls_bezctx;
            if prep_row_bc(bd) != 0 {
                let mut i: libc::c_int = (*bd).l;
                let ref mut fresh8 = (*((*bd).cd).offset(i as isize)).y3;
                *fresh8 = 0.0f64;
                let ref mut fresh9 = (*((*bd).cd).offset(i as isize)).x3;
                *fresh9 = *fresh8;
                let ref mut fresh10 = (*((*bd).cd).offset(i as isize)).y2;
                *fresh10 = *fresh9;
                let ref mut fresh11 = (*((*bd).cd).offset(i as isize)).x2;
                *fresh11 = *fresh10;
                let ref mut fresh12 = (*((*bd).cd).offset(i as isize)).y1;
                *fresh12 = *fresh11;
                (*((*bd).cd).offset(i as isize)).x1 = *fresh12;
                (*((*bd).cd).offset(i as isize)).ty = 'k' as i32 as libc::c_char;
                i += 1;
                i;
                (*bd).l = i;
            }
        } else if ((*bc).mark_knot).is_some() {
            ((*bc).mark_knot).expect("non-null function pointer")(bc, knot_idx);
        }
    }
}
