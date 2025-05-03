use ::libc;

use crate::{bezctx, spiro_cp, spiro_seg};
unsafe extern "C" {
    fn free_spiro(s: *mut spiro_seg);
    fn run_spiro0(
        src: *const spiro_cp,
        dm: *mut libc::c_double,
        ncq: libc::c_int,
        n: libc::c_int,
    ) -> *mut spiro_seg;
    fn spiro_to_bpath0(
        src: *const spiro_cp,
        s: *const spiro_seg,
        dm: *mut libc::c_double,
        ncq: libc::c_int,
        n: libc::c_int,
        bc: *mut bezctx,
    ) -> libc::c_int;
    fn spiroreverse(src: *mut spiro_cp, n: libc::c_int) -> libc::c_int;
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn SpiroCPsToBezier(
    mut spiros: *mut spiro_cp,
    mut n: libc::c_int,
    mut isclosed: libc::c_int,
    mut bc: *mut bezctx,
) {
    unsafe {
        SpiroCPsToBezier2(spiros, n, 0x400 as libc::c_int, isclosed, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn TaggedSpiroCPsToBezier(mut spiros: *mut spiro_cp, mut bc: *mut bezctx) {
    unsafe {
        TaggedSpiroCPsToBezier2(spiros, 0x400 as libc::c_int, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn SpiroCPsToBezier0(
    mut spiros: *mut spiro_cp,
    mut n: libc::c_int,
    mut isclosed: libc::c_int,
    mut bc: *mut bezctx,
) -> libc::c_int {
    unsafe {
        return SpiroCPsToBezier2(spiros, n, 0x400 as libc::c_int, isclosed, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn TaggedSpiroCPsToBezier0(
    mut spiros: *mut spiro_cp,
    mut bc: *mut bezctx,
) -> libc::c_int {
    unsafe {
        return TaggedSpiroCPsToBezier2(spiros, 0x400 as libc::c_int, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn SpiroCPsToBezier1(
    mut spiros: *mut spiro_cp,
    mut n: libc::c_int,
    mut isclosed: libc::c_int,
    mut bc: *mut bezctx,
    mut done: *mut libc::c_int,
) {
    unsafe {
        *done = SpiroCPsToBezier2(spiros, n, 0x400 as libc::c_int, isclosed, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn TaggedSpiroCPsToBezier1(
    mut spiros: *mut spiro_cp,
    mut bc: *mut bezctx,
    mut done: *mut libc::c_int,
) {
    unsafe {
        *done = TaggedSpiroCPsToBezier2(spiros, 0x400 as libc::c_int, bc);
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn SpiroCPsToBezier2(
    mut spiros: *mut spiro_cp,
    mut n: libc::c_int,
    mut ncq: libc::c_int,
    mut isclosed: libc::c_int,
    mut bc: *mut bezctx,
) -> libc::c_int {
    unsafe {
        let mut dm: [libc::c_double; 6] = [0.; 6];
        let mut s: *mut spiro_seg = 0 as *mut spiro_seg;
        if n <= 0 as libc::c_int || ncq < 0 as libc::c_int {
            return 0 as libc::c_int;
        }
        if n > 1 as libc::c_int && ncq & 0x800 as libc::c_int != 0 && spiroreverse(spiros, n) != 0 {
            return 0 as libc::c_int;
        }
        if isclosed != 0 {
            s = run_spiro0(spiros, dm.as_mut_ptr(), ncq, n);
        } else {
            let mut oldty_start: libc::c_char = (*spiros.offset(0 as libc::c_int as isize)).ty;
            let mut oldty_end: libc::c_char = (*spiros.offset((n - 1 as libc::c_int) as isize)).ty;
            (*spiros.offset(0 as libc::c_int as isize)).ty = '{' as i32 as libc::c_char;
            (*spiros.offset((n - 1 as libc::c_int) as isize)).ty = '}' as i32 as libc::c_char;
            s = run_spiro0(spiros, dm.as_mut_ptr(), ncq, n);
            (*spiros.offset((n - 1 as libc::c_int) as isize)).ty = oldty_end;
            (*spiros.offset(0 as libc::c_int as isize)).ty = oldty_start;
        }
        if !s.is_null() {
            if spiro_to_bpath0(spiros, s, dm.as_mut_ptr(), ncq, n, bc) != 0 {
                free_spiro(s);
                return 1 as libc::c_int;
            }
            free_spiro(s);
        }
        return 0 as libc::c_int;
    }
}
#[unsafe(no_mangle)]
pub unsafe extern "C" fn TaggedSpiroCPsToBezier2(
    mut spiros: *mut spiro_cp,
    mut ncq: libc::c_int,
    mut bc: *mut bezctx,
) -> libc::c_int {
    unsafe {
        let mut dm: [libc::c_double; 6] = [0.; 6];
        let mut s: *mut spiro_seg = 0 as *mut spiro_seg;
        let mut n: libc::c_int = 0;
        n = 0 as libc::c_int;
        while (*spiros.offset(n as isize)).ty as libc::c_int != 'z' as i32
            && (*spiros.offset(n as isize)).ty as libc::c_int != '}' as i32
        {
            n += 1;
            n;
        }
        if (*spiros.offset(n as isize)).ty as libc::c_int == '}' as i32 {
            n += 1;
            n;
        }
        if n <= 0 as libc::c_int || ncq < 0 as libc::c_int {
            return 0 as libc::c_int;
        }
        if n > 1 as libc::c_int && ncq & 0x800 as libc::c_int != 0 && spiroreverse(spiros, n) != 0 {
            return 0 as libc::c_int;
        }
        s = run_spiro0(spiros, dm.as_mut_ptr(), ncq, n);
        if !s.is_null() {
            if spiro_to_bpath0(spiros, s, dm.as_mut_ptr(), ncq, n, bc) != 0 {
                free_spiro(s);
                return 1 as libc::c_int;
            }
            free_spiro(s);
        }
        return 0 as libc::c_int;
    }
}
