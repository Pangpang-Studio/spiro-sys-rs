use spiro_sys_rs::{TaggedSpiroCPsToBezier, bezctx, spiro_cp};

use std::{
    fmt::Write as _,
    os::raw::{c_char, c_int},
};

macro_rules! spiro_cp {
    ({$x:literal, $y:literal, $ty:literal}) => {
        spiro_cp {
            x: $x as f64,
            y: $y as f64,
            ty: $ty as c_char,
        }
    };
}

trait BezCtx {
    fn moveto(&mut self, x: f64, y: f64, is_open: c_int);
    fn lineto(&mut self, x: f64, y: f64);
    fn quadto(&mut self, x1: f64, y1: f64, x2: f64, y2: f64);
    fn curveto(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64);
    fn mark_knot(&mut self, knot_idx: c_int);
}

#[repr(C)]
struct BezCtxAdapter<T> {
    pub raw: bezctx,
    pub data: T,
}

impl<T: BezCtx> BezCtxAdapter<T> {
    fn new(data: T) -> Self {
        BezCtxAdapter {
            raw: bezctx {
                moveto: Some(Self::_moveto),
                lineto: Some(Self::_lineto),
                quadto: Some(Self::_quadto),
                curveto: Some(Self::_curveto),
                mark_knot: Some(Self::_mark_knot),
            },
            data,
        }
    }

    unsafe extern "C" fn _moveto(bc: *mut bezctx, x: f64, y: f64, is_open: c_int) {
        unsafe {
            let this = &mut *(bc as *mut Self);
            this.data.moveto(x, y, is_open);
        }
    }

    unsafe extern "C" fn _lineto(bc: *mut bezctx, x: f64, y: f64) {
        unsafe {
            let this = &mut *(bc as *mut Self);
            this.data.lineto(x, y);
        }
    }

    unsafe extern "C" fn _quadto(bc: *mut bezctx, x1: f64, y1: f64, x2: f64, y2: f64) {
        unsafe {
            let this = &mut *(bc as *mut Self);
            this.data.quadto(x1, y1, x2, y2);
        }
    }

    unsafe extern "C" fn _curveto(
        bc: *mut bezctx,
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        x3: f64,
        y3: f64,
    ) {
        unsafe {
            let this = &mut *(bc as *mut Self);
            this.data.curveto(x1, y1, x2, y2, x3, y3);
        }
    }

    unsafe extern "C" fn _mark_knot(bc: *mut bezctx, knot_idx: c_int) {
        unsafe {
            let this = &mut *(bc as *mut Self);
            this.data.mark_knot(knot_idx);
        }
    }
}

#[derive(Debug, Default, Clone)]
struct TestBezCtx {
    buf: String,
}

impl TestBezCtx {
    const PRECISION: usize = 12;
}

impl BezCtx for TestBezCtx {
    fn moveto(&mut self, x: f64, y: f64, _is_open: c_int) {
        let p = Self::PRECISION;
        _ = writeln!(self.buf, "M {x:.p$} {y:.p$}");
    }

    fn lineto(&mut self, x: f64, y: f64) {
        let p = Self::PRECISION;
        _ = writeln!(self.buf, "L {x:.p$} {y:.p$}");
    }

    fn quadto(&mut self, x1: f64, y1: f64, x2: f64, y2: f64) {
        let p = Self::PRECISION;
        _ = writeln!(self.buf, "Q {x1:.p$} {y1:.p$}, {x2:.p$} {y2:.p$}");
    }

    fn curveto(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) {
        let p = Self::PRECISION;
        _ = writeln!(
            self.buf,
            "C {x1:.p$} {y1:.p$}, {x2:.p$} {y2:.p$}, {x3:.p$} {y3:.p$}"
        );
    }

    fn mark_knot(&mut self, _knot_idx: c_int) {
        _ = writeln!(self.buf);
    }
}

#[test]
fn spiro_to_beziers() {
    // Path from
    // https://github.com/fontforge/libspiro/blob/84ce4dfd24f0e3ee83589bbfb02723dff5c03414/tests/call-test.c#L278
    let mut path5 = vec![
        spiro_cp!({  0,   0, '{'}),
        spiro_cp!({100, 100, 'c'}),
        spiro_cp!({200, 200, '['}),
        spiro_cp!({300, 200, ']'}),
        spiro_cp!({400, 150, 'c'}),
        spiro_cp!({300, 100, '['}),
        spiro_cp!({200, 100, ']'}),
        spiro_cp!({150,  50, 'c'}),
        spiro_cp!({100,   0, '['}),
        spiro_cp!({  0,-100, ']'}),
        spiro_cp!({-50,-200, 'c'}),
        spiro_cp!({-80,-250, '}'}),
    ];

    let mut ctx = BezCtxAdapter::new(TestBezCtx::default());

    unsafe {
        TaggedSpiroCPsToBezier(path5.as_mut_ptr(), (&raw mut ctx).cast());
    }

    // You may verify the output at <https://svg-path-visualizer.netlify.app>.
    insta::assert_snapshot!(&ctx.data.buf);
}
