## spiro-sys-rs

This is a raw port of [`libspiro`] from C to Rust with the help of `c2rust`,
trying to mimic the public API of [`spiro-sys`] to avoid build-time quirks especially on Windows.

At the time of writing, this build is based on `libspiro`'s [1.5.0] release.

Please refer to [`libspiro`] for the licensing and authorship of the original code.

[`libspiro`]: https://github.com/fontforge/libspiro
[1.5.0]: https://github.com/fontforge/libspiro/releases/tag/20240903
[`spiro-sys`]: https://github.com/Pangpang-Studio/spiro-sys
