# DAMIT-convex
Asteroid light curve inversion code by Kaasalainen and Durech from DAMIT [website](https://astro.troja.mff.cuni.cz/projects/damit/pages/software_download).

## Description
The source codes of light-curve inversion routines together with brief manuals, example lightc urves, and the code for the direct problem are available. The code was developed by *Mikko Kaasalainen* in Fortran and converted to C by *Josef Durech*. There are two programs for light-curve inversion: **convexinv**, that optimizes all parameters and uses spherical harmonics functions for shape representation, and **conjgradinv**, that optimizes only shape and uses directly facet areas as parameters – it should be used at the final stage of the inversion process for 'polishing' the final shape model. You can also compute synthetic light curves with lcgenerator.

## License
This software is licensed under [CC Attribution 4.0 international License](https://creativecommons.org/licenses/by/4.0/legalcode).

## This fork (CheerfulUser/DAMIT-convex)
Changes to `convexinv` for the TESSELLATE asteroid shape pipeline. Without the new options and
without weights in the input, results are unchanged from upstream (outputs byte-identical, or
equal to ~1e-11 where summation order changed).

- **Weights**: a lightcurve header `n flag w [1]` gives the lightcurve a chi^2 weight `w`; with the
  fourth number `1`, every point line carries a ninth column, that point's weight.
- **Speed**: 11.6x faster per iteration (BLAS-batched derivatives and normal equations; linked
  against Accelerate on macOS, OpenBLAS elsewhere); `period_scan` ~10x faster.
- **Convergence**: with a stop condition below 1, the fit also ends once the damping exceeds 1e6
  (stalled at its minimum) instead of running to `MAX_N_ITER`.
- **Options**: `-e file` fitted parameters with 1-sigma uncertainties; `-c file` / `-i file` write /
  warm-start from the shape coefficients; `-y v0 free` YORP d(omega)/dt term (rad/day^2);
  `-r k n` Huber robust reweighting, `n` passes.
- **No point limits**: data arrays are allocated at run time (`MAX_N_OBS`, `POINTS_MAX` no longer
  apply); `MAX_LC` 10000.

### Build
`make` at the top level builds `convexinv/convexinv`, `convexinv/period_scan` and `minkowski`
(needs a C compiler, gfortran, and on Linux OpenBLAS).
