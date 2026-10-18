# What the data require of any background (from janus-fits)

1. Shape: H(z) must decelerate at z > ~0.7 and accelerate since. BAO with r_d free maps this directly.
   Flat LCDM: chi2 = 1422.5 on CC + DESI DR2 + Pantheon+. Target: within ~10 of this, at <= 4 shape parameters.
2. With baryons only (Om_b = 0.049) the data demand ~0.95 of additional positively gravitating energy with
   effective w ~ 0 early and ~ -0.7 late (closure test). Any "new component" is a replacement for that, not an addition.
3. Negative-density components of any equation of state are driven to zero. A twin sector must gravitate
   positively in the background if it is to do anything useful.
4. r_d must come out ~148 Mpc with H0 ~ 67-68 (BAO + CC, no CMB input); H0 ~ 73 from local calibration remains a
   tension for every model and is not a target.
5. The matter era must reach z > 15 (JWST), and the universe must be > 13 Gyr old.
6. Growth: fsigma8 shape close to LCDM (Gold 2018: chi2 ~ 12 for 20 dof) with sigma8 ~ 0.8.
7. Falsifiable predictions beyond the fitted data: redshift drift sign at z ~ 3; D_L(z) for sirens at z ~ 1.

# What is worth keeping from Janus
- The structure mechanism: a second population that clumps and pushes ordinary matter into walls (reproduced in N-body,
  linear growth rate matches two-fluid theory).
- Nothing from its background.

# Dark-sector fifth force: N-body implementation (sim_gpu2, A18)

- DM-DM Yukawa kernel 2 beta^2 k^2/(k^2 + 1/lambda^2) added to the PM solver (GPU path; both species positive,
  s is a species label). Validated against linear two-fluid theory with no free parameters:
  expanding runs (ratio 5.2, beta^2 = 0 and 0.3): growth rate of both species matches theory to ~2% over a
  50-fold expansion; the DM/baryon split matches to 2% through the linear phase and falls ~10% short once
  nonlinear coupling starts (t ~ 8-10 t_J). Static runs show a ~20% rate deficit common to both species
  (thermal velocities, collisionless damping); not relevant for expanding production runs.
- The constant amplitude offset (~0.6) is the expansion switch-on transient, as in the Janus runs.
