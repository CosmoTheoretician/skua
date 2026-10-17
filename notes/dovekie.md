# Running DES-Dovekie through cobaya's sn.desy5 class

cobaya adds the statistical diagonal (MUERR_FINAL^2) to the matrix it reads, so the file must be
systematics-only. The Dovekie release ships STAT+SYS.npz as the packed upper triangle of the
*inverse* full matrix. Conversion (see the git history for the script): unpack and symmetrise the
inverse, invert, subtract diag(MUERR^2), write in the covsys text layout (N, then N^2 values),
and write the HD in the DESY5 CSV layout with MUERR_FINAL = MUERR. By construction cobaya then
reconstructs the full Dovekie STAT+SYS matrix; the LCDM SN chi2 agrees with an independent
implementation to 0.7 units.
