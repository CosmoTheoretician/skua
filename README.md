# skua

A test rig for cosmological expansion histories, and what it says about the current evidence for
evolving dark energy. Successor to [janus-background-fits](../janus-background-fits), which ruled
out the Janus model; this repository asks the open question that followed: *where, if anywhere, do
the late-universe data lean away from ΛCDM, and how robust is that lean?*

Two instruments, calibrated against each other:

1. **A fast compressed rig** (`test_model.py`): any background model is a five-line file defining
   $E(z)$; it is fitted to cosmic chronometers, DESI DR2 BAO, a supernova sample (Pantheon+ or
   DES-Dovekie), optionally a compressed Planck likelihood ($R$, $\ell_A$, $\omega_b$) and the
   Gold-2018 $f\sigma_8$ compilation. Minutes per model.
2. **The full pipeline** (`yaml/`): the same question through `cobaya` + CAMB with the Planck 2018
   low-ℓ, high-ℓ (lite) and lensing likelihoods, the full DESI DR2 BAO likelihood with covariances,
   and the DES-SN5YR, DES-Dovekie or Pantheon+ supernova likelihoods. The machinery DESI uses.

## Validation

Every component reproduces a published number before it is used:

| component | this pipeline | published |
|---|---|---|
| compressed rig, flat ΛCDM, CC + BAO + Pantheon+ | H0 = 68.1, Ω_m = 0.304, r_d = 148.3 (free) | DESI DR2 BAO: Ω_m ≈ 0.297, h r_d ≈ 101.5 |
| compressed CMB at Planck best fit | r_s(z*) = 144.1, D_M(z*) = 13868 Mpc, R = 1.7496 | 144.4, 13870, 1.7502 |
| compressed rig + CMB, flat ΛCDM | H0 = 68.6, Ω_m = 0.304, r_d = 146.9 | Planck + BAO: ≈ 67.7, 0.31, 147.1 |
| growth, flat ΛCDM, Gold 2018 | χ² = 13.3 / 22 pts, σ8 = 0.76 at Ω_m = 0.30 | the usual mild S8 tension |
| full pipeline, flat ΛCDM, Planck + DESI DR2 + DES-SN5YR | H0 = 68.2, Ω_m = 0.304, r_d = 147.5, σ8 = 0.807 | DESI DR2 II: ≈ 68.0, 0.30, 147.5 |
| Dovekie covariance round-trip | ΛCDM SN χ² = 1635.1 (cobaya) vs 1634.4 (own code) | — |

## Results

### Background sweep without the CMB (`results/sweep.md`)

Ten model families fitted to CC + BAO + Pantheon+ (no CMB), ranked by AIC. Everything with
a Λ-like component lands within ΔAIC ≈ 3 of ΛCDM; DGP is out (+9.7); the no-dark-energy family
(matter + curvature, which contains Janus and Milne) is out by +319. The one lean: every model
that lets $w$ move finds $w_0 \approx -0.91$ with $w_a < 0$, at ≈ 2σ in each supernova sample.

### With the compressed CMB

| model (DES-Dovekie) | Δχ² vs ΛCDM | extra params | best fit |
|---|---|---|---|
| constant $w$ | −1.3 | 1 | $w = -0.97$ |
| decaying DM → Λ | −0.2 | 1 | ξ ≈ 0 |
| quintessence ($w \ge -1$ enforced) | −2.5 | 2 | sits on the $w = -1$ boundary |
| free spline in $\rho_{DE}(z)$ | −8.5 | 3 | $\rho_{DE}$ 9% higher at $z = 0.5$, lower before |
| CPL $w_0 + w_a(1-a)$ | **−9.3** | 2 | $(-0.83, -0.60)$ |
| interacting, $Q \propto \rho_{DE}$ | = wCDM | 2 | exactly degenerate with wCDM in the background |
| interacting, $Q \propto \rho_{DM}$ | −7.0 | 2 | $w = -0.94$, ξ = −0.004: a 3% early/late matter shift, no crossing |

The signal is present in BAO + CMB alone (Δχ² = −6.2 without any supernovae), it is in the
*evolution* of $w$ (a constant $w$ buys nothing), and it requires an effective crossing of
$w = -1$ around $z \approx 0.5$: forbidding the crossing removes ~7 of the 9 units. A single
minimally coupled scalar field cannot produce it. Growth data at Gold-2018 precision do not
distinguish the readings (all models give growth χ² ≈ 13.4).

### Full likelihood (cobaya + CAMB) and the calibration

| data | compressed rig | full likelihood | published (DESI DR2 II; DES-Dovekie) |
|---|---|---|---|
| BAO + CMB + Pantheon+ | −6.4, $(-0.86, -0.48)$ | **−8.3**, $(-0.81, -0.71)$ | 2.8σ, $(-0.84, -0.62)$ |
| BAO + CMB + DES-SN5YR (2024 release) | — | **−19.9**, $(-0.76, -0.84)$ | −21, 4.2σ, $(-0.75, -0.86)$ |
| BAO + CMB + DES-Dovekie (2025 recalibration) | −9.3, $(-0.83, -0.60)$ | **−12.4**, $(-0.81, -0.66)$ | ≈ −11, 3.2σ |

Δχ² is w₀wₐCDM minus ΛCDM at the best fit, two extra parameters. Best fits from `cobaya`'s
minimizer; run-to-run spread ≈ 1–2 units.

## Conclusions

1. **The DESI DR2 evolving-dark-energy result reproduces.** With the full likelihoods this
   pipeline gives Δχ² = −19.9 for Planck + DESI DR2 + DES-SN5YR and −8.3 with Pantheon+, against
   the published −21 (4.2σ) and 2.8σ, with best-fit $(w_0, w_a)$ matching to ~0.05.
2. **The compressed rig is a reliable guide.** On matched data it tracks the full likelihood to
   within ≈ 2–3 units of Δχ² (always on the low side: the distance-prior CMB understates the signal modestly, as expected).
3. **The single largest effect in the question is the DES supernova calibration.** The 2025
   Dovekie recalibration leaves the best fit almost unchanged, $(w_0, w_a) \approx (-0.76, -0.9)$,
   but reduces the evidence from Δχ² ≈ −20 to ≈ −12 with everything else identical. The
   signal's location is stable; its strength depends on how ~1,800 supernovae were calibrated.
4. **Where the signal lives.** It is in BAO + CMB alone at ≈ 2.5σ; supernovae pin its amplitude.
   It is a preference for evolution with an effective phantom crossing near $z \approx 0.5$, which
   excludes thawing quintessence and the simplest interacting models, and leaves: a genuinely
   crossing effective equation of state (multi-field or non-minimally coupled), an interacting
   dark sector with evolving $w$, or an unresolved tension between datasets.
5. **Status.** The evidence is 2.5–3.3σ depending on catalogue and statistic, and the field is
   split between a frequentist "3σ and growing" reading and a Bayesian "no evidence once the extra
   parameters are penalised and Dovekie is used" reading. Nothing in this repository settles that;
   it reproduces both sides' numbers with one instrument. DESI DR3 and a unified supernova
   compilation are the data that will.

## What the data require of any model (`notes/requirements.md`)

Decelerating expansion at $z > 0.7$ and gentle acceleration since; $r_d \approx 147$ Mpc with
$H_0 \approx 67$–68 from BAO + CMB; ≈ 0.95 of positively gravitating energy beyond baryons that
acts like matter early and like Λ late (no negative-density sector of any equation of state
substitutes for it); a matter era reaching $z > 15$; $f\sigma_8$ close to ΛCDM with σ8 ≈ 0.8;
and, if the DESI lean is physical, an effective $w$ that crosses −1 near $z \approx 0.5$.

## Reproducing

    pip install numpy scipy matplotlib emcee cobaya camb getdist
    ln -s ../janus-fits/data data                      # Pantheon+, DES-Dovekie files from the janus-fits setup
    python run_all.py                                  # background sweep, no CMB
    python test_model.py cpl des cmb growth            # any model, any data combination
    cobaya-install yaml/lcdm.yaml -p ~/cobaya_packages # Planck, DESI, DES, Pantheon+ data (~1 GB)
    cobaya-run yaml/w0wa.yaml -p ~/cobaya_packages     # full-likelihood best fit
    cobaya-run yaml/w0wa_dovekie.yaml -p ~/cobaya_packages

`models/` holds the background families; `yaml/` the cobaya runs (`lcdm`, `w0wa`, each with
`_pp` for Pantheon+ and `_dovekie` for the recalibrated DES sample). Dovekie is run through
cobaya's `sn.desy5` class with a converted dataset (`DESDovekie/`, see `notes/dovekie.md`).

## Sources

- DESI Collaboration, Abdul Karim M. et al. (2025), DR2 Results II, Phys. Rev. D 112, 083515 (arXiv:2503.14738).
- Planck Collaboration (2020), 2018 results VI; likelihoods via cobaya (Torrado & Lewis 2021).
- Lewis A., Challinor A., Lasenby A. (2000), CAMB.
- Scolnic D. et al. (2022) ApJ 938, 113 (Pantheon+); DES Collaboration (2024) ApJL 973, L14 (DES-SN5YR);
  Popovic B. et al. (2025/26), DES-Dovekie recalibration.
- Sagredo B., Nesseris S., Sapone D. (2018) PRD 98, 083543 (Gold fσ8 compilation).
- Chen L., Huang Q.-G., Wang K. (2019) JCAP 02, 028 (Planck 2018 distance priors).
- Aubourg É. et al. (2015) PRD 92, 123516 (r_d fitting formula).
- Moresco M. et al. (2012, 2016); Moresco M. (2015) (cosmic chronometers).
- Discussion of the evidence: arXiv:2511.10631 (Bayesian reading), arXiv:2508.10514 and
  arXiv:2609.10567 (parametrisation surveys), Efstathiou (2024/25).
