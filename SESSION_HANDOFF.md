# SESSION HANDOFF — RF ion-trap RNG research plan

**Updated:** 2026-10-28, Europe/Moscow  
**State (latest: M3a.3, 28 Oct 2026):** an exact symmetry audit passed and two octant P2 RF/DC solves completed on genuinely refined central meshes. DC vector-force change is still 2.51% between the solved levels (2.09% on generic off-axis points); fine full-domain versus fine octant DC forces differ up to 5.24%. The third, finer mesh is mesh-only because available RAM was about 1.28 GiB. **NO-GO remains for deterministic trajectory extension**. Older M3a short trajectories and every failed/superseded case are preserved. No new trajectory, stochastic dynamics, entropy, RNG or integrated-device demonstration was performed.  
**Resume rule:** M3a.3 is stopped after a numerical NO-GO. Do not run trajectories, stochastic/entropy simulations, RNG modeling, remote jobs, or further large P2 solves without explicit authorization. Recommended next bounded step is adaptive local P2 refinement/error estimation on the symmetric octant, subject to RAM; preserve fixed voltages and keep field uncertainty relative to stability margins explicit.

## Objective and required geometry

Determine whether a miniaturized silicon RF quadrupole can support a useful integrated RNG when confinement, noise, cooling, measurement, entropy, vacuum or gas, ion loading/recovery, fabrication, lifetime, total power and module size are evaluated together. A justified no-go is a valid result.

The required RF geometry is four separate branches in XY:
- two branches of \(x^2-y^2=+r_0^2\), electrically paired;
- two branches of \(x^2-y^2=-r_0^2\), paired and driven with opposite RF;
- extruded along z as a linear trap.
Additional DC structures provide axial confinement. Finite silicon thickness, finite length and truncation are later design variables; the ideal hyperbolas remain the reference.

## Recovery and state

When this resumed session started, SESSION_HANDOFF.md was absent from research/ion_trap_rng. The Russian plan and raw literature-result archive were present and read; this handoff was reconstructed from them and completed research.

The workspace contains older semiconductor-impedance and rectangular rail research. Its source/results and root literature audit describe a different geometry. Do not treat those outcomes as results for this hyperbolic RNG. The broader existing research files were not modified.

## Work completed

- Four-hyperbola parameterization, voltage convention, ideal quadrupole field and finite-domain boundary condition.
- Mathieu q/a convention, stability, pseudopotential, micromotion, and axial DC confinement with radial defocusing.
- Fixed-q scaling laws, energy/power scaling and analytic \(^{40}\mathrm{Ca}^+\) size estimates.
- Separate definitions for physical escape depth, safe-radius acceptance and full-RF escape.
- Air / controlled gas / sealed vacuum comparison and pressure–lifetime model.
- Candidate ions, source/loading/reload architectures and the gap between heater dimensions and total package.
- Entropy candidates, FDT and PSD conventions, deterministic/chaotic-motion caveats, image-current weighting fields and paired-electrode cancellation, optical readout, entropy conditioning and statistical tests.
- Whole-system energy/area/volume budget, CMOS TRNG comparison, proposed milestones M0–M8, Python layout, convergence/validation checklist and critical review.
- M2a finite 2D four-terminal hyperbola mesh, ideal quadrupole manufactured-boundary validation, four weighting potentials, differential and paired-electrode gradients, energy-derived capacitance per length, and lumped RC/RLC receiver scenarios.
- Curated source table with DOIs and scope limits; raw search notes retained separately.

## Main physics findings

1. Ideal single-ion RF motion is deterministic. A validated stochastic bath/readout mechanism is needed; chaos and PRNG seeds do not prove physical entropy.
2. \(V_{\rm RF}\) is the differential peak voltage between pairs; each pair is at \(\pm V_{\rm RF}/2\). With \(q=2|Q|V_{\rm RF}/(mr_0^2\Omega^2)\), \(\omega_r\simeq q\Omega/\sqrt8\) at small q.
3. For \(^{40}\mathrm{Ca}^+\), q=0.2 and 10 V differential peak, \(f_{\rm RF}\) is approximately 247 MHz at 10 µm, 2.47 GHz at 1 µm, 24.7 GHz at 100 nm, 247 GHz at 10 nm.
4. In that ideal example \(D_1=0.25\) eV and \(D_{0.5}=0.0625\) eV, about \(2.42k_BT\) at 300 K. \(D_\eta\) is a chosen-radius acceptance energy, not the finite trap escape depth.
5. A provisional \(20k_BT\) screen at \(\eta=0.5\) gives 82.7 V, but this is not a universal lifetime law.
6. Axial DC is essential and defocuses radial motion. Finite 3D geometry must be checked.
7. No nanometer silicon noise/emission threshold or pressure requirement is established here. A real noise–cooling–lifetime intersection is unknown.
8. Sealed-ion package precedents are cm³ scale; a micro-oven filament may already occupy \(90\times90\) µm². These do not prove mm³ integration.
9. M2a verified that paired right/left electrodes reduce the center weighting gradient to about 4.39e-4 of the single-electrode value; normalized differential sensing restores 0.39974/r0.
10. At the optional 1 mW / 100 kbit/s comparator scenario, energy is 10 nJ/bit, about \(8.6\times10^4\) above the reported 0.116 pJ/bit CMOS latch TRNG. This is a C-level comparison only, not a physical limit or predicted device result.

## Verified literature anchors

- Paul (1990), RF traps, DOI [10.1103/RevModPhys.62.531](https://doi.org/10.1103/RevModPhys.62.531); Leibfried et al. (2003), trapped-ion dynamics, [10.1103/RevModPhys.75.281](https://doi.org/10.1103/RevModPhys.75.281).
- Taylor et al. (2001), silicon MEMS QMS, [10.1116/1.1359172](https://doi.org/10.1116/1.1359172); Pau et al. (2006), doped-polysilicon cylindrical trap, [10.1103/PhysRevLett.96.120801](https://doi.org/10.1103/PhysRevLett.96.120801); Britton et al. (2009), degenerate Si traps, [10.1063/1.3254188](https://doi.org/10.1063/1.3254188); Mehta et al. (2014), CMOS-foundry surface trap, [10.1063/1.4892061](https://doi.org/10.1063/1.4892061).
- Turchette et al. (2000) noise heating, [10.1103/PhysRevA.61.063418](https://doi.org/10.1103/PhysRevA.61.063418); Brownnutt et al. (2015) surface noise, [10.1103/RevModPhys.87.1419](https://doi.org/10.1103/RevModPhys.87.1419); Sedlacek et al. (2018) distance scaling, [10.1103/PhysRevA.97.020302](https://doi.org/10.1103/PhysRevA.97.020302).
- Tian et al. (2014) scaling, [10.1002/jms.3343](https://doi.org/10.1002/jms.3343); Viehland & Goeringer (2004) collision kinetics, [10.1063/1.1691405](https://doi.org/10.1063/1.1691405); DeVoe (2009) nonthermal tails, [10.1103/PhysRevLett.102.063001](https://doi.org/10.1103/PhysRevLett.102.063001).
- Li et al. (online 2025), hyperbolic LIT simulation, [10.1021/jasms.5c00307](https://doi.org/10.1021/jasms.5c00307); Nguyen et al. (2025) trap scaling, [10.3390/e27060576](https://doi.org/10.3390/e27060576); Sterk et al. (2026) RF losses, [10.1116/5.0332336](https://doi.org/10.1116/5.0332336).
- Pironio et al. (2010) Bell-certified ion randomness, [10.1038/nature09008](https://doi.org/10.1038/nature09008); Um et al. (2013) contextuality RNG, [10.1038/srep01627](https://doi.org/10.1038/srep01627).
- Kumar et al. (2025) micro-oven, [10.1063/5.0271048](https://doi.org/10.1063/5.0271048); Schwindt et al. (2016) mini package, [10.1063/1.4948739](https://doi.org/10.1063/1.4948739); Thrasher & Schwindt (2025) sealed package, [10.1063/5.0243517](https://doi.org/10.1063/5.0243517).
- Ramo (1939) induced current, [10.1109/JRPROC.1939.228757](https://doi.org/10.1109/JRPROC.1939.228757); Wang & Wang (2013) ion-motion detection, [10.1088/0143-0807/34/3/787](https://doi.org/10.1088/0143-0807/34/3/787); Setzer et al. (2021) monolithic SPAD, [10.1063/5.0055999](https://doi.org/10.1063/5.0055999).
- NIST SP 800-90B, [10.6028/NIST.SP.800-90B](https://doi.org/10.6028/NIST.SP.800-90B); SP 800-22 Rev.1a, [10.6028/NIST.SP.800-22r1a](https://doi.org/10.6028/NIST.SP.800-22r1a); Wang et al. (2024) CMOS latch TRNG, [10.1109/TVLSI.2023.3328602](https://doi.org/10.1109/TVLSI.2023.3328602).

## Unresolved questions

Species-specific collision and charge-exchange rates; processed-silicon electric-field noise spectra; emission/breakdown for fabricated nanogaps; whether useful entropy, measurable readout and acceptable ion lifetime coexist; single-ion capture/reload probability; getter capacity/leak/outgassing; source contamination and total module budgets; final novelty after systematic review. Provisional benchmark values remain open to user adjustment.

## Files

- RESEARCH_PLAN_RU.md — full equations, model, gates and scientific review.
- RESEARCH_SUMMARY_RU.md — concise research summary.
- COMPUTATIONAL_ROADMAP.md — staged models and validation criteria.
- LITERATURE_AUDIT_RU.md — curated DOI/source-claim audit.
- LITERATURE_SESSION_NOTES.json — raw source discoveries, including duplicate/unverified leads.
- SESSION_HANDOFF.md — this record.

The root README and older research/LITERATURE_AUDIT.md are preserved, unmodified. The raw literature archive is preserved.

## Current validated state (27 October 2026)

User authorized M1 analytical code and a preliminary M1.5 screen; both bounded tasks are now complete. `RESEARCH_PLAN_RU.md`, `COMPUTATIONAL_ROADMAP.md`, and `RESEARCH_SUMMARY_RU.md` received targeted additions for early observability/entropy/cooling/lifetime/power screening, independent P/E/C criteria, entropy-source prioritization, and F/T/R/U classifications of miniaturization limits. Literature audit and raw notes were preserved; historical checkpoint text was not deleted.

M1 package: `pyproject.toml`, `configs/m1_baseline.yaml`, `src/iontrap_rng/{constants,geometry,ions,trap,floquet,scaling,screening}.py`, `tests/test_m1_analytic.py`, `scripts/run_m1_reference.py`, and `results/m1_reference.json`. 12 tests pass locally. Validation found and corrected a sign error in the radial electric-field components of the Laplace-consistent axial DC potential: E=-grad(Phi) points outward in x/y and confines along z. A finite-difference gradient and Laplace trace test now cover this. Floquet at q=.2: beta=0.1425513, determinant=0.99999999999985; q=0,a=0 is correctly treated as marginal free motion, not confinement; q=.91,a=0 is unstable.

At r0=1 µm, Vdiff=10 V peak, q=.2, axial fz=1 MHz, selected Ca+/Sr+/Yb+ outputs are RF ~2.472/1.667/1.185 GHz; small-q radial secular ~174.8/117.8/83.8 MHz; eta=.5 acceptance energy .0625 eV for all three. This is not the finite-trap escape depth.

M1.5: for Ca+, assumed A=100 nm and deff=1 µm give ideal image-current scale ~17.6 pA at ~175 MHz. A hypothetical 1 MΩ, 300 K, 1 kHz ENBW has 4.07 pA rms Johnson current. No weighting-field solution, pair-cancellation, RF feedthrough, or full receiver transfer function is included, so this is a sensitivity screen only. Conditional entropy rate remains unknown without a supported physical noise/collision model, quantization, correlations and side-information boundary.

For an explicitly assumed k=0.3–3×10^-15 m³/s at p=10^-8 Pa and 300 K, total collision rate is roughly 7×10^-4–7×10^-3 s^-1. Independent sparse Poisson event/no-event indicators have an ideal Hmin ceiling roughly 10^-3–10^-2 bit/s; not a device prediction. A 100 kbit/s binary event-indicator scenario requires Γ≈6.93×10^4/s, p≈0.096–0.96 Pa and ENBW ≥50 kHz; at the baseline 1 kHz bandwidth there are at most ~2k binary sample symbols/s. At 50 kHz the assumed 1 MΩ Johnson current noise rises to 28.8 pA rms versus the assumed 17.6 pA signal scale. The simple 99%/hour survival model requires loss probability <4.03×10^-11 per collision at the target rate. Applying the literature electric-field PSD envelope 10^-15–10^-6 V²m⁻²Hz⁻¹ to Ca+ at ~175 MHz yields ~8×10^-4–8×10^5 phonons/s and a no-cooling acceptance-energy scale ~10^8–10^-1 s. This nine-decade calculation extrapolates across different surfaces/distances/frequencies and is not a local Si prior. Cooling must offset net motional heating and can suppress the fluctuations being sampled while increasing power.

No primary entropy mechanism is selected: Johnson noise is the best bounded direct circuit control but may be sampled without an ion; collisions are the conditional ion-mediated candidate but trade event rate against RF heating and lifetime; surface charge fluctuations couple directly but the proposed silicon spectrum/correlation time is unbounded by available measurements. 100 kbit/s, 1 mW and package sizes are application scenarios only. P=physical, E=engineering, C=competitive; each boundary is tagged F/T/R/U.

## M2a completion update

User authorized only a bounded local readout study. `src/iontrap_rng/readout_fem.py` uses Gmsh and a custom SciPy P1 triangular Laplace FEM. The ideal reference used exact quadrupole Dirichlet data on an outer square. The physical finite-electrode case used separate truncated hyperbola branches with thickness 0.1r0, extent 1.5r0, in-plane cap radius 0.05r0, and a grounded enclosure at 6r0; it did not impose the ideal quadrupole at the outer boundary. The 2D capacitance is stored as F/m, with an explicit 10r0=10 µm length approximation only for receiver scenarios.

Five local mesh minima (h/r0=1/20, 1/40, 1/80, 1/160, 1/320) were run (finest: 49,135 nodes and 88,544 triangles), plus curve representation checks at 50/100/200 points per branch at h=r0/160. Finest selected outputs: ideal center Hessian error 0.173%; differential weighting gradient 0.399741/r0; electrode C'=7.29524e-11 F/m; differential-port C'=3.74953e-11 F/m. Changes from 1/160 to 1/320: differential gradient ~0.060%, electrode C' ~0.019%; higher multipole fit is less robust and is not treated as precise. **M2a.1 corrected the ideal field L2 metric**: prior 1.73% used centroid-only error integration, corrected exact element integration gives a ~2.46% plateau due to far-field coarsening. See correction below.

At Ca+ 174.784 MHz and assumed A=100 nm, the normalized differential current is ~7.03 pA peak. At 1 pF parasitic, 1 MΩ, 300 K, 2 nV/√Hz amplifier voltage noise, 2 fA/√Hz current noise and ENBW=1 kHz, unmatched parallel RC gives signal/rms noise ~0.10. A high-Q scenario can improve this but demands uncertain C/Q and resonator parameters. Conventional symmetric-pair pickup is geometrically suppressed; balanced differential sensing is a conditional candidate. No measured or device-specific receiver noise, RF feedthrough, physical entropy, RNG rate, or module power was established.

The package test suite now has **17 passing tests** (`python -m pytest tests -q`). M2a files added: `M2A_READOUT_REPORT_RU.md`, `configs/m2a_readout.yaml`, `src/iontrap_rng/readout_fem.py`, `tests/test_m2a_readout.py`, `scripts/run_m2a_readout.py`, `results/m2a_readout.json`, `results/figures/m2a_mesh_convergence.png`, and `results/figures/m2a_potential_weighting.png`. Updated: `RESEARCH_PLAN_RU.md`, `RESEARCH_SUMMARY_RU.md`, `COMPUTATIONAL_ROADMAP.md`, `pyproject.toml`, and this handoff. Literature audit/source notes and unrelated project files remain preserved.

### Recommended resume point

M2a is complete. For a future authorized stage, first specify or bound a concrete differential-sense topology: input capacitance and impedance, loaded Q, amplifier voltage/current noise, common-mode rejection, RF isolation/feedthrough, electrode length, and package parasitics. Only then decide whether one targeted 3D field/port calculation is justified. Do not infer 10 nm feasibility from continuum electrostatics: with the selected thickness ratio, the conductor is 1 nm thick at (r_0=10) nm, and continuum doping, conductivity, interface, and roughness assumptions are unresolved. No further computation is authorized by the current task.

## Files added or modified in M1/M1.5

- Modified: `RESEARCH_PLAN_RU.md`, `RESEARCH_SUMMARY_RU.md`, `COMPUTATIONAL_ROADMAP.md`, this handoff, `configs/m1_baseline.yaml`, `src/iontrap_rng/trap.py`, `src/iontrap_rng/screening.py`.
- Added: `tests/test_m1_analytic.py`, `tests/conftest.py`, `scripts/run_m1_reference.py`, `scripts/run_m15_screen.py`, `results/m1_reference.json`, `results/m15_screen.json`, `M1_ANALYTICAL_REPORT_RU.md`, and `results/figures/m1_analytical_screen.png`.
- Preserved: `LITERATURE_AUDIT_RU.md`, `LITERATURE_SESSION_NOTES.json`, the root literature audit, and unrelated impedance project content.

## Resume recommendation

This last recommendation from the pre-M2 checkpoint is superseded by the M2a completion update above. M2a is complete; the next step is to define receiver parameters before deciding whether one bounded 3D differential-port model is useful. No 3D FEM, broad sweep, stochastic trajectory campaign, or remote calculation has been run.

## M2a.1 completion update — 27 October 2026

The user authorized a bounded numerical-accuracy/readout architecture review, explicitly withholding full 3D FEM and stochastic work. Primary target remains a vacuum empty cavity inside a semiconductor/MEMS stack; the ion does not move in solid dielectric.

### Numerical audit

- Added `scripts/run_m2a1_bounded_audit.py`; outputs are in `results/m2a1_numeric_audit.json`.
- The initial M2a centroid-only field-L2 metric omitted within-triangle variation of the exact linear quadrupole gradient. M2a.1 corrected the element integral and reran only the ideal references: local hmin r0/40→r0/160 now gives 2.4644%→2.4559%; region rho<r0 gives 8.75%→8.16%. The global error plateau remains because far-field max element size is fixed at 0.4r0.
- Corrected ideal uniform reference at h/r0=0.4/0.2/0.1 gives global E-L2 3.319/1.664/0.833%, respectively (finest 16,962 nodes/33,442 triangles); observed scaling is roughly first order. No sub-0.1% claim. P2 has not been run; exact quadratic representation is a theoretical future benchmark only.
- At a finite-device hmin/r0=1/80, grounded outer boundary 4/6/10r0 barely changed center differential weighting gradient (~0.39908/r0; span below 0.003%) but changed per-electrode C′ by ~7.8%. This is not a converged package boundary study.
- Existing P1 finite geometry h refinement remains limited: between /160 and /320 differential center gradient changes ~0.060%, C′ ~0.019%, while higher multipole fit changes ~5.4% relative. Curve sampling 50/100/200 changed selected integrated observables <0.03%; no curved P2 boundaries or a posteriori residual estimator.

### Size and receiver screen

`scripts/run_m2a1_scale_screen.py` produced `results/m2a1_scale_screen.json`. For analytical Ca+, fixed q=0.2, Vdiff=10V peak, A/r0=0.1, and Lz=10r0: r0 10µm/3µm/1µm/300nm/100nm gives fsec 17.48/58.26/174.78/582.61/1747.84MHz; fRF 0.247/0.824/2.472/8.24/24.72GHz; Ramo peak current 0.703/2.344/7.033/23.445/70.335pA; approximate 2D differential trap C 3.75/1.125/0.375/0.1125/0.0375fF. These are consistent conditional scaling results, not physical device predictions. Ideal acceptance energy stays 2.5meV at eta=.1 and 62.5meV at eta=.5; the latter is only 2.42 kBT at 300K.

Published source anchors are deliberately not conflated into one receiver range: 7fF/193µm² 65nm charge buffer (Joshi et al., DOI 10.1109/ISCAS.2016.7539210); 1.02pF/low-kHz CMOS capacitance-sensor front end (Liu et al., DOI 10.1063/5.0186178); released CMOS-MEMS inductor 1.8nH at 1GHz, peak Q15 at 15GHz, SRF25.5GHz (Yang et al., DOI 10.3390/s111009798); cryogenic Penning amplifier Q up to3500, ~1mW near39MHz (Chaddha et al., DOI 10.1007/s10909-025-03309-z); superconducting Q~98004 near1MHz with ~0.8nV/√Hz amp noise (Wang et al., DOI 10.1007/s41365-025-01797-y).

The most relevant architecture evidence is the peer-reviewed 2025 room-temperature electron Paul-trap image-current experiment (Taniguchi & Noguchi, DOI 10.1103/6mpr-5s2h): separate cavity readout mode at619MHz and trap mode at3.105GHz; estimated loaded readout Q~1220; mode-selective ports suppressed trap-microwave leakage, but the current copper cavity required ~0.8–1.5W drive. This is evidence for modal separation as an architecture, not for Ca+, silicon hyperbolic electrodes, or a competitive module. Recommended future candidate is mode-selective secular readout plus balanced four-terminal differential pickup; no receiver SNR or entropy rate is claimed.

For fixed normalized motion, unresonated capacitive-node voltage is about 0.60–0.91µV with the 7fF CMOS sensor anchor, versus ~6.26–6.28nV with a 1.02pF input anchor, nearly independent of scale because current rises as frequency while capacitive impedance falls. Passive tuning L for a 7fF input decreases from ~7.7mH (10µm) to ~1.18µH (100nm); with 1.02pF it decreases from ~82.6µH to ~8.13nH. Small input capacitance improves direct voltage observability but makes a passive tuned circuit demand very large L. High Q narrows bandwidth and lengthens ringdown; it does not establish innovation or entropy rate.

### Recommended resume point

The following M2a.1 resume guidance is historical and was superseded by M2b authorization. The current resume point is the validation repair list at the end of this file. Do not use the present 3D fields for ion dynamics until local Hessian/multipole convergence and a physically grounded stack have been established.

M2a.1 files added: `M2A1_READOUT_REVIEW_RU.md`, `scripts/run_m2a1_bounded_audit.py`, `scripts/run_m2a1_scale_screen.py`, `results/m2a1_numeric_audit.json`, `results/m2a1_scale_screen.json`. Updated: `COMPUTATIONAL_ROADMAP.md`, `M2A_READOUT_REPORT_RU.md`, `src/iontrap_rng/readout_fem.py`, `tests/test_m2a_readout.py`, `results/m2a_readout.json`, and this handoff. No unrelated semiconductor impedance files were modified.

## M2b completion update — 27 October 2026

### Objective and geometry retained

The research objective is still to establish full-device physical, engineering and competitiveness limits for an integrated ion-trap RNG. The required RF geometry remains four distinct hyperbolic branches (x^2-y^2=+r_0^2) (right/left) and (x^2-y^2=-r_0^2) (top/bottom), independently addressed and driven as opposite pairs, extruded along z. This stage used (r_0=1\,\mu\mathrm m), Ca+, (L=10r_0), two axial annular DC rings and a grounded enclosure.

No published stack matching four volumetric hyperbolic silicon branches + supports/dielectrics + differential terminal routing was identified in the focused source check. Consequently `configs/m2b_preflight.yaml` declares a **vacuum-only, perfect-conductor** FEM preflight with analyst-selected (0.1r_0) backing, branch extent (1.5r_0), annular endcaps at (z=\pm5.5r_0), radii (0.5ldots1.5r_0), and an outer grounded box. These geometric ratios are hypotheses for a benchmark, not foundry rules. Dielectric supports, CMOS stack, finite conductivity, surfaces and RF contacts are omitted.

### Work and findings

- P1/P2 manufactured solution Φ=(x²−y²)/2 on ([-6,6]^2r_0), three uniform mesh sizes. P2 global E-L2 errors (2.55\cdot10^{-15}) through (2.90\cdot10^{-14}), Hessian relative errors (1.92\cdot10^{-14}) through (3.30\cdot10^{-13}). P1 global field error 8.333%, 4.167%, 2.083% and local ρ<r0 56.87%, 28.55%, 14.34%. Exact P1 element-moment integration independently agrees with FE quadrature within (4\cdot10^{-17}). This validates the manufactured solver path, not finite device fields.
- Three bounded 3D P1 tetra mesh levels calculated RF/DC solutions, six terminal weighting solutions and a grounded-reference Maxwell C matrix. On the two finest levels, differential center weighting-gradient magnitude changes ~0.21%, differential C changes 2.42%, matrix Frobenius norm changes 3.4%. RF Hessian component changes are 9–82%; axial DC Hessian changes are 24–60%; RF Hessian trace residual is ~9.9% on the fine level. The coarsest local fit is unresolved. Mesh deltas are not rigorous error bars.
- The completed 3D levels refine at the center but not along every electrode/end surface. A surface-distance refinement attempt exceeded the bounded local runtime and was stopped; no output from that attempt is used. `configs/m2b_preflight.yaml` records this limitation. Therefore the requested surface and end-fringing convergence remains incomplete.
- Fine-mesh differential weighting-gradient coefficient is 0.4135/r0 for the normalized unit-terminal weighting potential; differential capacitance is 5.36e−16 F relative to the declared grounded enclosure. For assumed A=0.1r0, Ca+ motion at 174.8 MHz gives a nominal Ramo current scale near 7.3 pA. This is a prescribed-trajectory estimate; receiver SNR, signal transfer, entropy and throughput are unknown.
- The figures are `results/figures/m2b_3d_geometry.png`, `m2b_3d_potential_slices.png`, and `m2b_p1_p2_convergence.png`. Reproducible data/config/scripts are enumerated below.
- `pytest -q`: 23 passed without warnings after removing a deprecated field-value accessor. No remote computation or unrelated impedance-project edit was made.

### Remaining tasks and exact resume instructions

1. Read `M2B_3D_FEASIBILITY_REPORT_RU.md` and current `COMPUTATIONAL_ROADMAP.md` before new work; historic M2a stop messages in older sections are superseded.
2. Repair local polynomial-fit diagnostics: report sample count, matrix rank/condition number and set a clear invalid flag for under-resolved central Hessian.
3. For the one selected geometry, refine locally at ion, hyperbolic boundaries and z ends; test P1/P2 or a second solver, true curved/isoparametric surfaces, local Hessian/multipoles and symmetry. Cross-check C from energy against flux/charge and reciprocity.
4. Seek a source-supported process/device stack. Unless one is found, retain only preflight status; do not describe this as an integrated CMOS trap.
5. After field convergence, request/confirm authorization before passing fields to dynamics. Entropy work also needs physical noise/collision model, cooling and receiver transfer/noise. Statistical tests alone do not establish entropy or security.

Files added: `M2B_3D_FEASIBILITY_REPORT_RU.md`, `configs/m2b_preflight.yaml`, `scripts/run_m2b_manufactured.py`, `scripts/run_m2b_3d_preflight.py`, `scripts/plot_m2b_3d_preflight.py`, `tests/test_m2b_manufactured.py`, `tests/test_m2b_preflight_results.py`, `results/m2b_manufactured_validation.json`, `results/m2b_3d_preflight.json`, `results/m2b_3d_fields.npz`, and three M2b figures. Updated: `RESEARCH_PLAN_RU.md`, `RESEARCH_SUMMARY_RU.md`, `COMPUTATIONAL_ROADMAP.md`, `configs/m2b_validation.yaml`, and this handoff. No unrelated semiconductor impedance research was modified.

## M2b.1 completion update — 27 October 2026

The user authorized only root-cause diagnosis and bounded local field validation, explicitly prohibiting ion trajectories, stochastic dynamics, entropy calculations, extensive sweeps and remote runs. The original four separate hyperbolic RF branches and empty vacuum ion volume remain unchanged. The only finite-geometry runs here use the pre-existing analyst-selected ideal-conductor geometry; these do not represent a qualified silicon/dielectric/device stack.

### Findings

- The legacy P1 RF trace residual near 9.93% is dominated by a poor local postprocessing estimate: quadratic fit on only 23 nodes inside `0.5r0`, normalized condition ~12.96, potential-fit residual 18.46%; changing weights and radii changes Hessian substantially. Direct P1 tetra Hessians are zero per element. A trace-free harmonic fit returning trace zero is constrained by construction and is not validation.
- The independent manufactured 3D harmonic quadratic test passes for P2 at floating-point precision for potential, field and Hessian at structured 2/4/8 cells per axis. P1 field errors converge first order but its local Hessian cannot be read directly from P1 elements.
- Two bounded finite-device P2 levels used 95,068 and 116,097 DOFs. Central `H_RF` tensor changes 0.0173%; its recovered trace fractions are `2.03e−5` and `1.70e−5`. This supports a local RF quadrupole estimate in the idealized model. It is not a complete convergence certificate because center and global sizes changed together, surfaces/end regions were not independently refined, and boundaries remain polygonal.
- The axial DC Hessian changes 97.2% in relative Frobenius norm; DC gradient recovery residuals are ~0.11–0.18. This is a hard numerical no-go for ion dynamics. The result is unresolved, not evidence of physical trap impossibility. Maxwell capacitance matrix changes 0.938%, `Cdiff` 0.90%, and differential weighting center gradient ~0.39%; all are only two-level ideal-conductor estimates.
- Old P1 Gaussian-sphere flux indicators range 0.43–9.52% on a 12×24 angular rule; no angular quadrature convergence was demonstrated. P2 Gaussian flux, strong-form elementwise Laplacian, higher multipoles, full symmetry norm and independent surface/end h-refinement remain incomplete. No multipole precision or global Laplace error bar is claimed.
- A postprocessor extension to add nodal-potential Hessian crosschecks was attempted but exceeded the local execution/automatic-review window before producing results. Existing successful artifacts are intact; the report marks that comparison unresolved.

### Validation and files

Full local suite: `..\\..\\.venv\\Scripts\\python.exe -m pytest -q` — 26 passed. Primary report: `M2B1_3D_VALIDATION_REPORT_RU.md`. New/modified code: `src/iontrap_rng/hessian.py`, `scripts/run_m2b1_hessian_diagnostics.py`, `scripts/run_m2b1_manufactured3d.py`, `scripts/run_m2b1_p2_physical_check.py`, `tests/test_m2b1_hessian.py`; the P2 geometry helper was minimally corrected to permit boundary-tolerance selection. New results: `results/m2b1_old_grid_diagnostics.json`, `results/m2b1_manufactured3d.json`, `results/m2b1_result_index.json`, `results/m2b1_p2_physical_check.json` (initial coarse run), and corrected `results/m2b1_p2_h*_hc*_v2.json`. The initial coarse/medium JSON outputs are retained but superseded because one weighting-gradient JSON key was mislabeled; use the v2 files. The separate P2 nodal-potential fit extension did not produce a result within the bounded execution window. Earlier M2b outputs remain unchanged. `COMPUTATIONAL_ROADMAP.md` and this handoff were updated. No unrelated impedance-project files were changed.

### Resume point

M2b.1 is complete and **ion dynamics is not authorized by its numerical gate**. Before dynamics, the next bounded validation should separately refine the central volume, hyperbolic electrode surfaces, annular DC surfaces and z-end regions; compare local P2 nodal-potential fits with gradient recovery; compute strong-form P2 `ΔΦ` on interior elements; refine Gaussian flux quadrature; inspect full terminal/electric-field symmetry; and converge axial DC curvature, multipoles, weighting gradient and capacitance. If convergence costs exceed local resources, record DOFs/memory/time and prepare a reproducible remote plan, without launching it. Only after <1% successive changes and independent PDE/flux checks agree should the user be asked to authorize ion-dynamics modeling. No FEM/trajectory work should resume without explicit user direction.





## M2b.2 completion update — 28 October 2026

The interrupted M2b.2 investigation resumed without deleting/replacing prior artifacts. Root cause: the generator set `Mesh.MeshSizeMin=0.18r0`, silently clamping the requested DC ring surface h values 0.15, 0.10 and 0.075r0. Corrected the global floor to the minimum active local target, disabled competing automatic Gmsh sizing heuristics, and retained independent center, RF surface, DC surface, axial end and far-field controls. Old clamped `m2b2_dc_surface_mesh_probe*` and associated P2 cases remain historical/ineffective for DC surface refinement.

Three mesh-only probes verified actual progressive refinement before P2 assembly: DC boundary triangles 2068/3756/5228; P2 DOFs 201361/208503/215157; surface median edge sizes 0.1804/0.1314/0.1125r0 for requested targets 0.15/0.10/0.075r0. Center h=0.18r0 and RF/end h=0.30r0 stayed fixed. Three effective P2 solutions were completed. At R=0.5r0 the DC diagonal Hessian per 1-V-per-ring excitation is (-9.392e-5,-9.279e-5,1.869e-4), (-9.383e-5,-9.343e-5,1.877e-4), and (-9.382e-5,-9.401e-5,1.880e-4) V/r0^2. Signs give radial defocusing and axial confinement as expected. Last-pair full-tensor change is 0.526%; regional nodal-potential and gradient-recovery estimates agree at roughly percent level. The exact axial harmonic P2 reference is reproduced at approximately 1e-14 relative Hessian error. Dense Gaussian flux diagnostic residuals are 0.382% (h=.10) and 0.051% (h=.075), with h=.10 finite-difference step comparison stable. RF Hessian, differential weighting gradient and Maxwell capacitance are also stable across these three cases. Full pytest suite: 26 passed.

**Decision:** numerical GO recommendation for a later, bounded ion-dynamics stage in this explicit ideal PEC/vacuum geometry only, after a user authorization and specification of a drive point. M2b.2 did not run trajectories or randomness. This does not validate an integrated silicon device: actual dielectric/material stack, finite conductivity, fabrication curvature, RF circuit, vacuum package, ion loading, lifetime and entropy remain unresolved. An asymptotic convergence/error certificate is still absent; .075 was not independently tested at a second angular quadrature order, and no separate center/end refinement campaign was run.

### M2b.2 files and resume point

- Report: `M2B2_DC_CONFINEMENT_VALIDATION_RU.md`
- Updated roadmap: `COMPUTATIONAL_ROADMAP.md`
- Gmsh generator: `scripts/run_m2b_3d_preflight.py` (corrected independent mesh-size controls)
- Mesh-only probes: `scripts/run_m2b2_dc_surface_mesh_probe.py`, `results/m2b2_effective_mesh_probe_h0p15/h0p1/h0p075.{json,npz}`
- P2 runs/fields: `scripts/run_m2b1_p2_physical_check.py`, `results/m2b2_p2_h0.75_hc0.18_dc0.15.json`, `results/m2b2_p2_effective_h0.75_hc0.18_dc0.1.json`, `results/m2b2_p2_effective_h0.75_hc0.18_dc0.075.json`, matching field NPZs
- Reference/crosschecks: `scripts/run_m2b2_dc_reference.py`, `scripts/analyze_m2b2_dc_saved_p2.py`, `scripts/run_m2b2_gaussian_flux.py`, corresponding `results/m2b2_dc_reference.json`, `results/m2b2_dc_saved_p2_crosscheck_h0p15_effective_h0p1_effective_h0p075.json`, `results/m2b2_gaussian_flux.json`
- Test command: `..\..\.venv\Scripts\python.exe -m pytest -q` (26 passed).

Next session should remain stopped until explicit authorization. If approved for dynamics, scope the first step narrowly: choose Ca+ operating RF/DC voltage and frequency, check q/a and Mathieu/Floquet margins from the validated fields, and then integrate a single deterministic ion with timestep/conservation convergence. Do not add stochastic baths or claim RNG entropy at that step. Preserve all old failed/interrupted outputs and do not modify adjacent impedance research.

## M3a completion update — 28 October 2026

M3a was authorized and completed locally. Exact inputs: M2b.2 `results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz` plus matching `results/m2b2_p2_effective_h0.75_hc0.18_dc0.075.json`. Geometry remains four finite hyperbolic RF branches, finite length 10r0, annular DC rings z=±5.5r0, grounded enclosure, empty vacuum, ideal conductors/polygonal boundaries. RF basis is 1-V differential terminal pattern; DC basis is 1 V on each ring. No analytic replacement was used for trajectories.

Reference point: Ca40+, r0=1um, VRF=10V differential peak, finite-field kRF=0.9998641 and |q|=0.2 gives fRF=2.4728311GHz; Hzz=1.879916e−4 V/r0^2 per 1V/ring gives VDC=0.0869788V per ring for the 1MHz local axial target. Full time-dependent FEM forces integrated for 20 RF periods (8.087ns), initial conditions near-center/radial/axial/combined and phase 0/π/2. Maximum excursion was 0.1532r0. No trajectory neared the selected absorbing cutoff 0.30r0. DOP853 and RK4, timestep tests, RF work/energy residual, ideal Mathieu monodromy convergence, and finite-field Floquet monodromy were saved. Floquet frequencies are radial 176.2509/176.2512MHz, axial 1.0000005MHz; multiplier modulus numerical error <1e−13. Deterministic pseudopotential minimum on a radius-.2r0 sphere is 0.325µeV along −z; this is not trap depth. Selected endpoint mesh sensitivity between DC surface h=.1/.075 is <1.5e−5r0.

Important limitations: pointwise RF gradient cross-mesh RMS difference is below 0.013% to R=.3r0; weak DC gradient relative difference is about 5–15% over R=.05–.3r0, with small absolute mismatch. P2 gradients jump at element faces. 20 RF periods cover only 0.008 axial secular period, so axial frequency is a local Floquet result, not a directly observed trajectory. No long-term heating/lifetime, escape saddle, gas collisions, cooling, stochastic source, entropy or RNG work was performed. Conditional deterministic GO only; stochastic/entropy stage is NO-GO until bounded multi-axial-period integration and force-recovery uncertainty are addressed.

### M3a files

- Report: `M3A_DETERMINISTIC_ION_DYNAMICS_RU.md`
- Configuration: `configs/m3a_ca40_reference.yaml`
- FEM point evaluator: `src/iontrap_rng/fem_fields.py`; exact polynomial test `tests/test_fem_fields.py`
- Scripts: `scripts/m3a_validate_field_region.py`, `scripts/m3a_mathieu_benchmark.py`, `scripts/run_m3a_deterministic.py`, `scripts/m3a_finalize_linear_modes.py`, `scripts/m3a_finalize_pseudopotential.py`, `scripts/m3a_mesh_trajectory_check.py`, `scripts/plot_m3a_results.py`
- Results: `results/m3a_field_region_audit.json`, `results/m3a_mathieu_benchmark.json`, `results/m3a_deterministic_results.json`, `results/m3a_trajectory_traces.npz`, `results/m3a_mesh_trajectory_check.json`, figure `results/figures/m3a_deterministic.png`
- Tests: 29 passed with `..\\..\\.venv\\Scripts\\python.exe -m pytest -q --basetemp results/pytest_m3a`.

Next step: stop and await user authorization. If a follow-on is authorized, restrict it first to a few-axial-period full-RF convergence test and smoother/recovered P2 force field, with explicit uncertainty; no stochastic dynamics, entropy estimation, RNG bit generation, remote computation or broad sweep until that gate passes and receives separate authorization.

## M3a.1 completion update — 28 October 2026

M3a.1 inspected the saved M3a inputs and all three effective M2b.2 P2 field sets. Geometry remained the four-branch hyperbolic RF structure with axial DC rings in an empty vacuum region and ideal-conductor boundaries. RF normalization is 1 V differential peak basis; DC basis is 1 V on each ring. M3a `VDC=0.08697882845 V/ring` was calibrated to 1 MHz using the h=.075 `Hzz`; it was not an independent mode prediction.

At the same fixed `VDC`, local Hessian frequencies are 997.220, 999.234 and 1,000.001 kHz for surface targets .15, .10 and .075r0. Direct DC-force mesh differences are 4.41% at x=+.1r0, 8.97% at x=−.1r0, 5.14% at z=+.1r0, 12.66% at z=−.1r0, 10.52% at z=+.2r0 and 3.96% at z=−.2r0. Fine-grid axial-force magnitudes at ±.1r0 differ by ~13.7%. Local potential fits have low potential residual but their gradients remain 4.02–6.50% RMS from elementwise P2 gradients inside .3r0; a reduced field is not validated for long runs.

**Decision:** do not run extended trajectories on current fields. This is a numerical stop, not a physical impossibility result. No 0.1/1/3/5 axial-period integration, new timestep study, escape result, stochastic force, entropy, or RNG calculation was performed. Recommend authorized independent central h-refinement with fixed physical DC voltage, force-symmetry checks, and direct/recovered-gradient comparison before reopening trajectory work. Local direct-P2 cost estimate is ~10 min/trajectory per axial period and ~50 min for five periods; these calculations were not launched.

### M3a.1 files

- Report: `M3A1_EXTENDED_DYNAMICS_VALIDATION_RU.md`
- Configuration and scripts: `configs/m3a1_force_gate.yaml`, `scripts/m3a1_dc_force_audit.py`, `scripts/m3a1_reduced_dc_field_audit.py`, `scripts/m3a1_field_cost_benchmark.py`, `scripts/plot_m3a1_dc_force_audit.py`
- Results and figure: `results/m3a1_dc_force_audit.json`, `results/m3a1_dc_force_audit.log`, `results/m3a1_reduced_dc_field_audit.json`, `results/m3a1_field_cost_benchmark.json`, `results/figures/m3a1_dc_force_audit.png`
- Test: `tests/test_m3a1_audit.py`; command `..\..\.venv\Scripts\python.exe -m pytest -q --basetemp results/pytest_m3a1` — 31 passed.

**Resume point:** wait for explicit authorization. First fix/characterize central field convergence; do not begin trajectories or stochastic/RNG modeling until the force gate passes. Preserve M3a artifacts and do not change the neighboring impedance project.

## M3a.2 completion update — 28 October 2026

The M3a.2 study reused the exact hcenter=.18r0 M2b.2 mesh/field (points and tetrahedra match exactly) and generated `.09/.06r0` central mesh probes before P2 solves. The actual local edge distributions became finer; DC surface `.075r0`, RF/end `.30r0`, geometry, tags, enclosure and fixed voltage were held constant. New `.09/.06` fields are saved; full fixed-coordinate vector forces, both error definitions, Hessian fits, gradient-based recovery, symmetry, point-to-face clearance and tag audits are in `results/m3a2_center_field_comparison.json`.

The median central edges inside `.15r0` were `.384/.202/.146r0`; P2 DOFs were 215,157/214,566/214,324. hcenter `.09→.06` force mismatch reached 10.55% vector and 4.78% magnitude across the axis sample set. At ±z=.1r0, the direct axial force magnitudes still differ by ~4.9%; potential parity residual is small, but gradient-parity RMS is 2.73%. The likely residual is unstructured mesh/P2 derivative asymmetry. Boundary tags are complete and symmetric by counts; no post-hoc symmetry was imposed. The `.06` P2 factorization took 112.2s; no 3D sweeps.

**Decision:** NO-GO for trajectory extension. Neither 0.1 nor 1 axial period was integrated. The M3a.1 vector-error table headings were corrected to distinguish vector difference from magnitude difference. No stochastic/entropy work was done.

### M3a.2 files

- Report: `M3A2_CENTRAL_FORCE_VALIDATION_RU.md`
- Config: `configs/m3a2_center_refinement.yaml`
- Scripts: `scripts/m3a2_probe_center_meshes.py`, `scripts/m3a2_solve_center_mesh.py`, `scripts/m3a2_analyze_center_fields.py`, `scripts/plot_m3a2_center_convergence.py`
- Mesh probes: `results/m3a2_center_mesh_hc0p18.json/.npz`, `...hc0p09.json/.npz`, `...hc0p06.json/.npz`; saved P2 fields: `results/m3a2_center_fields_hc0p09.npz` and `results/m3a2_center_fields_hc0p06.{npz,json}`
- Analysis/figures: `results/m3a2_center_field_comparison.{json,log}`, `results/figures/m3a2_center_force_convergence.png`, `results/figures/m3a2_axial_symmetry.png`

**Historical M3a.2 resume note — superseded by M3a.3 below.** At that point the next task was to improve center mesh quality/symmetry and force convergence; M3a.3 has since implemented and audited a symmetric octant, but the ion-dynamics gate remains NO-GO.

## M3a.3 completion update — 28 October 2026

**Decision: NO-GO for extending deterministic trajectories.** The physical CAD and paired RF/symmetric DC excitation patterns were verified invariant under x, y, and z reflection. The positive octant was solved with homogeneous Neumann conditions only on open vacuum sections of x=0/y=0/z=0; conductor cut faces were given terminal Dirichlet values, and the outer enclosure was grounded. All external boundary facets were classified on both solved and the mesh-only finest candidates. Individual weighting-potential problems were not assumed even.

The local sizing field was changed to have a central plateau plus smooth transition while retaining separately fixed DC/RF electrode surface, rail-end and far-field controls. The `.030`, `.018`, and `.014r0` center targets yield measured local tetrahedral edge median/p90 inside `.15r0` of `.03298/.04852`, `.02184/.03002`, and `.01787/.02344r0`, respectively. P2 DOFs are 33,014, 45,824 and 63,759. Only `.030` and `.018` received P2 solves. DC ring surface edges were unchanged between meshes: 1,346 unique local edges, median `.09798r0`, p90 `.12930r0`; thus central refinement is independently controlled.

At fixed `VDC=.08697882845042744 V` per ring and `VRF=10 V` differential peak, DC force change between the two solved octants reaches 2.512% vector/magnitude and 2.088% over the generic off-axis points. Fine full-domain versus fine octant DC force differs up to 5.243% vector and 3.416% magnitude. Full-domain DC potential reflection mismatch is at most `8.34e-8 V` per 1-V basis at selected reflected pairs, while direct P2 gradient mirror residual reaches 5.09%, pointing to derivative/postprocessing plus unstructured mesh effects; the precise split is unresolved. RF Hessian and RF field on relevant generic off-axis points are stable; relative RF metrics at near-null axial points are denominator-dominated.

Unconstrained octant DC quadratic fit at R=.35r0 gives normalized-basis Hessian `(-9.3488e-5,-9.4572e-5,+1.8745e-4) V/r0²`; from `.030→.018`, `Hzz` shifts 1.65%. Order-6 harmonic polynomial fits from fine octant potentials agree with direct P2 gradients on six held-out points within 0.096% DC and below `4e-5%` RF, but the fit enforces Laplace and shares the underlying FEM potential, so it is a smooth local representation test, not independent PDE validation.

Mesh-quality audit still finds central low-quality tetrahedra (quality <.1: 0/1/2 in `.030/.018/.014` meshes; `.014` central minimum mean-ratio `.0297`). Available physical RAM was ~1.28 GiB. The `.014` P2 factorization was not attempted; a simple fill estimate yields roughly 50–80 million LU nonzeros before peak work buffers. No trajectories were run. The octant gives substantial computational savings and can support adaptive P2, but the force-convergence gate remains failed.

### M3a.3 files and reproducibility

- Report: `M3A3_SYMMETRIC_FEM_VALIDATION_RU.md`
- Configuration: `configs/m3a3_symmetry_refinement.yaml`
- Scripts: `scripts/m3a3_symmetry_audit.py`, `scripts/m3a3_octant_mesh.py`, `scripts/m3a3_octant_p2.py`, `scripts/m3a3_compare_fields.py`, `scripts/m3a3_octant_recovery.py`, `scripts/plot_m3a3_diagnostics.py`
- Meshes/stats: `results/m3a3_octant_mesh_plateau030.{npz,json}`, `plateau018.{npz,json}`, `plateau014.{npz,json}`
- Boundary-tag audits: `results/m3a3_octant_boundary_plateau030.json`, `...plateau018.json`, `...plateau014.json`
- Solved fields: `results/m3a3_octant_fields_plateau030.npz`, `...plateau018.npz`
- Comparisons: `results/m3a3_full_octant_force_comparison.json`, `results/m3a3_octant_recovery.json`, symmetry audit `results/m3a3_symmetry_audit.json`; figure `results/m3a3_convergence.png`
- Code validation: six scripts pass `py_compile`; existing project suite is **33 passed** with `pytest -q --basetemp results/pytest_m3a3`.
- The initial `.030` transition mesh without central plateau and the early failed RF-end surface selection remain preserved under their earlier filenames/console history; use the `plateau*` results for this report.

## Stage R0 update — 29 October 2026 (latest)

### Scope and status

R0 prepared a candidate public scientific project, a Windows FEM workflow, bounded species comparison, and effects register. No remote FEM, new FEM solve, stochastic trajectory, entropy estimate or RNG bitstream was run. Current physics gate remains the M3a.3 numerical **NO-GO for extended ion dynamics**: resolved `.030`/`.018` P2 octant DC forces differ up to 2.512%; fine full-domain vs octant differs up to 5.243%. `.014` is a generated mesh only, not a solved field.

### Added R0 assets

- `README.md`, `.gitignore`, `requirements-lock.txt` (exact top-level versions from local Python 3.12.10; transitive/platform lock is incomplete).
- `src/iontrap_rng/species_screen.py`, `scripts/compare_species.py`, `tests/test_species_adaptive.py`; generated result: `results/r0_species_comparison.json`.
- `src/iontrap_rng/adaptive.py`: residual/flux-jump estimator assembly and Dörfler marking helper. Not a solver adapter or independent FEM error certificate.
- `PARTICLE_SPECIES_COMPARISON_RU.md`, `PHYSICAL_EFFECTS_REGISTER_RU.md`, `PHYSICAL_EFFECTS_REGISTER.json`.
- `REMOTE_FEM_EXECUTION_GUIDE_RU.md`, `configs/r0_remote_fem.json`, five scripts in `remote/`.
- `scripts/audit_publication.py`; candidate manifest `PUBLICATION_MANIFEST.json` is generated after final edits and excludes the raw literature search response dump.

### Analytical results and evidence

At ideal fixed (|q|=.2), 10V differential peak, (r_0=1\,\mu m), script RF frequencies are 2.473GHz Ca+, 1.667GHz Sr+, 1.196GHz Yb+, 5.208GHz Be+, 667.5GHz electron. All give the same idealized 0.25eV pseudopotential at (r=r_0). At fixed 100MHz, q-preserving voltage/depth are 16.35mV/0.409meV (Ca+) and 0.2245µV/5.61neV (electron). These are not finite geometry escape barriers. The 1µm plane image-charge estimate −0.36meV flags electron surface sensitivity but is not a complete potential. Conditional (g_w=.4/r_0) current estimates are not validated 3D receiver predictions. References include Taniguchi & Noguchi (2025), DOI 10.1103/6mpr-5s2h; NIST isotope tables; Brownnutt et al. DOI 10.1103/RevModPhys.87.1419; Flowers-Jacobs et al. DOI 10.1088/1681-7575/aa7b3f.

### Publication status / repository boundary

At the end of Stage R0, no `.git` existed at the project root or parents and GitHub CLI auth was unavailable. R0.1 resolved this by creating and pushing the standalone local project only. Project root remains separate from sibling/parent impedance materials. Publication audit excludes cache, `.venv`, pytest temp, and `LITERATURE_SESSION_NOTES.json` (unassessed redistribution rights of embedded search snippets); curated literature audit is kept. No license file or license grant included.

### Remote workflow

Target host Ryzen 5 7430U / ~16GiB, possibly shared. Config requests one sequential solve from existing `plateau014` octant mesh, single worker, at least 13GiB free RAM and 10GiB disk. This threshold is an operational guard, not a peak-memory guarantee. 50–80M LU entries for `.014` are extrapolated, not measured. No remote action initiated. Solver has no checkpoint within LU; run script defaults to dry-run and requires `-Execute` plus prior user approval. Runtime polling can stop only its own Python process at configured RSS/free-RAM thresholds; this is not hard memory isolation and a burst can outrun polling. P2 SPD CG+AMG is a possible lower-fill method but must first match SuperLU force/Hessian/flux on `.018` with true residual and measured memory; untested here. After approval run setup/test, preflight/dry-run; only one solve, then inspect actual fill and force-functional convergence.

### Resume instructions

1. Complete `python scripts/audit_publication.py`, inspect manifest/file list and verify excluded content.
2. Run `..\\..\\.venv\\Scripts\\python.exe -m pytest -q --basetemp results\\pytest_r0` and `python scripts/compare_species.py`; generated JSON is `results/r0_species_comparison.json`.
3. R0.1 created the public repository and verified a fresh clone; see the latest R0.1 section below for URL and SHA status.
4. Do not execute remote FEM until explicit authorization; after approval run `remote/setup_windows.ps1`, preflight and dry run, then at most one `.014` P2 solve if resource gates pass.
5. Do not start deterministic ion trajectories until DC force convergence is repaired; no stochastic/entropy/RNG work was authorized by R0.

### Validation and versioning details

The Stage R0 local regression suite passed **37 tests** using Python 3.12.10. Four PowerShell scripts parsed successfully. An archive snapshot also passed 37 tests before publication. During R0.1, an actual public `git clone` succeeded, its tracked inventory matched the local project (196 files), and all 37 tests passed in the clone. The remote initial SHA matched the original local commit. Git author metadata uses the connected account's GitHub noreply form; personal email is not in the commit object.

**Next step only after new user authorization:** one bounded adaptive P2 octant refinement with an independent force-functional/error estimator and improved central element quality. If that still cannot resolve DC-force uncertainty, separately scope a P3 or boundary-element benchmark. Do not restart trajectory dynamics or stochastic/entropy/RNG modeling until the PDE force gate and stability-margin uncertainty are acceptable. Preserve the M3a.2 and M3a.3 artifacts; do not modify the adjacent impedance research.

## Stage R0.1 — public publication completed (29 October 2026)

- Public URL: https://github.com/snejokzboxik/ion-trap-rng-feasibility
- Visibility verified as `PUBLIC` by `gh repo view`.
- Initial pushed scientific commit: `aa73cfb439f3f439bcbf2087ff983c53162c8b28`; remote `main` initially matched local HEAD.
- Fresh public clone: `research/_r0_1_public_clone_verify_20261029`; clone inventory exactly matched local tracked inventory (196 files), and 37/37 tests passed in the clone using the available Python 3.12.10 dependency environment. Clone working tree was clean.
- Publication manifest: 195 non-self files, full SHA-256 match; 0 secret-scan findings; no tracked file >50MB, largest 5,711,193 bytes; no PDF files. `LITERATURE_SESSION_NOTES.json` remains untracked/ignored, excluded for unresolved redistribution rights. No license file/grant added.
- A follow-up commit updates README, manifest generator/manifest, roadmap, summary, remote guide, and this handoff to reflect the publication. This follow-up was pushed and checked so local and remote `main` match at the end of R0.1.
- No FEM, remote compute, trajectories, stochastic dynamics, entropy extraction or RNG bit generation was run during R0.1.

For remote Windows setup, clone the public URL above, then from the clone root run `./remote/setup_windows.ps1 -DryRun`, `./remote/setup_windows.ps1`, `./remote/preflight.ps1 -DryRun`, and `./remote/run_fem.ps1 -DryRun`. A real FEM launch remains separately gated and was not authorized/launched in R0.1.

## Stage R0.2 — remote FEM workflow repair (29 October 2026)

The first published R0 remote workflow was incomplete: `run_fem.ps1` invoked a `preflight.ps1` that was absent from the published tree. R0.2 adds a fail-closed preflight and repairs end-to-end run and result verification. The mesh and configuration identities are SHA-256 checked; `.014` mesh coordinates are stored as `(N,3)` and tetrahedra as `(Nt,4)`, with the P2 DOF count independently reconstructed from unique edges. The solver now writes a distinct solve boundary audit so it cannot overwrite the preserved mesh-only audit.

`remote/preflight.ps1` enforces the configuration thresholds: 13 GiB minimum free physical RAM, 10 GiB minimum free disk, a nominal 5 GiB OS/other-workload reserve, 8 GiB maximum solver working set, one thread and one solve at a time. It inspects Windows/Python/package versions, config integrity, mesh/audit identity, CPU load, existing solver processes, and output write access. `remote/run_fem.ps1` defaults to safe dry-run and starts no solver unless `-Execute` is passed; it reruns real preflight before launch, refuses output overwrite, preserves logs/partials, and terminates only its own child through its process handle if the memory guard trips. It does not manage the separate Drosophila workload.

`remote/verify_results.ps1` calls `scripts/verify_remote_outputs.py`, which validates output readability, array finiteness and dimensions, mesh identity, P2 DOFs, audit status, and successful solver completion. This is an artifact-integrity gate, not a force-convergence certificate. A fixture test suite plus PowerShell `-TestMode` failure scenarios exercise success, missing prerequisites, insufficient RAM/disk, missing mesh, invalid config, failed solve, and incomplete outputs without starting FEM.

**R0.2 publication and validation completed:** repair commit `06af7aab3808625af7444af989b46a94a4958a95` was pushed to existing public `main`; verified repository visibility is `PUBLIC`, and remote `main` matched this commit. A genuine fresh clone was created at `research/_r0_2_public_clone_verify_20261029`; it had 201 tracked files, clean status, and all 45/45 Python regression tests passed. Local suite also passed 45/45. The PowerShell validation harness passed all 13 simulated host/runner scenarios, including normal safe dry-run, missing Python/dependency, low RAM/disk, absent mesh, invalid configuration, high CPU, existing solver, unwritable results, simulated solver failure and incomplete output. All PowerShell scripts passed parser validation. Publication audit includes 200 non-self files (52,165,150 bytes), reports zero secret-pattern findings, and continues to exclude `LITERATURE_SESSION_NOTES.json` because redistribution rights for raw search snippets are unresolved. The actual sandbox-host preflight was attempted in dry-run only and failed closed because CIM access is denied by this local execution sandbox; this is not a NINKEAR result. No FEM solve, remote computation, particle dynamics or RNG simulation was launched.

**Readiness decision:** the repaired Windows execution path is code- and test-validated but is not yet authorized or operationally cleared for an unattended overnight solve. Before any first run, NINKEAR must pass real setup/regression, `preflight.ps1 -DryRun`, and `run_fem.ps1 -DryRun` while the other workload is in its intended operating state. Keep the configured 13 GiB free RAM, 10 GiB free disk, 5 GiB reserve, 8 GiB worker limit and one-solve restriction; if they fail, do not lower them. The next authorized action can be the actual NINKEAR dry-run/preflight review, followed by a separate decision on whether to authorize one bounded `.014` solve.
