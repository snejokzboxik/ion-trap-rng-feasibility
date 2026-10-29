# Вычислительная дорожная карта: RF ion-trap RNG

**Статус (после M3a.2, 28 Oct 2026):** центральное уточнение `.18→.09→.06r0` реально уменьшило локальный mesh size и улучшило зеркальность potential/gradient. Тем не менее прямые DC силы между hcenter=.09/.06 сетками расходятся до 10.55% по вектору и 4.78% по модулю; осевая зеркальная ошибка остаётся около 4.9% при ±.1r0. Поэтому ни 0.1, ни один осевой период не интегрировались. Текущий numerical gate остаётся **NO-GO для trajectory extension и stochastic/entropy работ**.

## Предлагаемый benchmark полного устройства

Сравнительные scenarios: combinations useful output 10 kbit/s, 100 kbit/s, 1 Mbit/s and continuous power 0.1, 1, 10 mW; die/package/startup/survival defined per stated workload. Это comparator points, не физические пределы. Например, 1 mW/100 kbit/s = 10 nJ/bit; это C-level benchmark, не запрет P-level feasibility.

## Этапы и go/no-go

| Stage | Работа | Scope | Gate |
|---|---|---|---|
| 0. Требования и источники | Заморозить operating envelope и uncertainty | Source registry, без solver | Каждый critical параметр имеет источник/range; иначе только conditional conclusion |
| 1. Аналитика (M1, выполнено) | Геометрия, Mathieu/Floquet, q–frequency–voltage–depth–field, axial DC | Python package + 12 unit tests | Analytical references and unit checks pass; this is not RNG proof |
| 1.5. Early observability/entropy/system screen (предварительно выполнен) | Readout, noise/collision heating, cooling, lifetime and energy | Baseline + explicit sensitivity assumptions; unknowns remain unknown | At least one conditional interval survives; no entropy claim without supported source model |
| 2a. Bounded 2D readout (выполнено) | Ideal manufactured Laplace test + finite four-hyperbola terminals, weighting potentials, C' and Ramo current | 5 local h levels r0/20–r0/320, geometry sampling 50/100/200; P1 triangles | Finite weighting gradient changes 0.060% from /160 to /320; electrode C' changes 0.019%. Corrected global ideal E-L2 is ~2.46%; local high multipoles not precision-validated |
| 2b. 3D sense/electrode (bounded preflight выполнен; field validation failed) | Manufactured P1/P2 check; one idealized vacuum geometry with axial rings, terminal weighting fields and C matrix | Local tetrahedral FEM only; no qualified materials, port, or full-wave model | Do not advance to dynamics until Hessian/multipole, curved-boundary and energy/flux checks pass and a source-grounded support/dielectric stack is specified |
| 3. Finite 3D trap | Truncated hyperbolic electrodes, axial DC, rounded edges | Baseline и symmetry control | Multipole residual ≤1% при \(\rho≤0.5r_0\), converged saddle и mode stability |
| 4. EM/material/thermal | Conductivity, contacts, dielectric, RF circuit, heat | EQS first; selected EM check | Model change ≤5% в decision observables; margins robust |
| 5. Динамика/lifetime | Full RF, collisions, stray fields, shutdown/reload | Selected scenarios | Timestep convergence, survival CI проходит target |
| 6. Readout/entropy | Weighting fields, circuit/ADC, backaction, controls | M2a has field and lumped-readout screens; full entropy remains future work | Conditional entropy lower bound поддерживает useful rate; M2a SNR is assumption-conditional and not entropy evidence |
| 7. Loading/package | Micro-oven/ionizer, getter, cavity, startup energy | Один module budget + sensitivity | Startup, availability, area, volume, power проходят benchmark |
| 8. Feasibility map | Adaptive, uncertainty-aware map | Только после validation | Robust feasible region есть; иначе публикуется dominant no-go |

Эти tolerances относятся к verification солвера, не к универсальным physical requirements. 20 \(k_BT\) — sensitivity screen, не lifetime law. \(10^7\) V/m не является допустимым silicon field limit. Pressure определяется species/gas-specific rates and loss probabilities.

## Три независимых уровня успеха и классификация границ

- **P (physical):** stable ion plus physically supported stochastic source, measurable ion-mediated component and defensible conditional entropy/lifetime model.
- **E (engineering):** integrated loading/restart, cavity, process flow, thermal/power, package and reliability within declared system boundary.
- **C (competitive):** common-boundary comparison to CMOS TRNG by useful entropy rate, joule/bit, footprint and reliability/security overhead.

Каждую min-size boundary помечать F (fundamental при оговорённых предположениях), T (technology-dependent), R (application requirement-dependent) или U (model/data unresolved), указывать uncertainty interval/source. Доминирующий U оставляет вывод conditional; он не превращается в universal minimum.

## M1.5 equations / screening point

Image-current scale $I\sim |Q|\omega A/d_{\mathrm{eff}}$ uses $E_w\sim1/d_{\mathrm{eff}}$ pending an actual weighting-field solve; Johnson current screen $i_n=\sqrt{4k_BT B/R}$. For weak one-sided field PSD, $\dot n=Q^2S_E/(4m\hbar\omega)$, requiring active cooling/removal to balance absorbed noise/RF power in a stationary source. Gas rate $\Gamma=pk/(k_BT_g)$; ideal sparse independent Poisson event-indicator entropy ceiling is $\Gamma/\ln2$ bit/s only. These are model-conditional scales, not device output.

Current Ca+ point (r0=1 µm, q=.2, Vdiff=10 V): RF ~2.47 GHz and radial secular ~175 MHz; A=100 nm and deff=1 µm yield ~17.6 pA idealized image-current scale. A 1 MΩ, 300 K, 1 kHz ENBW gives 4.07 pA rms Johnson current, but weighting cancellation, impedance transfer, RF pickup and GHz front-end remain unknown. Literature noise amplitudes span about nine decades across other trap geometries/conditions and are not directly transferable to micron silicon at 175 MHz.

Для бинарных collision event/no-event samples ENBW 1 kHz даёт не более 2 k symbols/s. Условный 100 kbit/s binary stream требует минимум 50 kHz ENBW; при том же 1 MΩ thermal current noise составляет 28.8 pA rms, выше nominal 17.6 pA image-current scale. Для Poisson event model этого также требует collision rate ~6.93×10^4 s⁻¹, pressure ~0.096–0.96 Pa при assumed $k=0.3$–$3\times10^{-15}$ m³/s. Это показывает competing requirements, но не bound для multi-bit analog readout.

M2a resolved a bounded weighting-field question and tested solver/mesh convergence. The paired-terminal center gradient cancels to 4.39e-4 of the one-terminal gradient; a normalized differential port gives 0.39974/r0. At the Ca+ baseline, assumed A=100 nm implies ~7.03 pA peak at 174.78 MHz. With assumed 1 pF parasitic and stated amplifier noise, 1 MΩ parallel RC gives ~0.10 peak/rms at 1 kHz ENBW. A high-Q parallel resonator can raise this conditional ratio but requires uncertain input capacitance, loaded Q and inductance; feedthrough and front-end power are unknown. These screens supersede the scalar effective-distance image-current estimate for that selected geometry, but they still do not establish measurement feasibility.

**M2a gate result (historical):** at that time only a conditional 3D model was recommended and no authorization existed. This was superseded by the explicit M2b authorization. Conventional paired RF terminal direct sensing remains a no-go at the ideal symmetric center. 2D electrostatics cannot resolve surface spectra, collision lifetime, front end, loading or total power.

## M2a.1 — audit readout and mesh accuracy (bounded, complete)

- **Mesh finding:** the original centroid-only error metric underestimated the exact element-integrated P1 field L2 error by omitting within-element variation of the analytical linear field. Corrected local-series error plateaus near 2.46%; lowering local hmin from r0/40 to r0/160 changes 2.464%→2.456%. Uniform ideal reference gives 3.319%, 1.664%, 0.833% at h/r0=0.4, 0.2, 0.1. Do not claim <0.1%; no P2 FEM or rigorous a posteriori estimator has been implemented. The quadratic manufactured ideal solution is exactly representable by P2 in theory, which is a future validation case, not a completed P2 result.
- **Physical enclosure:** with hmin/r0=1/80, moving the grounded outer square from 4r0 to 10r0 changed the local differential weighting gradient by <0.003%, but the 2D per-electrode capacitance by about 7.8%. Package geometry and dielectric stack remain unresolved.
- **Receiver evidence:** published anchors span 7 fF input capacitance for a 65-nm sensor buffer and 1.02 pF total input for a separate CMOS capacitance sensor; these are dissimilar devices, not a design range for the ion trap. Peer-reviewed cryogenic Penning receivers report Q up to 3500 near 39 MHz/~1 mW; a 2025 room-temperature Paul electron image-current experiment used separate 619-MHz readout and 3.105-GHz trapping cavity modes, with hundreds of mW input in its present cavity coupling. These are topology/parameter analogues, not proof of CMOS Ca+ readout.
- **Candidate:** distinct readout cavity/mode plus balanced differential pickup is the one conditional architecture for a future port model. At constant q,Vd,A/r0=0.1, Ca+ analytical scaling across r0=10 µm…100 nm gives fsec=17.5 MHz…1.75 GHz, trap RF=247 MHz…24.7 GHz, current=0.703…70.3 pA and ideal 2D trap capacitance=3.75…0.0375 fF for Lz=10r0. These are reference-model scalings only; material/drive losses and actual motion are unknown.
- **Historical next gate:** at M2a.1 close, the above 3D model awaited authorization and source-bounded stack/receiver parameters. The later M2b authorization superseded that stop and produced the idealized preflight documented below. Its results fail Hessian convergence; material/receiver feasibility remains open.

Reproducible M2a.1 outputs: `results/m2a1_numeric_audit.json`, `results/m2a1_scale_screen.json`; scripts: `scripts/run_m2a1_bounded_audit.py`, `scripts/run_m2a1_scale_screen.py`; report: `M2A1_READOUT_REVIEW_RU.md`.

At r0=10 nm and fixed backing ratio 0.1, silicon thickness would be 1 nm. That marks an unresolved/technology-sensitive continuum, doping, oxide and interface problem, not a mesh-resolution result or universal physical size boundary. The r0/320 element size is numerical; it does not assert an atomically meaningful material mesh.

## M2b — historical bounded 3D preflight (superseded by M2b.1 audit)

Manufactured square Φ=(x²−y²)/2, ([-6,6]^2r_0), exact Dirichlet data, P1/P2, 12/24/48 cells per side: P1 global E-L2 = 8.333/4.167/2.083%, local ρ<r0 = 56.87/28.55/14.34%; P2 errors for potential/field/Hessian are at numerical roundoff (<(3.3\cdot10^{-13}) relative). The independent exact P1 element-moment calculation matches FE quadrature to <(4\cdot10^{-17}). Linear polygon sampling sensitivity from M2a remains a shape sampling check, not curved-boundary validation.

The authorized 3D run is an idealized **vacuum/perfect-conductor preflight**, because published work did not ground the exact four-branch hyperbolic silicon/support/CMOS stack. It includes four finite hyperbolic branches, L=10r0, analyst-assumed backing thickness (0.1r_0), two annular DC rings at z=±5.5r0, an enclosure, six individual terminal weighting solutions, and a grounded-reference Maxwell matrix. P1 tetrahedral meshes use ((h_{max},h_c)/r_0=(.90,.25),(.75,.18),(.65,.13)). Finest (C_{diff}approx5.36\cdot10^{-16}) F and (|\nabla\phi_{w,d}|approx0.4136/r_0) V(^{-1}).

**3D validation is not passed.** Between the two finest grids the weighting-gradient magnitude changes ~0.21%, Cdiff ~2.42%, and matrix norm ~3.4%; RF Hessian components change 9–82%, axial DC Hessian components 24–60%, and the RF Hessian trace residual is ~9.9%. These are grid deltas, not certified errors. The coarsest center fit is under-resolved. Consequently no ion dynamics or entropy use of these fields is justified. No dielectric distortion, material loss/heating, receiver transfer, RF mode or feedthrough was calculated. See `M2B_3D_FEASIBILITY_REPORT_RU.md`.

### Gate before dynamics

1. Add automated rank/conditioning and minimum-sample checks to central polynomial fitting.
2. Refine the one declared geometry locally at the ion, curved electrode surfaces, and z ends; compare P1/P2 or an independent solver and curved/isoparametric boundaries.
3. Cross-check the Maxwell C matrix with surface charge/flux and reciprocity; verify Laplace residual, RF trace/Hessian and symmetry.
4. Establish a source-supported electrode/support/dielectric/ground stack; if impossible, retain explicit preflight status and do not claim engineering feasibility.
5. Only after the 3D Hessian/multipole gate passes may dynamics begin. Entropy requires a separate physical stochastic source and receiver/cooling budget.

Files: `configs/m2b_validation.yaml`, `configs/m2b_preflight.yaml`, `scripts/run_m2b_manufactured.py`, `scripts/run_m2b_3d_preflight.py`, `scripts/plot_m2b_3d_preflight.py`, `results/m2b_manufactured_validation.json`, `results/m2b_3d_preflight.json`, `results/m2b_3d_fields.npz`, `results/figures/m2b_*.png`, `tests/test_m2b_manufactured.py`, `tests/test_m2b_preflight_results.py`.

## M2b.1 — Hessian diagnosis and bounded P2 recovery (complete; dynamics gate remains closed)

- The archived P1 10% RF trace came from a local quadratic fit over 23 nodes (`R=.5r0`, normalized condition number 12.96, relative fit residual 18.46%). P1 tetrahedra have zero elementwise Hessian, so their “Hessian” is necessarily recovered/postprocessed. Radius, weights, and fit order changed the result. Harmonic-constrained fits force zero trace by construction and are not evidence of FEM harmonicity.
- Added `src/iontrap_rng/hessian.py` for scaled rank-aware QR/SVD least squares and P1 gradient recovery; `scripts/run_m2b1_hessian_diagnostics.py` preserves old-grid diagnostics separately.
- Independent 3D quadratic harmonic test: P2 reproduces potential, electric field and Hessian at roundoff on 2/4/8 structured cells per axis. P1 field L2 error reduces 54.25%→27.12%→13.56%; its element Hessian is undefined as a smooth observable. This validates the solver/extraction path only for a quadratic exact reference.
- Two finite-geometry P2 levels (95,068 and 116,097 DOFs) produce central RF tensor change 0.0173%, and RF trace fractions 2.03e−5 and 1.70e−5. This is encouraging for the local RF quadrupole, but both global and center sizes changed together; no surface/end independent refinement or curved P2 CAD is represented.
- Axial DC Hessian changes 97.2% in Frobenius norm between those P2 levels, and DC gradient-fit residual remains ~0.11–0.18 around the center. This fails the requested convergence gate. Maxwell matrix change is 0.938%, and differential capacitance changes 0.90%; only two levels and ideal conductors. Differential weighting-gradient estimate changes ~0.39%, preliminary.
- Old P1 sphere-flux indicators are 0.43–9.52% depending on field/radius using 288 angular points; quadrature is not converged, so these remain diagnostic only. P2 sphere flux, strong-form residual, multipoles and independent full-domain symmetry remain undone. No universal FEM uncertainty interval can yet be assigned.

**Decision:** numerical no-go for ion dynamics, not a physical no-go for the trap. Retain ideal RF electrostatics as a conditional local reference. Next authorized work should independently refine central volume, hyperbolic surfaces, ring surfaces and axial ends; compute nodal-potential and gradient-recovery Hessians, elementwise P2 Laplacian and converged Gaussian flux; verify reflection symmetry, higher multipoles, weighting gradient and C matrix. Continue dynamics only after key observables change <1% on successive independent refinements and independent error checks agree. No source-qualified dielectric stack exists, so do not infer CMOS-ready integration.

Report: `M2B1_3D_VALIDATION_REPORT_RU.md`. Reproduction artifacts: `scripts/run_m2b1_manufactured3d.py`, `scripts/run_m2b1_hessian_diagnostics.py`, `scripts/run_m2b1_p2_physical_check.py`, `src/iontrap_rng/hessian.py`, `tests/test_m2b1_hessian.py`, `results/m2b1_manufactured3d.json`, `results/m2b1_old_grid_diagnostics.json`, and `results/m2b1_p2_h*_hc*_v2.json`. Full test suite: 26 passed.

## Python layout proposal

~~~text
research/ion_trap_rng/
  configs/ data/ results/ tests/ scripts/
  configs/m2a_readout.yaml
  src/iontrap_rng/  # M1: geometry.py, trap.py, floquet.py, scaling.py, screening.py, ions.py, constants.py
  scripts/run_m1_reference.py
  tests/test_m1_analytic.py
  results/m1_reference.json
  results/figures/m1_analytical_screen.png
  src/iontrap_rng/readout_fem.py  # M2a finite 2D Laplace / weighting-field / lumped port screens
  scripts/run_m2a_readout.py
  tests/test_m2a_readout.py
  results/m2a_readout.json
  # 3D FEM, material/RF loss, stochastic dynamics and entropy are future work
~~~

Candidate stack: Python 3.12, NumPy/SciPy, pytest, Gmsh, scikit-fem. Проверить versions, Windows support, licenses and solver limits before pinning or installing. Use SI internally, immutable hashed configs and complete manifests. Do not change the root .venv/requirements belonging to adjacent impedance work. Full Maxwell should use a vetted solver, not a new production solver built from scratch.

## Required validation

- Geometry: four separate nets, hyperbola residual, normals, clearances, no intersections.
- Fields: manufactured harmonic solution vs physical curved-boundary geometry; full-domain vs symmetry; exterior-boundary, h/p, fit-radius and Hessian convergence.
- Dynamics: Mathieu/Floquet reference, RF phase, timestep at \(T_{\rm RF}/16,32,64\), energy/work balance.
- EQS/circuit: analytic capacitor and RC slab, reaction vs surface flux vs energy, power balance, skin/displacement checks.
- Thermal: analytic slab and heat conservation.
- Stochastic: Langevin FDT/equipartition; collision kinematics/rate; PSD convergence; lifetime rare-event benchmark.
- Readout: analytic weighting fields, paired-electrode cancellation, current conservation, RLC/Johnson reference, feedthrough and ADC.
- Entropy: biased IID, Markov, periodic, colored-Gaussian, deterministic controls; ion-present/absent and bath-only comparison.
- System: wall-power accounting, startup/reload energy, package boundary and false-ready rate.

Proceed from selected validated points to any adaptive map. Convergence only supports the stated model; it does not experimentally validate unknown silicon noise, hermeticity, loading reliability or security.



## M2b.2 — status and gate update (28 October 2026)

**This section supersedes the earlier M2b.1 recommendation that the DC Hessian was unresolved and that the geometry could not enter ion dynamics.** Inspection found `Mesh.MeshSizeMin=0.18r0` clamped all requested DC surface sizes below 0.18r0. The generator now sets the global floor to the smallest requested local target and uses independent center, RF-surface, DC-ring-surface, axial-end and far-field controls. Historical clamped probes and P2 runs remain in `results/` and must not be labelled effective DC refinement.

Mesh-only probes verified monotonic DC triangle counts and shrinking local surface edge distributions for targets 0.15/0.10/0.075r0, while center h=0.18r0 and RF/end h=0.30r0 stayed fixed. P2 DC Hessian tensor change between the last two levels is 0.526%; independent regional potential-value fits and gradient recovery agree at approximately the percent scale. The exact axial harmonic manufactured solution passes to roundoff; free regional Hessian trace and Gaussian-sphere flux cancellation are small. The local validation suite remains 26 passed.

**Numerical recommendation:** GO only for a bounded, explicitly authorized ion-dynamics stage using this ideal PEC/vacuum geometry and declared drive voltages. This is not an engineering feasibility result. No trajectories were run in M2b.2. Further geometry/stack claims remain limited by polygonal boundaries, incomplete independent refinement of the center and electrode ends, absent silicon/dielectric specifications, no independent PDE solver, and no demonstrated complete device. See `M2B2_DC_CONFINEMENT_VALIDATION_RU.md` for tables and exact artifacts.

### Resume rules

1. Preserve the M2b.1 and M2b.2 result files, including ineffective and interrupted cases. Use only the `m2b2_effective_mesh_probe_*` and matching `m2b2_p2_*effective*` files for the repaired DC-refinement comparison.
2. M2b.2 did not run particle trajectories, collisions, stochastic noise, entropy or RNG output. Start none of these without explicit authorization.
3. If dynamics is later authorized, first fix a physical RF/DC operating point and ion species; independently verify stability against timestep/phase and interpolation error. Keep it a limited model; no parameter campaign.
4. Keep semiconductor integration classified as unverified. The present FEM uses an empty vacuum domain and ideal conductor boundary values; it does not model doped silicon, oxide, package, surface charge or RF loss.

## M3a — bounded deterministic reference case (complete; 28 Oct 2026)

- Inputs: saved M2b.2 `hDC=0.075r0` P2 FEM fields for the unchanged four hyperbolic RF branches, annular DC electrodes, finite length `10r0`, empty vacuum, PEC walls. Field amplitudes mean 1 V differential RF pattern and 1 V on each DC ring; spatial coordinate is `r/r0`, gradient converted by `1/r0`.
- Reference: (^{40}\mathrm{Ca}^+), `r0=1µm`, `VRF=10V` differential peak, finite Hessian-calibrated `fRF=2.472831GHz` for `|q|=0.2`, `VDC=0.086979V` on each ring for `fz,target=1MHz`.
- Validation: exact P2 quadratic gradient unit test; ideal Mathieu monodromy RK4 convergence; finite-field symmetrized-Hessian Floquet modes, DOP853/RK4 comparisons, timestep sensitivity, RF work balance. 29 tests pass, including M3a result checks.
- Short trajectories: four initial displacements and phases 0/π/2, 20 RF cycles (8.087ns), max radius `0.1532r0`; two selected trajectories compared across DC surface h=.10/.075 and endpoint differences <1.5e−5r0. Figure `results/figures/m3a_deterministic.png`.
- Linear modes: radial `176.2509/176.2512MHz`, axial `1.0000005MHz`; finite-field Floquet multipliers have modulus error <1e−13. Ideal q=.2 benchmark is `176.2527MHz`. No nonlinear long-term stability claim.
- Field warning: RF gradient pointwise agreement is <0.013% RMS through R=.3r0. DC gradient difference between final two meshes is roughly 5–15% relative in central samples (absolute ≤~1.8e−6 per normalized coordinate per volt); elementwise P2 gradients jump. Hessian mesh convergence and trajectory agreement support this short RF-dominated test, not precision axial phase over many cycles.
- Pseudopotential: safe-sphere minimum at `.20r0` is `0.325µeV` toward −z; radial x cut at `.20r0` is ~10meV. These are local acceptance/cut values, not actual escape barriers.

**Decision:** conditional GO for this deterministic numerical reference case only. **NO-GO for stochastic entropy/RNG simulations** until a bounded multi-axial-period full-RF integration demonstrates phase/energy stability and the P2 gradient/interpolation uncertainty is controlled. Continue no collisions, baths, cooling, lifetime, RNG or remote work unless separately authorized. Report and scripts: `M3A_DETERMINISTIC_ION_DYNAMICS_RU.md`, `scripts/run_m3a_deterministic.py`, `scripts/m3a_finalize_linear_modes.py`, `scripts/m3a_finalize_pseudopotential.py`, `scripts/m3a_mathieu_benchmark.py`, `scripts/m3a_validate_field_region.py`, `scripts/m3a_mesh_trajectory_check.py`.

## M3a.1 — фиксированная DC сила и допуск длительных trajectories (28 Oct 2026)

**Результат:** NO-GO для расширенного интегрирования на текущем FEM поле. Сравнение P2 сеток hDC=.15/.10/.075r0 при общем `VDC=0.08697882845042744 V/ring` даёт локальные axial Hessian frequencies 997.220/999.234/1,000.001 kHz. Но прямые DC forces на конечных осевых/радиальных смещениях расходятся между двумя последними mesh на 3.96–12.66%. На fine сетке `|Fz(+.1r0)|` и `|Fz(−.1r0)|` отличаются примерно на 13.7% в номинально симметричном устройстве.

Полиномиальный скалярный potential fit порядков 2–4 остаётся на 4.02–6.50% RMS от elementwise P2 gradient в `0.30r0`; reduced force пока не является проверенной основой длительных trajectories. Частота 1 MHz была закалибрована по Hessian h=.075, поэтому не считается независимой. M3a.1 не запускал интеграции на 0.1/1/3/5 axial periods; никаких новых данных о confinement/lifetime/escape нет.

**Следующий рекомендуемый шаг после отдельной авторизации:** удержать фиксированное VDC, геометрию и terminal normalization, независимо уточнить central volume минимум на трёх уровнях, проверить реальную сеточную статистику до solve, затем проверить силу/зеркальную симметрию и direct-vs-recovered gradients. Только после прохождения этого gate выполнить staged long-time integration. Оценка прямого P2 benchmark: ~10 min на одну 1 MHz осевую моду и ~50 min на пять периодов одной траектории (80 RK4 steps/RF period); это cost estimate, не запущенные расчёты. Полный отчёт и воспроизведение: `M3A1_EXTENDED_DYNAMICS_VALIDATION_RU.md`, `scripts/m3a1_dc_force_audit.py`, `scripts/m3a1_reduced_dc_field_audit.py`, `scripts/m3a1_field_cost_benchmark.py`, `tests/test_m3a1_audit.py`.

## M3a.2 — центральная сила и зеркальная симметрия (28 Oct 2026)

На неизменной геометрии и одном fixed `VDC=0.08697882845 V/ring` построены mesh-only проверки hcenter=.18/.09/.06r0 при DC-surface h=.075r0 и неизменных RF/end h=.30r0. Фактический median edge length внутри .15r0 упал .384→.202→.146r0; оценка DOF осталась около 214–215k. Использованы базовое ранее решённое поле hcenter=.18r0 и только два новых P2 solve.

Локальные Hessian Hzz между `.09→.06r0` изменяются на 0.014% при fit radius .35r0 и на 0.92% при .50r0, но pointwise forces не сходятся: max relative vector difference 10.55%, max relative magnitude difference 4.78% по выбранным 22 точкам. DC potential-even residual уменьшается, gradient parity RMS остаётся 2.73%; осевой |Fz| mismatch при ±.10r0 снизился с 13.7% на старой сетке до примерно 4.9%. Вероятный источник — unstructured mesh asymmetry и P2 gradient recovery; terminal counts и tags симметричны, но независимого PDE solver proof нет.

**Decision: NO-GO для 0.1/1 axial-period trajectories.** Условная интеграция не запускалась, как требует validation gate. Следующая работа — улучшить central edge resolution/symmetry и доказать force convergence при фиксированном voltage; stochastic/entropy/RNG остаются закрыты. Полный отчёт и artifacts перечислены в `M3A2_CENTRAL_FORCE_VALIDATION_RU.md`.

## M3a.3 — symmetry-reduced FEM and force convergence (28 Oct 2026; current status)

**This section supersedes the M3a.2 recommendation that the next numerical work should start by repairing central mesh sizing.** The M3a.3 audit confirmed exact CAD and drive-pattern reflection symmetry, implemented octant natural-Neumann conditions only on open vacuum plane facets, and assigned conductor cross-sections the proper terminal Dirichlet values. Every external facet on two solved meshes was classified; unclassified count is zero. No weighting-field problem was reduced to an octant.

Three true central mesh levels were generated with unchanged DC/RF electrode-surface sizes and outer geometry. Actual median/P90 tetrahedral edge lengths in `r<0.15r0` are `.03298/.04852`, `.02184/.03002`, and `.01787/.02344r0` for targets `.030/.018/.014`; the last case is mesh-only. P2 DOFs are 33,014 / 45,824 / 63,759. The two solved octants retain the fixed M3a.2 voltage `VDC=.08697882845042744 V/ring` and `VRF=10 V differential peak`.

DC force changes `.030→.018` are max 2.512% over the fixed-point set and 2.088% over six generic off-axis samples. The fine full-domain vs fine octant DC force differs up to 5.243% (3.416% in magnitude). Unconstrained DC Hessian fits shift Hzz by 1.5–1.7% between octants; RF quadrupole Hessian and generic-point RF field are highly stable. Full-domain scalar potential mirror residual is small (`<=8.34e-8 V` per 1-V DC basis at tested reflected pairs), but direct P2 gradient mirror mismatch remains up to 5.09%, consistent with elementwise gradient jumps plus unstructured mesh asymmetry; the causes are not fully separable with the same P2 solver.

**Decision: NO-GO for deterministic trajectory extension.** The 1% force criterion is not met, and several low-quality tets remain in the central volume. The `.014` solve was not attempted because system available RAM was about 1.28 GiB and factor-memory extrapolation could approach/exceed that after matrix/factor workspaces. Harmonic potential fits reproduce held-out direct gradients well, but are derived from the same FEM potential and do not establish PDE convergence. Next method: one bounded adaptive P2/error-estimator study on the octant; if independent force uncertainty remains inadequate, plan either P3 or BEM as a separate benchmark. Do not run those methods in parallel or proceed to ion trajectories, stochastic dynamics, entropy or RNG until separately authorized and the force gate passes. Full details/artifacts: `M3A3_SYMMETRIC_FEM_VALIDATION_RU.md`, `configs/m3a3_symmetry_refinement.yaml`, scripts `m3a3_*`, results `m3a3_octant_*` and `m3a3_full_octant_force_comparison.json`.

## Stage R0 — reproducible repository and physical feasibility audit (29 Oct 2026)

R0 is preparation only; no remote FEM or stochastic campaign was run. The bounded species module (`src/iontrap_rng/species_screen.py`) compares 40Ca+, 88Sr+, 171Yb+, 9Be+ and electron under fixed-q/V scaling and 100MHz/1GHz RF-cap scenarios. It is an ideal Mathieu/pseudopotential estimate, not an extension of the finite 3D trap. `src/iontrap_rng/adaptive.py` provides residual/flux-jump estimator assembly and Dörfler marking helpers; solver-specific residual and flux extraction and an independent goal-error certificate remain unimplemented.

The effects register (`PHYSICAL_EFFECTS_REGISTER.json`) covers FDT/JN noise, thermal motion and cooling, surface noise/charge, gas collisions, image charge, micromotion/technical noise, ion loss, field emission, Si transport, oxide charging, fabrication tolerances, thermal/mechanical/magnetic/quantum effects, readout backaction, loading and total power. Unknown surface/circuit/package properties remain unresolved, not negligible. Particle-free receiver noise is a mandatory baseline before any ion-mediated entropy claim.

Remote scripts/config were prepared but not executed. First conditional job: one sequential solve on existing `.014` octant P2 mesh; config requires >=13 GiB free RAM and >=10 GiB disk, single worker. This is a preflight guard, not proof that peak usage remains below the limit. M3a.3 remains numerical NO-GO for dynamics: `.030→.018` DC force changes up to 2.512%, and fine full-vs-octant difference reaches 5.243%. Do not enter dynamics until adaptive/independent force validation passes.

At the end of R0, public GitHub publication was blocked by an invalid CLI token and repository-creation limitations of the connector. This was resolved during R0.1: the public repository is now available at `https://github.com/snejokzboxik/ion-trap-rng-feasibility`; the initial commit was pushed and a fresh clone plus 37-test run succeeded. The raw `LITERATURE_SESSION_NOTES.json` remains excluded because redistribution rights for embedded excerpts are unassessed. No license file/grant is assigned.

Solver note: saved M3a.3 P2 matrix after grounded Dirichlet elimination is expected SPD, so SciPy CG with an SPD preconditioner is a candidate, not a demonstrated replacement. Smoothed aggregation AMG could reduce direct LU fill, but hierarchy memory and P2-force accuracy need a controlled `.018` cross-check (true residual, force/Hessian/flux, peak RSS, runtime) before any `.014` run; no new solver installed or benchmarked in R0. See SciPy CG and PyAMG official docs linked in the Windows guide. `adaptive.py` provides dimensional residual/jump terms and Dörfler marking, but force-functional/dual-weighted error estimation still needs a validated FEM adapter.

## Stage R0.1 — public GitHub publication and clone validation (29 Oct 2026)

The repository `snejokzboxik/ion-trap-rng-feasibility` is PUBLIC. The original clean local commit `aa73cfb439f3f439bcbf2087ff983c53162c8b28` was pushed to `main` and confirmed equal to `git ls-remote origin refs/heads/main`. `gh repo view` confirmed public visibility. A genuine new public clone was created under the separate `research/_r0_1_public_clone_verify_20261029` directory; it contained the same 196 tracked files as the local repo and passed all 37 tests. Clone worktree was clean. The publisher scan had 0 sensitive-pattern findings, manifest SHA/size inventory matched all 195 non-self files, no tracked file exceeded 50MB (largest 5,711,193 bytes), and no PDFs were included. Raw `LITERATURE_SESSION_NOTES.json` remains local, ignored and excluded (SHA256 `1edd316966462575960e7cd5e6851416e34a9d2b6d0e840f9b5bc362a39c8670`) because third-party snippet redistribution rights are unresolved. No license file/grant was assigned. No FEM/remote/stochastic/RNG calculation was launched.

The manifest and handoff were updated after the initial push to reflect publication. These documentation updates are included in the follow-up commit, which is pushed and verified against `origin/main`; preserve the initial published commit as scientific provenance.
