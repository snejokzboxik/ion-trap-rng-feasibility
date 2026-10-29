# M1 и предварительный M1.5: аналитическая проверка

**Статус:** пакет M1 и небольшой M1.5 sensitivity screen выполнены локально. FEM, полные particle trajectories, seeded entropy simulations, крупные sweeps и remote вычисления не запускались. Experimental TRNG claims не делаются.

## 1. Реализовано

Пакет src/iontrap_rng содержит SI-формулы: четыре отдельные ветви \(x^2-y^2=\pm r_0^2\), branch parameterization/normals/residual; ideal quadrupole potential and fields с \(V_{\rm RF}\) как differential peak (пары \(\pm V/2\)); Mathieu \(a_i,q_i\), axial DC quadrupole and radial defocusing; Floquet fundamental-matrix integration; small-q secular frequencies; peak-field pseudopotential and safe-radius acceptance energy; fixed-q scaling exponents; image-current, Johnson, FDT heating, gas collision, survival and joule/useful-bit screens; nominal \(^{40}\mathrm{Ca}^+\), \(^{88}\mathrm{Sr}^+\), \(^{174}\mathrm{Yb}^+\). Isotope masses use mass-number approximations.

Reproduction commands:

~~~powershell
.\.venv\Scripts\python.exe research\ion_trap_rng\scripts\run_m1_reference.py
.\.venv\Scripts\python.exe research\ion_trap_rng\scripts\run_m15_screen.py
~~~

The test suite reports **12 passed**. Floquet integrator comparison uses DOP853 and RK45; it does not validate full stochastic ion dynamics.

## 2. M1 validation results

| Test | Result |
|---|---|
| Four branches satisfy their hyperbola equation | Maximum absolute residual \(<2\times10^{-25}\,\mathrm{m^2}\) in test at \(r_0=3.2\) µm |
| Outward normals | Unit length and point from central vacuum region toward each branch |
| Ideal potential/field | Finite-difference Laplace check; at \(x=r_0\), RF peak face potential is \(+V_{\rm diff}/2=5\) V |
| Mathieu coefficients | Reconstructed \(|q_x|=0.2\), \(q_y=-q_x\) |
| Floquet, \(a=0,q=0.2\) | Stable; \(\beta=0.1425513265\); \(\det M=0.99999999999985\); multipliers lie on unit circle within numerical tolerance |
| Small-q reference | \(q/\sqrt2=0.1414214\); Floquet \(\beta\) differs by ~0.80%, consistent with a low-q asymptotic |
| Stability references | \(q=.91,a=0\) unstable; \(q=0,a=0\) is marginal free motion, not confinement |
| Independent integration methods | DOP853 and RK45 discriminants agree within \(2\times10^{-9}\) test tolerance |
| Axial DC | Finite-difference gradient agrees with \(-\nabla\Phi\); radial field points outward; Laplace Hessian trace is zero; radial frequency falls when \(f_z\) is added |
| Pseudopotential | Direct peak-field energy agrees with algebraic safe-radius acceptance energy |

### Исправленная математическая ошибка

В первой версии axial_dc_field_v_per_m радиальные знаки были ошибочны. Для
\[
\Phi_{\rm ax}=\frac{m\omega_z^2}{2Q}\left[z^2-\frac{x^2+y^2}{2}\right]
\]
получается \(E_x=+m\omega_z^2x/(2Q)\), \(E_y=+m\omega_z^2y/(2Q)\), \(E_z=-m\omega_z^2z/Q\). Радиальные составляющие направлены наружу и дают defocusing. Это исправлено и теперь покрыто gradient/Laplace tests.

## 3. Selected fixed-q ion comparison

Inputs: \(r_0=1\) µm, \(V_{\rm RF}=10\) V differential peak, \(q=.2\), \(U=0\), \(f_z=1\) MHz, \(T=300\) K, \(\eta=.5\). Radial values use small-q pseudopotential with axial defocusing, not exact Floquet exponents with \(a_z\ne0\).

| Ion | \(f_{\rm RF}\) | \(f_{r,x}\approx f_{r,y}\) | \(D_{\eta=.5}\) | \(x_{\rm rms}\), equilibrium-only |
|---|---:|---:|---:|---:|
| \(^{40}\mathrm{Ca}^+\) | 2.472 GHz | 174.78 MHz | 0.0625 eV | 227 nm |
| \(^{88}\mathrm{Sr}^+\) | 1.667 GHz | 117.84 MHz | 0.0625 eV | 227 nm |
| \(^{174}\mathrm{Yb}^+\) | 1.185 GHz | 83.80 MHz | 0.0625 eV | 227 nm |

Equal \(x_{\rm rms}\) follows from fixed-q/fixed-voltage scaling \(\omega_r^2\propto1/m\), not from a demonstrated shared ion temperature. \(D_\eta\) is a selected-radius RF pseudopotential, not finite-trap escape depth. Baseline \(D_\eta/k_BT=2.418\); no lifetime follows without a distribution, cooling model, finite saddle and collision/noise physics.

For Ca+, fixed q=.2 and 10 V, r0=10 µm, 1 µm, 100 nm, 10 nm give RF frequencies 247 MHz, 2.472 GHz, 24.72 GHz, 247.2 GHz. At fixed q,V each decade of miniaturization multiplies RF frequency by 10, holds this ideal acceptance energy constant, and multiplies \(V/r_0\) by 10. This is a similarity law, not proof of fabricated operation.

## 4. M1.5: observability, entropy, cooling, lifetime, energy

### Image-current screen

For Ca+, \(f_r\simeq174.8\) MHz, assumed secular amplitude \(A=100\) nm and effective weighting distance \(d_{\rm eff}=1\) µm:
\[
I_{\rm pk}\sim |Q|\omega_r A/d_{\rm eff}\simeq17.6\ {\rm pA}.
\]
For assumed \(R_{in}=1\) MΩ, \(T=300\) K and ENBW \(B=1\) kHz, Johnson/Norton noise is \(\sqrt{4k_BT B/R}=4.07\) pA rms; simple peak/rms ratio is ~4.32. These are scalar sensitivities, not receiver SNR measurements.

Assumed ranges \(A=10\)–500 nm, \(d_{\rm eff}=0.5\)–5 µm, \(R=10\) kΩ–1 MΩ, ENBW 1 kHz–1 MHz produce ideal current scales 0.352–176 pA, resistor-noise scales 4.07 pA–1.29 nA and extreme ratios 0.00027–43. These corners are independent and not a design envelope. Ramo weighting field, paired-electrode cancellation, RF feedthrough, actual high-frequency impedance, amplifier/noise matching, resonant circuit, ADC/quantization and backaction are absent. ENBW assumes a narrowband receiver centered near 175 MHz; no such receiver is designed yet. Readout feasibility remains unresolved.

### Три источника энтропии

| Source | Analytical screen | Interpretation |
|---|---|---|
| Johnson–Nyquist | 4.07 pA rms in the illustrative 1 MΩ/300 K/1 kHz circuit. Ion response depends on impedance-to-field transfer and Ramo coupling. | Most bounded model/control. Circuit entropy is directly measurable without an ion; trapping is not inherently an advantage. FDT requires noise and dissipation together. |
| Ion-buffer-gas collisions | Assumed $k=0.3$–$3\times10^{-15}$ m³/s, $p=10^{-8}$ Pa, 300 K gives $\Gamma=7.24\times10^{-4}$–$7.24\times10^{-3}$ s$^{-1}$. Independent sparse Poisson event/no-event min-entropy ceiling is $\Gamma/\ln2=0.00104$–0.01045 bit/s. | Ion-mediated stochastic candidate, but low rate at this pressure. A 100 kbit/s binary event-indicator ceiling requires $\Gamma\approx6.93\times10^4$ s$^{-1}$, $p\approx0.096$–0.957 Pa and ENBW >=50 kHz. At that minimum bandwidth, nominal 1 MΩ Johnson noise is 28.8 pA rms, yielding only 0.61 for the assumed 17.6 pA signal/noise scale before omitted effects. For 99% hourly survival, independent loss probability must be <4.03×10⁻¹¹ per collision. Collision/RF energy exchange, cooling and charge exchange remain unmodeled; this is not a no-go theorem. |
| Silicon surface charge / electric-field fluctuations | External literature summary spans roughly \(S_E=10^{-15}\)–\(10^{-6}\,\mathrm{V^2m^{-2}Hz^{-1}}\), mainly 0.1–20 MHz and 30–3500 µm. Substitution at Ca+ ~175 MHz gives \(\dot{\bar n}=8.34\times10^{-4}\)–\(8.34\times10^5\) s\(^{-1}\) and motional input power \(9.66\times10^{-29}\)–\(9.66\times10^{-20}\) W. | Deliberate extrapolation, not a local micron Si/oxide prior. PSD amplitude alone does not provide correlation time or conditional min-entropy rate. |

For a binary event/no-event output sampled at ENBW B=1 kHz, no more than 2B=2 k symbols/s are available regardless of collision rate. Multi-bit analog samples may exceed one bit/sample but require quantizer, distribution and conditional side-information analysis. For that surface-noise sensitivity, no-cooling \(D_\eta/P_{\rm ion}\) is \(1.04\times10^8\)–0.104 s. \(D_\eta\) is only a safe-radius acceptance proxy, not an escape barrier. A stationary motion distribution requires cooling/removal to balance absorbed motional power plus collision/RF heating. Cooling may suppress the sampled fluctuations and costs system power. No primary source is selected: Johnson is the first control; gas and surface noise remain conditional candidates.

The assumed gas coefficient is a sensitivity input, not a measured cross-section for Ca/Sr/Yb with a named neutral. Device-specific \(S_E(f,d,T)\), surface correlation time/length, differential collision cross-section, charge exchange and loss probability remain unbounded.

### First-order system power

Total module power remains unknown: RF source/resonator, ion source/ionizer, cavity/getter, cooling, high-frequency readout/ADC, controller and health checks all count. An assumed 1 mJ startup once per hour adds 0.278 µW average.

| Useful rate | 0.1 mW | 1 mW | 10 mW |
|---:|---:|---:|---:|
| 10 kbit/s | 10 nJ/bit | 100 nJ/bit | 1 µJ/bit |
| 100 kbit/s | 1 nJ/bit | 10 nJ/bit | 100 nJ/bit |
| 1 Mbit/s | 0.1 nJ/bit | 1 nJ/bit | 10 nJ/bit |

These are comparator arithmetic, not predicted device energy. The reported CMOS latch point is 0.116 pJ/bit under its published boundary; a common-boundary audit is still needed. 100 kbit/s and 1 mW are competitiveness (C) scenarios, not physical (P) thresholds.

## 5. Boundary classes and conclusions

- **F, ideal model:** Mathieu/Floquet stability for fixed ideal coefficients; it supplies no direct geometric minimum after nondimensionalization.
- **R:** \(D_\eta\ge\chi k_BT\) depends on temperature, survival/workload and selected \(\chi\); it is not a lifetime law.
- **T + U:** emission/breakdown and RF delivery depend on silicon doping, oxide/adsorbates, edge radius, nanogap process and RF impedance. No validated field limit for this geometry is available.
- **U:** local surface PSD/correlation, gas-specific collision/loss, finite-hyperbola image Green function and receiver coupling are unbounded.
- **E-level T/R/U:** loading, sealed cavity, restart and module budget depend on process and use case; no full-module precedent defines a minimum.

No universal minimum size is reported. A thermal/field inequality provides only a scenario boundary after field allowance and lifetime proxy are stated. The results neither prove nanometer confinement impossible nor support a practical nanometer RNG.

## 6. Is M2 (2D FEM) justified?

M1 validation passes. M1.5 does **not** pass physical entropy-rate or readout feasibility, but it does not exclude all mechanisms: Johnson remains a plausible control, gas collisions trade entropy rate against pressure/loss/cooling, and the surface branch is dominated by local-data uncertainty.

M2 is not justified as a broad feasibility campaign: the ideal hyperbolic field is analytic, and field convergence cannot resolve surface PSD, collision lifetime, cooling, entropy, electronics, loading or total power. One small FEM benchmark could be justified narrowly as solver/geometry verification after selecting a finite truncated-hyperbola design and a weighting-field/readout question. It should use exact analytic exterior boundary data and explicit mesh/domain convergence and must not be called a P/E/C pass. No FEM starts here; separate approval is required.

## 7. Literature anchors used

- Brownnutt et al., *Ion-trap measurements of electric-field noise near surfaces*, DOI [10.1103/RevModPhys.87.1419](https://doi.org/10.1103/RevModPhys.87.1419): measured noise spans, multiple physical mechanisms and geometry dependence; the reported range is not transferred as a device prior.
- Turchette et al., trapped-ion electrical-noise heating, DOI [10.1103/PhysRevA.61.063418](https://doi.org/10.1103/PhysRevA.61.063418): weak-noise heating relation under the stated one-sided PSD convention.
- Viehland & Goeringer, RFQ collision kinetics, DOI [10.1063/1.1691405](https://doi.org/10.1063/1.1691405), and DeVoe, DOI [10.1103/PhysRevLett.102.063001](https://doi.org/10.1103/PhysRevLett.102.063001): buffer-gas damping/heating and nonthermal energy distributions.

## 8. Recommendation

M2 is conditionally useful only for a bounded numerical verification or a specified weighting-field geometry after one concrete finite electrode design exists. If that design question is not yet specified, first constrain receiver weighting/transfer model and device-specific noise/collision inputs. M1 and this preliminary M1.5 do not authorize further FEM; wait for separate approval.

Figure: [M1 analytical scaling and surface-noise sensitivity](results/figures/m1_analytical_screen.png). Numerical outputs: results/m1_reference.json and results/m15_screen.json.
