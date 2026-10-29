# M3a — детерминированная динамика одного иона в 3D RF trap

**Результат:** ограниченные короткие траектории (^{40}\mathrm{Ca}^+) в сохранённых P2 FEM-полях M2b.2 прошли численные проверки времени и интегратора и остались в центральной исследованной области. Локальный Floquet анализ показывает устойчивые радиальные и осевую моды выбранного reference drive. Это **не** подтверждение долгосрочного физического удержания или функционирующего RNG. Следующий stochastic/entropy этап не выполнялся.

## 1. Проверка FEM входа и валидированная область

Для расчётов использованы только сохранённые effective P2 поля M2b.2:

- `results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz`;
- соответствующее `results/m2b2_p2_effective_h0.75_hc0.18_dc0.075.json`.

Идентификатор модели в JSON проверен скриптом: 215,157 P2 DOF; DC поверхность `0.075r0`, центр `0.18r0`, RF поверхности и осевые концы `0.30r0`, дальняя область до `0.75r0`. Это finite-length модель четырёх отдельных гиперболических RF ветвей длиной `10r0`, кольцевых axial DC electrodes при `z=±5.5r0`, grounded outer enclosure и пустой вакуумной области. Электроды идеальные проводники; поверхности CAD полигональные. Это та же геометрия, что M2b.2, не новый идеализированный аналитический заменитель.

Нормировка восстановлена из linear superposition boundary solves в `scripts/run_m2b1_p2_physical_check.py`:

- RF поле `0.5*(right+left−top−bottom)` даёт `+0.5 V` на одной паре и `−0.5 V` на другой, то есть 1 V differential. Поэтому физическая амплитуда 10 V differential реализуется умножением RF FEM potential на 10 V.
- DC поле `dc_plus+dc_minus` даёт 1 V на каждом из двух DC колец при нуле на прочих терминалах. `Vdc` ниже — напряжение **на каждом кольце**.
- FEM coordinates заданы как ξ=r/r0, поэтому для физического поля ∇Φ=(1/r0)∇ξΦ. В JSON потенциалы хранятся как V на один вольт соответствующего terminal pattern.

Точечный evaluator реализован в `src/iontrap_rng/fem_fields.py`. Его P2 potential/gradient проверены на точном квадратичном трёхмерном manufactured solution; тест даёт аналитический потенциал и градиент с точностью порядка roundoff. P2 potential непрерывен, P2 elementwise gradient в общем случае скачет на гранях тетраэдров. Это ключевое численное ограничение траекторной force interpolation.

### Проверка поля в объёме траекторий

`scripts/m3a_validate_field_region.py` сравнил поля на 256 детерминированных Sobol-точках внутри шаров радиуса `0.05`, `0.10`, `0.20`, `0.30r0` между M2b.2 DC-surface mesh levels `0.10r0` и `0.075r0`. Центр, RF и end refinement у этих сеток одинаковы.

| radius | RF ∇Φ RMS difference / fine RMS | DC ∇Φ RMS difference / fine RMS | абсолютная RMS разность DC (V на normalized coordinate, на 1 V terminal pattern) |
|---:|---:|---:|---:|
| 0.05r0 | 0.0124% | 14.62% | 7.65e−7 |
| 0.10r0 | 0.00684% | 7.00% | 7.31e−7 |
| 0.20r0 | 0.00412% | 5.43% | 1.11e−6 |
| 0.30r0 | 0.00358% | 5.79% | 1.79e−6 |

Высокий относительный DC процент вблизи нулевого поля обусловлен малым знаменателем, но в области конечного смещения остаётся заметная относительная погрешность слабой DC силы. Это не отменяет 0.526% изменения регионального Hessian между последними M2b.2 сетками, но показывает, что mesh-only comparison кольцевых поверхностей не является полной pointwise оценкой axial force. Следовательно, здесь надёжнее всего линейная кривизна около центра и RF-dominated radial motion; субпроцентная точность слабой DC силы во всём объёме не доказана. Все trajectories ниже остались внутри `0.16r0`, то есть внутри проверенного объёма `0.30r0`.

## 2. Физические параметры и уравнения

Ион: (^{40}\mathrm{Ca}^+), (m=39.96259098u=6.63594436\times10^{-26}\,\mathrm{kg}), (Q=e); (r_0=1\,\mu\mathrm m). Из finite-field radial Hessian `kRF=Hxx=0.9998641` на вольт differential использовано:

\[
q_{\rm eff}=\frac{2 Q V_{\rm RF} k_{\rm RF}}{m r_0^2\Omega^2},\qquad
\Omega=1.5537256\times10^{10}\;\mathrm{rad/s},\quad
f_{\rm RF}=\frac{\Omega}{2\pi}=2.4728311\;\mathrm{GHz},
\]

для (V_{\rm RF}=10\,\mathrm V) differential peak и (q_{\rm eff}=0.2). DC axial curvature `Hzz=1.879916e−4 V/r0²` на 1 V **на каждом** кольце даёт (1.879916\times10^8\,\mathrm{V/m^2}) при (r_0=1\,\mu\mathrm m). Для локальной целевой (f_z=1\,\mathrm{MHz}):

\[
V_{\rm DC}=\frac{m(2\pi f_z)^2 r_0^2}{Q H_{zz,\rm FEM}}=0.086979\,\mathrm V
\]

на каждом кольце. Это калибровка по local Hessian, не экспериментально установленное напряжение и не гарантированная частота полного nonlinear trap.

Интегрировалось исходное time-dependent force equation по FEM potentials:

\[
m\ddot{\mathbf r}=-Q\nabla\left[V_{\rm RF}\phi_{\rm RF}(\mathbf r/r_0)\cos(\Omega t+\varphi)
+V_{\rm DC}\phi_{\rm DC}(\mathbf r/r_0)\right].
\]

С (\tau=\Omega t/2), (\boldsymbol\xi=\mathbf r/r_0):

\[
\boldsymbol\xi''=-\gamma\left[V_{\rm RF}\nabla_\xi\phi_{\rm RF}\cos(2\tau+\varphi)
+V_{\rm DC}\nabla_\xi\phi_{\rm DC}\right],\qquad
\gamma=\frac{4Q}{m r_0^2\Omega^2}=0.0400054\;\mathrm V^{-1}.
\]

Интегратор использовал P2 potential gradient, не pseudopotential. RF work included:

\[
\epsilon=\tfrac12|\boldsymbol\xi'|^2+
\gamma(V_{\rm RF}\phi_{\rm RF}\cos(2\tau+\varphi)+V_{\rm DC}\phi_{\rm DC}),\qquad
\frac{d\epsilon}{d\tau}=-2\gamma V_{\rm RF}\phi_{\rm RF}\sin(2\tau+\varphi).
\]

Тем самым RF drive может обмениваться энергией с частицей; сохранение mechanical energy не ожидалось.

## 3. Локальная устойчивость и Floquet режимы

Regional Hessians из M2b.2 были симметризованы перед линеаризацией; относительная antisymmetric часть до симметризации: RF `2.71e−6`, DC `2.47e−3`. Для уравнения вариаций использовалась матрица

\[
\mathbf\xi''+\gamma[V_{\rm RF}H_{\rm RF}\cos(2\tau)+V_{\rm DC}H_{\rm DC}]\mathbf\xi=0,
\]

а 6×6 monodromy интегрировалась на периоде `π` в τ DOP853. Результат:

| mode | Floquet β | secular frequency | multiplier modulus |
|---|---:|---:|---:|
| axial | 0.00080879 | 1.0000005 MHz | 1 within 1e−13 |
| radial 1 | 0.14254986 | 176.250856 MHz | 1 within 1e−13 |
| radial 2 | 0.14255016 | 176.251230 MHz | 1 within 1e−13 |

Finite FEM idealized radial values agree with ideal (a=0, |q|=0.2) Mathieu Floquet (176.252675\,\mathrm{MHz}) within `0.0011%`; small-q approximation (q f_{RF}/(2\sqrt2)=174.856\,\mathrm{MHz}) is lower by about `0.8%`. This is a local linear-stability result. Multipliers are symplectic/unit-modulus within numerical precision, but it does not prove finite-amplitude or long-term confinement.

The RF null estimated from the transverse Hessian lies about `(1.82, −0.61) pm` from the geometric center; axial RF field has near-zero axial curvature so a 3D RF-null inversion in z is ill-conditioned and was not used. Local pseudopotential equilibrium is approximately `(1.84, −0.63) pm, z=−1.79 nm`; the predicted peak excess micromotion from its radial displacement is only `0.194 pm`. These estimates inherit FEM gradient recovery and geometry asymmetry and must not be treated as device alignment tolerances.

## 4. Интеграторы, шаг, RF фаза и trajectories

Идеальный Mathieu reference (x''+2q\cos(2\tau)x=0) проверен отдельной монодромией. RK4 сравнен с DOP853 за один RF period:

| RK4 steps/RF period | relative monodromy Frobenius error | β error | det(M)−1 |
|---:|---:|---:|---:|
| 20 | 2.60e−6 | 4.44e−6 | −9.91e−7 |
| 40 | 1.60e−7 | 2.81e−7 | −3.12e−8 |
| 80 | 9.91e−9 | 1.76e−8 | −9.77e−10 |
| 160 | 6.18e−10 | 1.10e−9 | −3.06e−11 |

В конечной FEM модели интегрированы следующие начальные положения с нулевой скоростью, фазы RF `0` и `π/2`, время `20 RF periods = 8.087 ns`, либо выход при `|r|>0.30r0`. RK4 выполнен при 40/80/160 шагов за период для near-center phase 0 и при 80 шагах для остальных условий; DOP853 (`rtol=2e−10`, `atol=2e−12`, `max_step=T_RF/80`) — для phase 0.

| начальное ξ | RF phase | max |r|/r0 (RK4) | max |r|/r0 (DOP853) | исход |
|---|---:|---:|---:|---|
| (0.02,0,0) near-center | 0 | 0.0200 | 0.0200 | остаётся внутри |
| (0.02,0,0) near-center | π/2 | 0.0383 | — | остаётся внутри |
| (0.08,0,0) radial | 0 | 0.0800 | 0.0800 | остаётся внутри |
| (0.08,0,0) radial | π/2 | 0.1532 | — | остаётся внутри |
| (0,0,0.10) axial | 0 | 0.1000 | 0.1000 | остаётся внутри |
| (0,0,0.10) axial | π/2 | 0.1000 | — | остаётся внутри |
| (0.05,0.04,0.08) combined | 0 | 0.1025 | 0.1025 | остаётся внутри |
| (0.05,0.04,0.08) combined | π/2 | 0.1281 | — | остаётся внутри |

Окно 20 RF периодов содержит около 1.43 radial secular period, но только `0.0081` axial period при 1 MHz. Поэтому короткая axial trajectory не измеряет осевую частоту; её частота получена из local Floquet analysis. Ни один расчёт не приближался к FEM escape boundary, и тест не оценивает ионную lifetime.

В FEM trajectory сравнения RK4 80 steps/period с DOP853: максимальная разность конечной позиции среди phase-0 запусков `1.54e−7 r0`; при шаговом тесте near-center case ошибки относительно DOP853 составили `1.54e−7`, `7.77e−9`, `9.56e−9 r0` при 40/80/160 шагах. Остаток energy/work balance не превысил `3.1e−10` в dimensionless energy для RK4 и `9.7e−13` для DOP853. Плато около `1e−8r0` и немонотонность последних двух шагов согласуются с влиянием кусочно-непрерывного FEM градиента; не объявляется формальный четвёртый порядок этой конечной сеточной force.

При одинаковых условиях между mesh `hDC=0.10` и `0.075r0` конечная позиция отличается на `≈7.1e−6r0` near-center и `≈1.4e−5r0` combined. Это очень мало для коротких RF-dominated траекторий, но не проверяет накопленную фазовую ошибку на многих axial periods. Максимальные начальные RF-period внутрипериодные RMS-колебания для `0.08r0` radial case порядка `7.9e−3r0` (≈7.9 nm); эта оценка включает secular drift и не является чистым спектральным выделением intrinsic micromotion. leading-order (qA/2=8\,\mathrm{nm}) даёт ожидаемый peak-scale.

## 5. Псевдопотенциал и выход

По конечному RF полю использована формула

\[
U_{\rm eff}(\mathbf r)=\frac{Q^2|\mathbf E_{\rm RF,peak}(\mathbf r)|^2}{4m\Omega^2}+Q\Phi_{\rm DC}(\mathbf r).
\]

В шаре радиуса `0.20r0`, просканированном по 31 полярному ×64 азимутальным направлениям, минимальная разность `Ueff−Ueff(center)` равна `3.2494e−7 eV = 0.325 µeV` вдоль оси `−z`. Radial x-cut при `0.20r0` даёт около `10 meV`. Это **safe-radius local acceptance energy**, зависящая от выбранной искусственной сферической границы. Она не является фактической trap depth. Потенциальные седла у RF-электродов, осевых торцов и grounded enclosure не оценивались: поля и геометрия возле поверхностей/endpoints недостаточно сходились для escape-path вывода.

Рисунок траекторий и local potential cuts: `results/figures/m3a_deterministic.png`.

## 6. Неопределённости и научные ограничения

**Пройдено:** identity и normalization saved FEM fields; точный P2 gradient benchmark; короткие детерминированные trajectory comparisons двух интеграторов; timestep convergence в идеальном Mathieu limit; work accounting; local finite-field Floquet unit-modulus modes; confinement внутри radius `0.16r0` на коротком RF окне.

**Условно/не пройдено:** DC pointwise gradients между сетками отличаются на 5–15% относительно слабой DC силы в проверенной сфере; RF градиенты согласуются значительно лучше. Геометрический центр и endpoint/rf mesh не уточнялись в этой трёхсеточной последовательности; оценки local Hessian зависят от fit window. P2 gradient discontinuous, поэтому numerical trajectory has mesh-imprinted force changes. Накопление ошибки за много RF/axial periods, RF heating/heating-free reference, collision lifetime и nonlinear escape не проверялись. Local pseudopotential не заменяет full RF escape analysis.

**Не моделировалось:** collisions, stochastic forces/noise, cooling, loading/restart, material conductivity, dielectric/support/package stack, RF circuit, readout, entropy, RNG throughput. Детерминированная trajectory не является источником random bits.

**Рекомендация:** M3a deterministic short-run — conditional GO для локальной идеализированной модели. Следующий этап должен быть численно ограниченной проверкой накопленной фазовой/энергетической ошибки хотя бы на нескольких axial periods и улучшением/сглаживанием P2 force recovery с center-region refinement; только затем стоит вводить физически заданное охлаждение/шумы. Нынешний результат — **NO-GO для entropy или RNG simulations**, поскольку axial-cycle lifetime и стохастическая модель здесь не валидированы. Запуск stochastic этапа не выполнялся.

## 7. Воспроизводимость

- Параметры: `configs/m3a_ca40_reference.yaml`.
- P2 evaluator: `src/iontrap_rng/fem_fields.py`; exact quadratic test: `tests/test_fem_fields.py`.
- Пространственная field audit: `scripts/m3a_validate_field_region.py`, `results/m3a_field_region_audit.json`.
- Mathieu reference: `scripts/m3a_mathieu_benchmark.py`, `results/m3a_mathieu_benchmark.json`.
- Trajectories/work: `scripts/run_m3a_deterministic.py`, `results/m3a_deterministic_results.json`, `results/m3a_trajectory_traces.npz`.
- Mesh sensitivity: `scripts/m3a_mesh_trajectory_check.py`, `results/m3a_mesh_trajectory_check.json`.
- Linear modes and pseudopotential postprocessing: `scripts/m3a_finalize_linear_modes.py`, `scripts/m3a_finalize_pseudopotential.py`.
- Figure: `scripts/plot_m3a_results.py`, `results/figures/m3a_deterministic.png`.
- Tests: `..\\..\\.venv\\Scripts\\python.exe -m pytest -q --basetemp results/pytest_m3a` — **29 passed**.
