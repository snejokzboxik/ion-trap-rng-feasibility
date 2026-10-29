# Пределы миниатюризации кремниевой гиперболической RF ion trap как источника энтропии

**Research roadmap v1.4 · 28 октября 2026 · статус: M1–M2b.2 и короткий M3a завершены; M3a.1 force-validation gate завершился NO-GO для длительной динамики; stochastic/entropy modelling не запускалось**

Документ сохраняет исходную программу M1–M8 и исторические результаты предыдущих этапов. Короткий M3a остаётся воспроизводимым историческим reference-case, но свежий M3a.1 fixed-voltage audit ставит numerical NO-GO на длительные траектории из-за DC force convergence. Полного длительного удержания, стохастической энтропии, RNG, реального silicon stack или инженерной интеграции не установлено. Более ранние этапные статусы сохранены ниже как история; актуальное решение задаёт последний addendum.

## Исторический addendum M3a — 28 октября 2026

## Актуальный addendum M3a.1 — численная проверка DC силы

M3a.1 проведён локально на сохранённых трех эффективных P2 FEM сетках. При одном и том же `VDC=0.08697882845042744 V` на каждом кольце локальные осевые частоты из Hessian равны 997.220, 999.234 и 1,000.001 kHz для DC-surface targets 0.15, 0.10 и 0.075r0. Это значение напряжения ранее было откалибровано по последнему Hessian, поэтому близость частот не является независимым подтверждением модели.

Проверка конечных смещений выявила 4–13% различие прямой P2 DC restoring force между двумя последними сетками, до 12.66% при z=−0.10r0; fine mesh также демонстрирует приблизительно 13.7% асимметрию |Fz| при z=±0.10r0 в номинально симметричной геометрии. Потенциальные polynomial reductions степеней 2–4 сохраняют градиентную структуру, но дают 4.02–6.50% RMS расхождение градиента с элементным P2 полем внутри 0.30r0. Для длительной динамики они не валидированы.

**Текущий gate: NO-GO для длительных trajectories и stochastic/RNG этапов.** Это numerical blocker, не физическое опровержение ловушки. Не запускались 0.1/1/3/5 axial-period траектории; предыдущие M3a короткие результаты сохранены и относятся только к окну 8.087 ns. Первый рекомендуемый следующий расчёт — три уровня центральной h-refinement при неизменных геометрии, фиксированном `VDC`, RF, DC surface/end resolution, затем независимая force/symmetry верификация. См. [M3a.1 report](M3A1_EXTENDED_DYNAMICS_VALIDATION_RU.md).

### Краткое историческое резюме M3a

Проведён ограниченный deterministic (^{40}\mathrm{Ca}^+) расчёт в сохранённой P2 FEM геометрии M2b.2 (r0=1 µm; 10 V differential RF peak; q=0.2; обе DC обкладки по ≈0.08698 V для локальной axial частоты 1 MHz). FEM градиент проверен exact P2 quadratic case. Траектории на 20 RF периодах, selected radial/axial/combined initial conditions и RF phases, остались в центральном радиусе ≤0.154r0. RK4/DOP853, шаговая сходимость, энергообмен с RF и finite-field Floquet проверены. Получены локальные mode frequencies около 176.251 MHz радиально и 1.0000005 MHz осево.

Этот результат условен: P2 element-gradient разрывен; между последними сетками слабая DC сила имеет относительное RMS различие порядка 5–15% в центральной сфере, хотя absolute error мала и Hessian меняется на 0.526%. Наблюдение 8.1 ns не покрывает даже 1% axial period. Safe-radius pseudopotential difference 0.325 µeV при 0.2r0 — не trap depth. Surface/end escape, heating over long duration, collisions, cooling and entropy не проверялись. **Conditional GO только для локальной deterministic reference calculation; NO-GO для stochastic entropy/RNG следующего этапа до axial-cycle validation и улучшения/оценки force recovery.** Подробности и воспроизводимые файлы: `M3A_DETERMINISTIC_ION_DYNAMICS_RU.md`.

## 1. Исследовательский вопрос и границы вывода

Главный вопрос: существует ли непустое множество параметров **полного устройства**, допускающее удержание одного иона, наблюдаемую физическую энтропию, повторную загрузку, приемлемую энергию и реалистичное изготовление?

Предварительная оценка: микрометровая ловушка физически правдоподобна как исследовательская система; работоспособный нанометровый RNG при комнатной температуре пока не обоснован. Опасно сочетание ограниченной глубины, RF frequency, поверхностей, слабого single-ion readout и размеров вспомогательных узлов. Первый ограничивающий механизм нужно установить, а не назначить заранее.

Вычислительная работа может дать воспроизводимые **условные границы** и отрицательные результаты для заданных моделей. Она не подтверждает экспериментально спектр шума новой поверхности, герметичность технологии, долговечность эмиттера или физическую случайность образца. Положительный результат — model-based feasibility, а не демонстрация TRNG или сертификация.

| Статус | Значение | Пример |
|---|---|---|
| Known physics | Проверенная физика в области применимости | Laplace equation, Mathieu stability |
| Approximation | Упрощение, требующее контроля | pseudopotential при малом q |
| Hypothesis | Проверяемое предположение | полезная измеряемая энтропия движения одного иона |
| Unverified assumption | Неподтверждённый входной параметр | patch noise при 100 nm; герметичная полость 1 mm³ |

### Уровни успеха и классы границ

Результаты классифицируются по трём независимым уровням: **P — physical feasibility** (устойчивая локальная конфигурация иона и наблюдаемая физическая стохастическая компонента с обоснованной условной entropy model); **E — engineering feasibility** (автономные загрузка/перезапуск, упаковка, изготовление, тепло и питание при явно заданном workload); **C — practical competitiveness** (сопоставление полной энергии на полезный бит, rate, footprint, reliability и security-модели с CMOS TRNG). P не подразумевает E; E не подразумевает C. Все числовые цели из этого плана — сценарии проекта, не универсальные физические границы.

Каждая заявленная граница миниатюризации получает один или несколько ярлыков: **F** — фундаментальная при записанных предположениях; **T** — технологически зависимая; **R** — зависит от целевых требований/workload; **U** — model-dependent или пока не ограниченная данными. Указывать следует значение или диапазон, применимость модели, источник неопределённости и чувствительность решения. Универсальный минимальный размер из одной расчётной точки не выводится.

В рабочей папке уже есть исследования semiconductor impedance и прямоугольных rails; README использует продольную ось x. Новая задача требует hyperbolic electrodes и ось z. Старые fields и mesh convergence нельзя переносить автоматически. Потенциально переиспользуемы только независимо проверенные material, power-balance и RF-transport routines. При проверке папка не определилась как Git repository; подготовка Git остаётся отдельным организационным действием.

## 2. Четыре самостоятельные ветви гиперболы

### 2.1. Идеальная опорная геометрия

\(r_0>0\) — расстояние от оси до вершины RF electrode. В XY-plane задаём:

\[
\begin{aligned}
\Gamma_{x+}&:(x,y)=r_0(\cosh u,\sinh u),\\
\Gamma_{x-}&:(x,y)=r_0(-\cosh u,\sinh u),\\
\Gamma_{y+}&:(x,y)=r_0(\sinh u,\cosh u),\\
\Gamma_{y-}&:(x,y)=r_0(\sinh u,-\cosh u),\quad u\in\mathbb R.
\end{aligned}
\]

Первые две удовлетворяют \(x^2-y^2=r_0^2\), вторые — \(x^2-y^2=-r_0^2\). Центральная связная вакуумная область:

\[
\mathcal D_\infty=\{(x,y):|x^2-y^2|<r_0^2\}.
\]

Область имеет диагональные каналы; это не круглая апертура. Круг \(\rho<r_0\) целиком находится внутри. Размер \(2r_0\) — межвершинный зазор.

Противоположные электроды RF-paired:

\[
\Phi|_{\Gamma_{x\pm}}=+\frac{U+V_{\rm RF}\cos\Omega t}{2},
\quad
\Phi|_{\Gamma_{y\pm}}=-\frac{U+V_{\rm RF}\cos\Omega t}{2}.
\]

**Voltage convention:** \(V_{\rm RF}\) — peak differential voltage между парами; peak amplitude одного электрода относительно common mode равна \(V_{\rm RF}/2\). Это устраняет ошибку в два раза в Mathieu q и в четыре раза в pseudopotential. Опорный случай пользователя: \(U=0\).

В свободном объёме \(\nabla^2\Phi=0\), точное решение:

\[
\Phi_{\rm id}=\frac{U+V_{\rm RF}\cos\Omega t}{2r_0^2}(x^2-y^2).
\]

На искусственной внешней границе finite numerical domain сначала задаётся **это решение**. Заземление внешнего круга даёт другую задачу. Для неограниченного домена дополнительно ограничиваем допустимое поведение на бесконечности; одних значений на четырёх ветвях недостаточно для произвольных растущих harmonic functions.

### 2.2. Finite silicon electrodes

Ограничиваем \(|u|\le u_{\max}\) и экструдируем вдоль \(-L/2\le z\le L/2\). Четыре лицевые поверхности сохраняют гиперболическую форму. Материал находится за ними, вне центральной полости.

Для тонких тел удобна normal-offset parameterization:

\[
\mathbf R_j(u,s,z)=(\mathbf r_j(u)+s\mathbf n_j(u),z),
\quad 0\le s\le t_j(u,z),
\]

где \(\mathbf n_j\) направлена из вакуума в silicon. Например,

\[
\mathbf n_{x+}=\frac{(\cosh u,-\sinh u)}{\sqrt{\cosh2u}}.
\]

Остальные тела получают отражениями/поворотами. Обязательны проверки self-intersections, connectedness и electrical isolation. Радиус кривизны \(R_c=r_0(\cosh2u)^{3/2}\). Для простой normal-offset оболочки \(t<r_0\) — удобное достаточное локальное ограничение; более толстый backing строится отдельно.

Design variables: r0, u_max, L, thickness_profile, backing_width, edge_radius, anchor_geometry, feed_geometry, contact_impedance, dc_gap, oxide_thickness, surface_roughness, dopant_density.

Толщина и backing — независимые параметры оптимизации. Предпочтительны тонкие тела с локальными anchors, а не массивные wedges. Усечения скругляются с заданным edge_radius: максимум поля на математически остром крае mesh-dependent и непригоден для оценки breakdown.

Oxide и dielectric supports окружают **пустую полость**, необходимую для движения иона. Ион в сплошном SiO₂ не является свободным trapped ion. Native oxide, interface charge и dielectric losses включаются явно.

Cylindrical rods, rectangular rails и ring-and-endcap traps допустимы только как литературные/валидационные контроли, не как основная архитектура.

### 2.3. Метрики quadrupole accuracy

В центральном 2D сечении единичный RF basis potential:

\[
\phi_{\rm RF}=c_0+\sum_{n\ge1}(\rho/r_0)^n
[c_n\cos n\theta+s_n\sin n\theta].
\]

При differential voltage 1 V идеал: \(c_2=1/2\), остальные коэффициенты нулевые. При точной quadrupole symmetry допустимы \(n=2,6,10,\ldots\). Нежелательные \(n=1,3,4\) позволяют диагностировать нарушения симметрии. В 3D нужны harmonic-polynomial expansion, Hessian и z-dependent coefficients: 2D expansion недостаточно.

Для рабочей области \(\mathcal R:\rho\le\eta r_0\), например \(\eta=0.3,0.5\):

\[
\epsilon_E=
\left[\frac{\int_{\mathcal R}|\mathbf E_{\rm num}-\mathbf E_{\rm ref}|^2dV}
{\int_{\mathcal R}|\mathbf E_{\rm ref}|^2dV}\right]^{1/2}.
\]

Pointwise relative error в RF null не определена. Там оцениваем absolute field residual, null displacement и Hessian. Ошибку quadrupole amplitude отделяем от nonquadrupole residual.

## 3. Governing equations и axial confinement

### 3.1. Mathieu model

Для заряда \(Q=Ze\), массы m, без collisions и surface forces:

\[
\ddot x+\frac{Q(U+V_{\rm RF}\cos\Omega t)}{mr_0^2}x=0,\quad
\ddot y-\frac{Q(U+V_{\rm RF}\cos\Omega t)}{mr_0^2}y=0,\quad
\ddot z=0.
\]

При \(\tau=\Omega t/2\) используем convention

\[
u_i''+[a_i-2q_i\cos2\tau]u_i=0,
\]

\[
a_x=\frac{4QU}{mr_0^2\Omega^2},\quad a_y=-a_x,\qquad
q_x=-\frac{2QV_{\rm RF}}{mr_0^2\Omega^2},\quad q_y=-q_x.
\]

Знак q меняется с RF phase; далее \(q=|q_x|\). При \(a=0\) первая область идеальной undamped Mathieu stability заканчивается около \(q=0.908\). Начальный диапазон \(0.1\le q\le0.3\) выбран для контролируемой pseudopotential approximation, а не как предел trapping. Окончательный контроль — Floquet multipliers и full-RF trajectories. Основа: [Leibfried et al., 2003](https://doi.org/10.1103/RevModPhys.75.281).

### 3.2. Дополнительная axial DC structure

Сравниваем: (i) по четыре небольших DC pads у каждого торца; (ii) короткие DC segments с тем же hyperbolic facing profile, отделённые зазорами от центральных RF electrodes. В первом engineering reference четыре центральных RF bodies непрерывны. DC structures дополняют, а не заменяют их.

Вблизи центра:

\[
\Phi_{\rm ax}=\frac{m\omega_z^2}{2Q}
\left(z^2-\frac{x^2+y^2}{2}\right).
\]

Laplace equation требует **radial defocusing**. Convention для geometry factor:

\[
\omega_z^2=\frac{2Q\kappa U_{\rm end}}{mz_0^2};
\]

\(\kappa\) вычисляется из 3D fields. После добавления axial potential:

\[
\begin{aligned}
a_x&=\frac{4QU}{mr_0^2\Omega^2}-\frac{2\omega_z^2}{\Omega^2},\\
a_y&=-\frac{4QU}{mr_0^2\Omega^2}-\frac{2\omega_z^2}{\Omega^2},\\
a_z&=\frac{4\omega_z^2}{\Omega^2},\qquad q_z=0.
\end{aligned}
\]

При \(U=0\), малом q:

\[
\omega_r^2\simeq\frac{q^2\Omega^2}{8}-\frac{\omega_z^2}{2}>0.
\]

В finite geometry RF fringing и asymmetry дают axial RF components; \(q_z=0\) больше не гарантировано. Определяем все normal modes, RF null и полную 6D monodromy matrix. DC alone не обеспечит устойчивое удержание свободного иона во всех направлениях этой архитектуры после shutdown. Standby RF возможно, но требует энергии.

### 3.3. Pseudopotential, micromotion, глубина

Для peak RF field amplitude и separation of timescales:

\[
\Psi=\frac{Q^2|\mathbf E_{\rm RF}|^2}{4m\Omega^2},\qquad
U_{\rm eff}=\Psi+Q\Phi_{\rm DC}+U_{\rm image}.
\]

В идеале \(|\mathbf E_{\rm RF}|=V_{\rm RF}\rho/r_0^2\),

\[
\Psi=\frac{Q^2V_{\rm RF}^2\rho^2}{4m\Omega^2r_0^4},
\qquad \omega_{r,\rm RF}=\frac{q\Omega}{\sqrt8}.
\]

Intrinsic micromotion amplitude относительно secular amplitude порядка \(q/2\). Stray field смещает ион на \(x_s\simeq QE_s/(m\omega_r^2)\), создавая excess micromotion порядка \(qx_s/2\). RF phase imbalance также создаёт micromotion.

Разделяем три величины:

- pseudopotential_depth: lowest escape saddle finite 3D \(U_{\rm eff}\) относительно минимума, с отдельным учётом absorbing surfaces;
- acceptance_energy: энергия до выбранной безопасной границы;
- dynamic_escape_threshold: фазозависимая граница удержания full-RF dynamics.

В бесконечном harmonic pseudopotential нет собственного finite escape saddle. Для screening полезна энергия на \(\rho=\eta r_0\):

\[
D_\eta=\Psi(\eta r_0)=
\frac{Q^2V_{\rm RF}^2\eta^2}{4m\Omega^2r_0^2}
=\frac{q|Q|V_{\rm RF}\eta^2}{8}.
\]

Это radial RF acceptance без axial defocusing, image attraction и imperfections. Она не равна рассчитанной глубине реального устройства и не гарантирует lifetime.

## 4. Scaling laws и начальный nanometer screening

### 4.1. Dimensionless parameters

При \(\tilde{\mathbf r}=\mathbf r/r_0,\ \tilde t=\Omega t\) кроме \(a_i,q_i\) нужны:

\[
\frac L{r_0},\frac t{r_0},\frac g{r_0},\frac{r_{\rm edge}}{r_0},
\frac{\omega_z}{\Omega},\frac\gamma\Omega,\frac{m_{\rm gas}}m,
\frac D{k_BT},\frac{\hbar\omega_r}D,\frac{x_{\rm zpf}}{r_0},
\frac{\ell_{\rm rough}}{r_0},\frac{\ell_{\rm patch}}{r_0},
\frac{\Omega\varepsilon}{\sigma},\Omega R_{\rm feed}C_{\rm load},
\frac t\delta,\frac{\Omega L_{\rm circuit}}c.
\]

g — fabrication gap; \(\gamma\) — momentum damping rate;
\(\delta=\sqrt{2/(\mu\sigma\Omega)}\). Для stochastic dynamics добавляем \(\dot E_{\rm noise}/(D\Omega)\), correlation time и full noise spectrum. White-noise intensity нельзя считать неизменной при scaling.

### 4.2. Семейство scaling laws

При геометрическом подобии \(r_0\propto s,\ V_{\rm RF}\propto s^p\) и fixed q:

\[
\Omega,\omega_r\propto s^{p/2-1},\qquad
D_\eta\propto s^p,\qquad E_{\rm surface}\propto s^{p-1}.
\]

| Сохраняемое условие | p | RF frequency | Acceptance energy | Surface field |
|---|---:|---:|---:|---:|
| Voltage и q | 0 | \(s^{-1}\) | constant | \(s^{-1}\) |
| Surface field и q | 1 | \(s^{-1/2}\) | s | constant |
| Frequency и q | 2 | constant | \(s^2\) | s |

При fixed voltage и frequency получаем \(q\propto s^{-2}\): уменьшение размера выводит из исходной области stability. Сохранение dimensionless trajectory не сохраняет barrier относительно температуры.

Для подобных проводящих тел \(C\propto s,\ R\propto\rho_{\rm Si}/s\). С peak phasors:

\[
P_{\rm ohmic}\sim\frac12R(\Omega CV)^2\propto s^{3p-1},
\quad
P_{\rm diel}\sim\frac12\Omega CV^2\tan\delta_d\propto s^{5p/2}.
\]

При fixed voltage/fixed q идеализированные ohmic losses растут как \(s^{-1}\). Circulating reactive power не равна потреблению источника; нужны resonator quality factor и RF generator efficiency. Pad capacitance, contacts, feedlines, skin effect и dispersive conductivity нарушают простые laws. Thermal resistance подобного тела обычно \(\propto s^{-1}\), но suspended MEMS требует собственной heat equation.

Power/geometry scaling уже исследовались: [Nguyen et al., 2025](https://doi.org/10.3390/e27060576), [Sterk et al., 2026](https://doi.org/10.1116/5.0332336). Это не самостоятельная новизна проекта.

### 4.3. Оценки порядка величин — без FEM

Пример: \(^{40}\mathrm{Ca}^+\), масса \(40u\), q=0.2, differential peak voltage 10 V, U=0, без axial defocusing:

| \(r_0\) | \(f_{\rm RF}\) | \(f_{r,\rm RF}\) | Field у вершины | \(x_{\rm rms}\), 300 K |
|---:|---:|---:|---:|---:|
| 10 µm | 247 MHz | 17.5 MHz | 1 MV/m | 2.27 µm |
| 1 µm | 2.47 GHz | 175 MHz | 10 MV/m | 227 nm |
| 100 nm | 24.7 GHz | 1.75 GHz | 100 MV/m | 22.7 nm |
| 10 nm | 247 GHz | 17.5 GHz | 1 GV/m | 2.27 nm |

\[
f_{\rm RF}=\frac1{2\pi r_0}\sqrt{\frac{2|Q|V_{\rm RF}}{mq}},
\qquad x_{\rm rms}=\sqrt{\frac{k_BT}{m\omega_r^2}}.
\]

Температура здесь — предполагаемая secular thermal equilibrium; её существование ещё нужно обосновать cooling model. Эти строки не являются рабочими точками устройства. Во всех строках \(D_1=0.25\) eV, а \(D_{0.5}=0.0625\) eV \(\simeq2.42k_BT\) при 300 K.

Предварительное условие \(D_\eta\ge\chi k_BT\) даёт

\[
V_{\rm RF}\ge\frac{8\chi k_BT}{q|Q|\eta^2}.
\]

При \(\chi=20,\eta=0.5,T=300\) K нужно 82.7 V в идеализированном примере. **20 \(k_BT\) — screening margin**, не универсальный lifetime law. RF heating может создавать nonthermal tails.

Если \(E_{\max}=g_EV_{\rm RF}/r_0\), где \(g_E\ge1\) описывает finite geometry, то

\[
r_0\gtrsim\frac{8g_E\chi k_BT}
{q|Q|\eta^2E_{\rm allow}}.
\]

Для **иллюстративного**, не установленного для silicon, \(E_{\rm allow}=10^7\) V/m, \(g_E=1\) получаем \(r_0\gtrsim8.27\) µm. Это условная граница конкретного теплового режима. Иной field limit, cooling, q или lifetime меняют её. При 4 K thermal часть уменьшилась бы в 75 раз, но cryogenics меняет system footprint и power.

### 4.4. Нарушения classical similarity

1. **Surface noise:** наблюдаемое в некоторых microtraps scaling близко к \(d^{-4}\), но не является законом для silicon при 100 nm. Нужны patch correlation length и frequency dependence. [Sedlacek et al., 2018](https://doi.org/10.1103/PhysRevA.97.020302).
2. **Image attraction:** у conducting plane \(U_{\rm image}=-Q^2/(16\pi\varepsilon_0d)\). Для e: 0.00360 eV по модулю при 100 nm, 0.0360 eV при 10 nm. Это одна плоскость, не сумма четырёх гипербол. Симметрия обнуляет force в центре, но не destabilizing curvature. Для точной геометрии нужен Green function; self-energy содержит factor 1/2.
3. **Quantum motion:** \(x_{\rm zpf}=\sqrt{\hbar/(2m\omega_r)}\). При fixed voltage/fixed q отношение \(x_{\rm zpf}/r_0\propto s^{-1/2}\). При 100 nm в примере \(x_{\rm zpf}\simeq0.269\) nm: ground-state size сам по себе ещё не запрещает trapping. Проверяются \(\hbar\omega/D\), surface tunneling и quantum noise.
4. **Emission/breakdown:** microscopic enhancement, oxide и adsorbates важнее средней \(V/g\). Paschen law нельзя продолжать в vacuum/nanogap без emission model. Есть экспериментальные данные 20–800 nm gaps, но они не задают допустимое RF поле нашей поверхности. [Wang et al., 2022](https://doi.org/10.1116/6.0002160).
5. **Silicon:** при thickness порядка depletion/screening length или mean free path continuum conductivity требует пересмотра. Для degenerate doping классический Debye length может быть неверен; нужны Thomas–Fermi screening, carrier relaxation и dispersive response.

**Предварительный вывод:** Mathieu confinement на nanometer scale математически не исключён. Комнатный fully integrated RNG на десятках nm пока не обоснован. Ожидаемый конфликт depth–field–frequency–surface–readout является проверяемой гипотезой, не universal no-go theorem.

**Классификация границ miniaturization:** идеальная Mathieu/Floquet stability (F при заданных полях, q и отсутствии потерь) сама по себе не задаёт нижнего размера; thermal acceptance (D_\eta/\chi k_BT) — R, где \(\chi\) является выбранным требованием, а не физической константой; field emission, breakdown, conductivity, roughness и contacts — T с пока широкими, конкретно не установленными для fabricated silicon nanogaps диапазонами; image attraction и quantum/surface corrections — F для заданного conductor model, но численная применимость к четырём finite branches является U до Green-function/material model; noise spectrum, collision lifetime, stochastic first passage — U до получения соответствующих \(S_E(f,d,T)\), cross-sections и tails. Поэтому выражения \(r_{0,\min}\) в этом разделе — **conditional scenario boundaries**, не universal physical minimum.

## 5. Материалы, поля, RF losses и тепло

Последовательно используем три уровня:

1. **Electrostatic basis fields:** \(\nabla\cdot(\varepsilon\nabla\phi_j)=0\), unit voltage на каждом отдельном conductor. Superposition даёт RF, DC и weighting fields.
2. **Electroquasistatic finite conductivity:** при phasor convention \(e^{i\Omega t}\),
   \[
   \nabla\cdot[(\sigma+i\Omega\varepsilon)\nabla\tilde\phi]=0.
   \]
   На interfaces соблюдаются potential continuity и normal total-current continuity; contacts получают реальные impedances. Equipotential RF faces теперь результат расчёта, а не навязанное свойство silicon.
3. **Dispersive/full EM:** при нарушении критериев выше — frequency-dependent \(\sigma(\Omega,T,N_D)\), eddy currents/full Maxwell, circuit ports и dielectric dispersion. Нанометрический зазор не делает весь GHz feed circuit автоматически quasistatic.

Диагностические screening limits: \(\Omega\varepsilon/\sigma<0.1\), \(\Omega R_{\rm feed}C_{\rm load}<0.1\), \(t/\delta<0.3\), \(\Omega L_{\rm circuit}/c<0.1\). Это ориентиры выбора model hierarchy, не универсальные sufficient conditions. Расчёт большей сложности нужен, если ошибка EQS превышает 5% в требуемых observables.

Peak-phasor power и steady heat equation:

\[
P_{\rm Si}=\frac12\int_{\rm Si}\Re\sigma\,|\tilde{\mathbf E}|^2dV,
\quad
P_{\rm diel}=\frac{\Omega}{2}\int\varepsilon'\tan\delta_d
|\tilde{\mathbf E}|^2dV,
\]

\[
-\nabla\cdot[k(T)\nabla T]=p_{\rm Si}+p_{\rm diel}+p_{\rm other}.
\]

В transient loading добавляется \(\rho_m c_p\partial_tT\). Проверяем thermal anchors, electrode deformation, temperature-dependent resistivity и положительную electrothermal feedback. Обязательно считаем contacts, dielectric loss в supports, generator, resonator и sensing circuit.

Для field emission применяем сначала Fowler–Nordheim/Murphy–Good model как screening с uncertain barrier, oxide и enhancement factor; для semiconductor нужны carrier-supply ограничения. Рассчитываем total emission current и induced charging, а не только threshold field. Breakdown limit — envelope по материалу, gap, pressure, частоте и duty cycle; «универсальная прочность vacuum» не вводится.

Проверяем механическую жёсткость тонких electrodes, Maxwell stress порядка \(\varepsilon_0 E^2/2\), thermal stress, roughness и DRIE scalloping. Изменение geometry под voltage может испортить quadrupole раньше разрушения. Свободно задаваемый lithographic feature size не доказывает возможность получить suspended four-branch structure с нужными contacts и surface quality.

## 6. Operating environment: A, B, C

### 6.1. Общая collision model

Для газа j:

\[
n_j=\frac{p_j}{k_BT_g},\qquad
\nu_j=n_j\langle\sigma_j(v_{\rm rel})v_{\rm rel}\rangle,
\]

\[
\lambda_{\rm loss}=\sum_jn_jk_{{\rm loss},j}
+\lambda_{\rm surface}+\lambda_{\rm technical}.
\]

Нельзя приравнивать total collision rate к loss rate: elastic collisions могут охлаждать, нагревать или почти не влиять. Для charge exchange, chemical reactions и cluster formation нужны отдельные channels и literature-based cross sections. Full-RF collision simulation использует differential scattering, RF phase и mass ratio. Основа: [Viehland & Goeringer, 2004](https://doi.org/10.1063/1.1691405).

Langevin approximation для области частых малых momentum transfers:

\[
m\ddot{\mathbf r}=-Q\nabla\Phi-m\gamma\dot{\mathbf r}
+\mathbf F_{\rm th},
\quad
\langle F_{{\rm th},i}(t)F_{{\rm th},j}(t')\rangle
=2m\gamma k_BT_g\delta_{ij}\delta(t-t').
\]

\(\gamma\) и diffusion связаны fluctuation–dissipation theorem, не выбираются независимо ради желаемой энтропии. White Gaussian model проверяется против Boltzmann/Monte Carlo и не используется для редких больших collisions.

### 6.2. Configuration A — ambient air

При 1 atm collision-free Mathieu diagram недостаточна. При illustrative mobility \(\mu_{\rm ion}=2\times10^{-4}\) m²/(V·s), характерной по порядку для малых molecular air ions, для массы \(40u\):

\[
\gamma=\frac{|Q|}{m\mu_{\rm ion}}\simeq1.2\times10^{10}\ {\rm s}^{-1},
\quad D_{\rm diff}=\frac{\mu_{\rm ion}k_BT}{|Q|}
\simeq5.2\times10^{-6}\ {\rm m^2/s}.
\]

Это **оценочное допущение**, не измеренная mobility Ca+ в воздухе. Chemistry может быстро изменить species и mobility. В отсутствие confinement diffusion time через 100 nm порядка \(r_0^2/(2D_{\rm diff})\sim1\) ns; это не рассчитанный trapped lifetime.

В high-frequency averaging для constant drag damping RF pseudopotential подавляется:

\[
\Psi_\gamma\simeq
\frac{Q^2|\mathbf E_{\rm RF}|^2}{4m(\Omega^2+\gamma^2)}.
\]

Формула применима при малом вынужденном смещении и separation of timescales; её нельзя считать универсальной заменой kinetic model. В overdamped limit нужны drift–diffusion/Fokker–Planck и absorbing surfaces. Учитываются water clusters, oxygen, charge exchange, discharge и surface neutralization.

**Приоритет:** дешёвый analytical rejection test. Данные о trapping dust/microparticles в air не подтверждают trapping одного atomic ion. Ambient configuration сохраняется в сравнении, но не является начальным engineering reference.

### 6.3. Configuration B — controlled gas

Наиболее естественный первый buffer gas — He из-за малого mass ratio для Ca+, Sr+, Yb+; окончательный выбор требует collision/reaction data. Возможны Ne и mixtures. Даже elastic collisions в RF trap вызывают micromotion interruption и nonthermal power-law tails: thermal energy alone недостаточна для оценки потерь. [DeVoe, 2009](https://doi.org/10.1103/PhysRevLett.102.063001).

Искомое окно давления ограничено снизу необходимым damping/entropy refresh и сверху heating, chemistry, loss и readout broadening. Pressure не оптимизируется независимо от trap depth и loading contamination. Нужны stationary distribution и first-passage times, а не только effective temperature.

### 6.4. Configuration C — sealed vacuum microcavity

Первый engineering reference: MEMS/wafer-bonded cavity вокруг trap, без vacuum для всего processor. На стадии 0 pressure рассматривается как параметр, например от \(10^{-9}\) до \(10^2\) Pa для C/B, плюс отдельный ambient endpoint. Это диапазон аналитического screening, не заказ большой сетки.

Из требуемой survival probability при constant Poisson hazard:

\[
S(t)=e^{-\lambda_{\rm loss}t},\qquad
p_{\max}=\frac{-k_BT_g\ln S_*}{k_{\rm loss}t_*},
\]

для одного gas channel без прочих потерь. Например, \(S_*=0.99\) за 1 h требует \(\lambda\le2.79\times10^{-6}\) s\(^{-1}\). При **предполагаемом** \(k_{\rm loss}=10^{-18}\) m³/s это даёт \(p_{\max}\simeq1.16\times10^{-8}\) Pa (\(1.16\times10^{-10}\) mbar). Изменение k_loss на три порядка меняет допустимое pressure на три порядка. Это пример зависимости, не требование универсального UHV. Non-Poisson heating требует другой survival model.

Для каждого gas species решаем gas-load balance:

\[
V_c\frac{dp_j}{dt}=
Q_{{\rm leak},j}+Q_{{\rm outgas},j}
+Q_{{\rm source},j}+Q_{{\rm permeation},j}
-S_jp_j,
\]

где Q измеряются в Pa·m³/s, S — m³/s. Getter имеет конечную capacity и species-dependent speed, не удаляет одинаково все газы; He особенно проблематичен для обычного getter. При similar cavity surface/volume растёт как \(1/s\), поэтому начальный outgassing pressure rise ухудшается при scaling. Включаем bake/activation temperature, bonding seals, contaminant inventory, thermal cycles и долгосрочную permeation.

Есть экспериментальный precedent sealed package около 1 cm³ с linear Paul trap, Yb sources, getter и He buffer gas; это не nanocavity: [Schwindt et al., 2016](https://doi.org/10.1063/1.4948739). Для 3 cm³ packages сообщены trapping lifetimes до 50 дней после хранения sealed до 10 лет: [Thrasher & Schwindt, 2025](https://doi.org/10.1063/5.0243517). Такие результаты подтверждают принцип автономной упаковки, но не её scaling до mm³.

| Environment | Потенциальная польза | Главный риск | Начальное действие |
|---|---|---|---|
| A: air | Простая внешняя среда | Collisions, chemistry, discharge, surface loss | Mobility/kinetic upper-bound screening |
| B: controlled gas | Damping и entropy refresh | RF heating, tails, газовая purity | Найти pressure–mass-ratio window |
| C: sealed vacuum | Меньше collisions и химии | Getter, leakage, слабый entropy refresh | Основной system reference |

## 7. Ион, загрузка, shutdown и recovery

### 7.1. Species shortlist

| Species | Преимущества | Издержки |
|---|---|---|
| \(^{40}\)Ca+ | Стандартный reference, известные transitions/collisions | Источник Ca, photoionization, cooling/repump optics |
| Sr+ | Известная ion photonics, возможна optical detection | Больше масса, отдельные sources/lasers |
| Yb+ | Micro-oven и sealed-clock precedents, heavy ion полезен с He | UV optics, isotope/chemistry, source contamination |
| Mg+/Be+ | Higher secular frequency при тех же voltage/geometry | Deep UV; handling/process constraints |
| Molecular/air ions | Потенциально простая ionization газа | Химическая нестабильность и трудное species-selective readout |
| Electron | Лёгкий заряд, иные image-current возможности | Это другая архитектура, резко иные RF/EM и surface limits |

Первый comparative analytical screening: Ca+, Sr+, Yb+; Ca+ — convention reference, не заранее выбранный победитель. Multiply charged ions могут повысить force, но ухудшают source/charge-exchange problem и не считаются бесплатным улучшением.

### 7.2. Две loading architectures

**L1 — neutral micro-oven + photoionization.** Отдельная thermally isolated reservoir chamber, shutter/collimator, backside или axial loading aperture, ionization light в capture volume, затем cooling. Source material и evaporated deposit не должны перекрывать electrode gaps и менять work function.

Есть primary demonstration micro-heater до \(90\times90\) µm², response neutral fluorescence — десятки ms при milliwatt-scale power. Это **не доказательство** single-ion capture за те же ms, не размер всего source assembly и не lifetime test на 10 лет: lifetime источника в работе оценивался расчётно. [Kumar et al., 2025](https://doi.org/10.1063/5.0271048).

**L2 — neutral reservoir + field-emission electron-impact ionizer.** MEMS emitter, gate, electron energy control, collector и shielding. Убирает часть photoionization optics, но добавляет high voltage, electron-current noise, contamination, surface charging и finite emitter lifetime. Ionizer желательно размещать в отдельной loading zone; затем ion transport в RNG zone. Возможность fully integrated single-ion architecture — hypothesis.

Laser ablation — контрольный, но слабый кандидат для плотного CPU module из-за pulsed energy, debris и optical overhead. Внешний macroscopic source допустим для comparison, но такая система не называется fully integrated RNG.

### 7.3. Capture и управление загрузкой

Для dilute neutral volume model:

\[
R_{\rm load}=\int n_a(\mathbf r)\,
\Gamma_{\rm ion}(\mathbf r)\,
P_{\rm cap}(\mathbf r,\mathbf v,\varphi_{\rm RF})\,
f(\mathbf v)\,d^3r\,d^3v.
\]

В photoionization saturation и electron-impact models \(\Gamma_{\rm ion}\) рассчитывается отдельно. \(P_{\rm cap}\) определяется initial phase-space, RF turn-on ramp, actual depth, collisions и cooling. Нельзя заменять его единицей.

Для Poisson loading \(P(N_{\rm captured}\ge1)=1-e^{-R_{\rm load}t}\) — только первая оценка. Для **ровно одного** иона нужны detection feedback, выключение source и учёт multiple capture. При fixed neutral density/ionization rate capture volume примерно \(\propto r_0^3\), что ухудшает rate при miniaturization; focused injection может изменить scaling, но требует отдельной optics model.

State machine:

OFF → RF_RAMP → LOAD → VERIFY_ONE_ION → THERMALIZE → RNG_READY;
loss/health failure → OUTPUT_INHIBIT → RELOAD;
repeated failure → FAULT.

Учитываем false ion/no-ion decisions, multiple ions, dark molecular ions, timeout, источник в off-state, restart после outages и depletion inventory за \(10^4\) или больше cycles. Биты во время uncertain occupancy не выдаются.

### 7.4. Как loading увеличивает размер

Минимальная heater area опубликованного примера — 8100 µm² = 0.0081 mm², без reservoir, anchors, shielding, aperture и ionization. Если bare nano RF structure занимает 1–100 µm², один heater даёт увеличение area порядка 80–8000 раз; это **условное сравнение**, не layout calculation.

Полный source assembly разумно сначала параметризовать в sub-mm/mm range, а не присваивать ему nanometer footprint. Для module оценивать отдельно reservoir volume, optical path, heat isolation, routing, getter и seal width. Demonstrated cm³ sealed packages — технологический precedent; 1 mm³ package для этой схемы пока assumption.

Source energy example: budget 20 mW × 50 ms = 1 mJ на heating pulse. Это иллюстрация, не измеренная энергия single-ion initialization. При повторном запуске для каждых 256 bits только этот pulse дал бы 3.9 µJ/bit. Реальная startup energy дополнительно включает ionization, cooling и verification.

## 8. Entropy model: что является источником случайности

### 8.1. Полная particle model

\[
m\ddot{\mathbf r}=
-Q\nabla[\Phi_{\rm RF}(\mathbf r,t)+\Phi_{\rm DC}(\mathbf r)]
-\nabla U_{\rm image}
+\mathbf F_{\rm gas}
+Q\mathbf E_{\rm stochastic}
+\mathbf F_{\rm backaction}.
\]

Допускаем deterministic RF amplitude/phase drift и stochastic components с раздельным происхождением. Surface collisions — absorbing/reaction boundaries. Coulomb interactions вводятся только для explicit multi-ion control, не для замены single-ion reference.

| Candidate entropy | Физическая основа | Model | Ограничение и приоритет |
|---|---|---|---|
| Johnson–Nyquist circuit noise | Диссипативные резистивные элементы при T>0; FDT задаёт спектр при известном impedance | Thevenin/Norton network + electrode transfer/weighting fields; затем ion response | **Приоритет 1 для раннего quantitative control**: самая ограниченная физическая модель, но измеряемый circuit noise сам по себе не преимущество ловушки; слабая motional coupling может требовать охлаждения/усиления и нагревать ион |
| Stochastic ion–buffer-gas collisions | Случайные времена/параметры столкновений с молекулярными скоростями и differential cross-sections | Boltzmann collision operator или validated collision MC; Langevin только в проверенном пределе | **Условный кандидат 1 для ion-mediated entropy**: прямая связь с движением и потенциально сильный импульсный сигнал; rate и энергия обмена зависят от газа/массового отношения; RF micromotion может нагревать, давать tails и потери |
| Silicon surface charge / electric-field fluctuations | Зарядовые дефекты, adsorbates, interface traps, fluctuating patches | Измеренный one-sided $S_E(f)$, spatial correlations, telegraph/patch model только при данных | **Условный кандидат 2**: прямое силовое воздействие на ион и известная связь с motional heating, но спектр конкретного doped-Si/oxide при GHz secular frequency и nm–µm расстояниях не известен; максимальный риск heating/lifetime и высокая неопределённость |
| Quantum state preparation/measurement | Born rule и spontaneous emission | Master equation + quantum trajectories + detector POVM | Чёткое quantum происхождение, но lasers/reset/readout overhead |
| Nonlinear chaos + bath | Amplification и mixing physical noise | Full RF, nonlinear fields, Lyapunov exponents + noise | Chaos alone не создаёт условную энтропию |
| Numerical roundoff/solver error | Алгоритм | Только numerical error model | Никогда не entropy source |

В ideal single-ion quadrupole equations линейны и deterministic. Для известного initial state и drive неопределённости нет. Chaotic trajectory при известном exact initial state также deterministic. Positive Lyapunov exponent характеризует sensitivity, но не доказывает physical randomness или cryptographic unpredictability.

Приоритеты ниже — для **моделирования и controls**, не утверждение что один источник уже выбран. Сначала оценить Johnson-noise feedthrough напрямую при ion-absent и circuit-only controls: если он доминирует, ион может быть избыточен. Затем газ как наиболее понятный путь к реальным случайным импульсам, если pressure/loss окно не пусто. Surface noise оставить отдельным conditional branch: опубликованные ion-heating datasets охватывают примерно $10^{-15}$–$10^{-6}\,\mathrm{V^2m^{-2}Hz^{-1}}$, но в основном MHz и ion–surface distances десятки микрометров–миллиметры; это **не** prior для нашего Si trap на 1 µm/GHz. Обзор подчёркивает множественность noise mechanisms и геометрическую зависимость; экспериментальные distance exponents различаются. Ни частоту корреляции, ни Hmin rate нельзя получить из одной PSD amplitude без time-domain/conditional-source model. См. Brownnutt et al., DOI [10.1103/RevModPhys.87.1419](https://doi.org/10.1103/RevModPhys.87.1419); для gas kinetics — Viehland & Goeringer, DOI [10.1063/1.1691405](https://doi.org/10.1063/1.1691405).

### 8.2. PSD conventions и thermal consistency

Используем one-sided PSD \(S_E^{(1)}(f)\) по frequency в Hz, units (V/m)²/Hz, \(f>0\), и variance \(\int_0^\infty S(f)\,df\). В weak-noise harmonic limit:

\[
\dot{\bar n}=\frac{Q^2S_E^{(1)}(\omega_r/2\pi)}
{4m\hbar\omega_r},\qquad
\dot E=\frac{Q^2S_E^{(1)}(\omega_r/2\pi)}{4m}.
\]

RF sideband noise/parametric heating добавляются отдельно. Это не complete loss law и не универсальная entropy formula. Для circuit bath:

\[
S_V^{(1)}(f)=4k_BT\,\Re Z(f)
\]

в classical regime; при \(\hbar\omega\sim k_BT\) выбираем quantum spectrum с явно указанным ordering, emission/absorption и detailed balance. Не подставляем symmetric quantum noise механически в classical heating formula.

Surface model:

\[
S_E(f,d)=S_0(d_0/d)^\beta(f_0/f)^\alpha
\]

используется только в измеренном диапазоне или как **сценарная гипотеза** с cutoff и correlation length. Для degenerate silicon опубликовано около \(10^{-10}\) (V/m)²/Hz при d=40 µm, f=1.125 MHz: [Britton et al., 2009](https://doi.org/10.1063/1.3254188). Это anchor конкретного образца, а не свойство всех silicon surfaces и не основание экстраполировать на 100 nm.

Особый риск: noise bath без damping/cooling нагревает ион до loss и не образует stationary RNG. Необходимо обосновать heat removal через gas, resistive cooling, laser cooling или feedback, включая их noise, power и backaction.

### 8.3. Честная трактовка вычислительной энтропии

В simulation все stochastic drivers генерируются seeded PRNG для воспроизводимости. Если seed и алгоритм известны, сам simulated bitstream не является новым physical entropy source. Из simulation можно оценить **условную энтропию гипотетического физического устройства при принятой bath model**.

Нужно отдельно сравнить:
(i) trajectory/readout с physical bath;
(ii) отключённый bath при тех же numerical settings;
(iii) ion absent, тот же amplifier;
(iv) deterministic injected interference;
(v) amplifier noise-only;
(vi) measurement backaction-only.

Если ион не добавляет наблюдаемой энтропии относительно прямого измерения того же noise bath, это отрицательный результат для ion-based RNG advantage, даже если полный выход выглядит случайным.

Trapped-ion quantum RNG уже существуют: [Pironio et al., 2010](https://doi.org/10.1038/nature09008), [Um et al., 2013](https://doi.org/10.1038/srep01627), [Liu et al., 2025](https://doi.org/10.1038/s41586-025-08737-1). Bell/contextuality и processor protocols не переносятся автоматически на один classical motional sensor.

## 9. Реальный readout и преобразование в bits

### 9.1. Image-current architecture

Для каждого sensing terminal решаем weighting problem: \(\phi_{w,k}=1\) на этом terminal, 0 на остальных при отключённых независимых источниках; \(\mathbf E_{w,k}=-\nabla\phi_{w,k}\), units m\(^{-1}\). По Shockley–Ramo theorem:

\[
i_k(t)=Q\dot{\mathbf r}(t)\cdot\mathbf E_{w,k}[\mathbf r(t)]
\]

с выбранной terminal sign convention. Это induced current, не ток поглощения иона. Основа: [Ramo, 1939](https://doi.org/10.1109/JRPROC.1939.228757).

**Критическая симметрия:** при жёстком электрическом объединении противоположных electrodes weighting potential пары even относительно центра; её gradient в центре нулевой. Linear dipole pickup подавляется, остаётся quadratic response и harmonics. Нельзя считать single-electrode current измеренным, если circuit суммирует противоположные electrodes.

Сравниваем два варианта:

- отдельные axial DC sensing pads с differential pickup, RF filtering и ненулевым weighting gradient вдоль z;
- differential sensing противоположных RF electrodes через bias tees/diplexers, сохраняющие одинаковый RF drive на паре, но разделяющие sensing ports на secular frequency.

Реальная сеть меняет boundary conditions. Нужны multiport impedance/capacitance matrix, feed imbalance и влияние sense impedance на trapping. Axial pickup может иметь худший effective distance, чем r0; нельзя подставлять r0 автоматически.

Для harmonic motion amplitude A:

\[
I_{\rm pk}\sim\frac{|Q|\omega_s A}{d_{\rm eff}},
\qquad d_{\rm eff}^{-1}=|\mathbf E_w(\mathbf r_c)|.
\]

Пример \(Q=e,\ f_s=1\) MHz, A=0.1 µm, \(d_{\rm eff}=1\) µm даёт около 0.10 pA. Эти параметры иллюстративны и не совпадают с GHz example выше.

Для идеального resistor 1 MΩ при 300 K current-noise density \(\sqrt{4k_BT/R}\simeq0.129\) pA/\(\sqrt{\rm Hz}\); в 1 kHz band это около 4.1 pA rms до учёта resonator response. Higher impedance не бесплатна: capacitance ограничивает bandwidth, amplifier voltage noise преобразуется во входной ток, возникает backaction.

RF feedthrough: всего 1 fF coupling при 100 MHz и 10 V даёт \(\Omega CV\simeq6.3\) µA peak. Он находится на другой frequency, но перегрузка и нелинейное down-conversion делают rejection и dynamic range обязательными constraints. Простое сравнение RMS signal с шумом без front-end model недостаточно.

### 9.2. Readout transfer model

\[
v_{\rm ADC}(t)=(h_i*i_{\rm ion})(t)
+v_{\rm amp}(t)+v_{\rm feedthrough}(t)+v_{\rm interference}(t),
\]

\[
y_n=\mathcal Q_{\rm ADC}\{v_{\rm ADC}(nT_s+\delta t_n)\}.
\]

Включаем transimpedance, resonator bandwidth, anti-alias filter, amplifier current/voltage noise и correlation, finite ADC resolution, saturation, offsets, sampling jitter и RF phase locking. Impedance создаёт damping/backaction и thermal force noise по FDT; sensor не является невозмущающим.

\[
{\rm SNR}_B=
\frac{\int_B S_{v,\rm ion}(f)\,df}
{\int_B S_{v,\rm nonion}(f)\,df}
\]

служит screening metric; она не равна entropy. Для narrow-band noisy oscillator entropy refresh ограничивается envelope relaxation и measurement innovation rate, а не только большой secular frequency.

### 9.3. Optical architecture

Альтернатива: internal-state quantum preparation → rotation → fluorescence readout → reset, либо photon-emission timing. Нужны optical Bloch/master equations, collection efficiency, detector efficiency, background, dark counts, dead time, afterpulsing и ion-loss discrimination. Photon process не считается Poisson на всех timescales из-за antibunching и internal dynamics.

SPAD, monolithically integrated with surface trap, уже демонстрировал ion/no-ion fidelity 0.99 с average window 7.7 ms: [Setzer et al., 2021](https://doi.org/10.1063/5.0055999). Это precedent integrated sensing, **не** доказательство bit-rate нашего RNG и не готовая integrated laser/source system.

Для fluorescence:

\[
R_{\rm det}=\eta_{\rm collection}\eta_{\rm detector}R_{\rm scatter}
+R_{\rm background}+R_{\rm dark}.
\]

Подсчёт background/dark entropy как ion entropy запрещён. При optical quantum branch power включает lasers, stabilization, repumping, modulation и reset.

### 9.4. Entropy accounting и output words

Для discrete raw symbol X:

\[
H(X)=-\sum_xp_x\log_2p_x,\quad
H_\infty(X)=-\log_2\max_xp_x.
\]

Marginal entropy недостаточна. Используем conditional block min-entropy \(H_{\min}^{\epsilon_s}(X^N|E)\), где E содержит RF phase, drift, known electronics, past observations и выбранную attacker side information. Classical guessing probability:

\[
p_{\rm guess}(X|E)=\sum_e p(e)\max_xp(x|e).
\]

При quantum side information нужна соответствующая quantum conditional min-entropy; classical NIST estimators сами её не доказывают. Определяем threat/trust model до выбора estimator. Shannon entropy и variance не заменяют min-entropy.

Для properly seeded universal-hash extractor консервативный ориентир:

\[
m_{\rm out}\le H_{\min}^{\epsilon_s}(X^N|E)
-2\log_2(1/\epsilon_{\rm ext}),
\]

с точными conventions/security constants выбранной theorem. Seed независим от raw source; его роль и стоимость учитываются. Hash/conditioning не создаёт недостающую энтропию. Для NIST constructions применяем ограничения vetted conditioning components отдельно; произвольный hash не объявляется сертифицированным extractor.

\[
R_{\rm out}=\mathcal A\frac{m_{\rm out}}{T_{\rm acquisition}+T_{\rm processing}},
\qquad
E_{\rm useful}=
\frac{\int P_{\rm total}(t)dt}{N_{\rm accepted\ output\ bits}}.
\]

\(\mathcal A\) — operating availability, если downtime не включён уже в denominator. Двойной учёт downtime исключается. Health-test rejections уменьшают useful rate; processing throughput не увеличивает entropy rate.

8-, 32-, 64-, 128-, 256-bit words — packaging. Они требуют соответствующего **суммарного entropy budget**, а не иной trap geometry. При \(R_{\rm out}=10^5\) bit/s без block/startup overhead время накопления равно 80 µs, 0.32 ms, 0.64 ms, 1.28 ms, 2.56 ms. Реальный latency включает extractor block, restart и buffering. 256-bit output из 16 bits entropy не становится 256-bit entropy.

### 9.5. Statistical plan

Оцениваем raw data до conditioning: marginal distributions, autocorrelation/PSD, conditional prediction, nonlinear dependence, RF-synchronous bias и restart correlation. Нулевая autocorrelation не доказывает independence.

- SP 800-90B: noise-source model, IID/non-IID choice, conservative entropy estimators, restart assessment, Repetition Count Test, Adaptive Proportion Test и startup/continuous health monitoring. Для formal-like dataset планируем не менее \(10^6\) последовательных raw samples и restart matrix 1000×1000; точные процедуры/пороги проверяются по действующей публикации и errata.
- SP 800-22: battery на достаточно длинных raw/conditioned sequences; длины и число sequences выбираются под каждый test. Учитываются multiple tests и expected false rejection. Короткий 256-bit block непригоден для всей battery.
- Повторяем analysis при temperature, voltage, pressure, gain, species, aging scenarios; сохраняем rejected runs, а не только «удачные» sequences.
- Synthetic known-bias, Markov и deterministic controls проверяют сам pipeline.

Официальные источники: [SP 800-90B](https://csrc.nist.gov/pubs/sp/800/90/b/final), [SP 800-22 Rev.1a](https://csrc.nist.gov/pubs/sp/800/22/r1/upd1/final). На дату обзора страница 90B отмечает errata 2025 года, 800-22 — план revision; [NIST publications](https://csrc.nist.gov/Projects/random-bit-generation/publications) также указывает final 800-90C 2025 года. Проверку актуальной редакции повторить перед публикацией.

Statistical tests не доказывают cryptographic security. Для simulation успешный NIST result проверяет совместимость simulated process с выбранными tests, но не существование physical entropy в изготовленном устройстве.

## 10. Полная semiconductor integration и честное сравнение

### 10.1. System architecture

Полный module включает RF silicon electrodes, axial/compensation/sensing electrodes, cavity, supports, getter, source reservoir, ionizer, RF generator/resonator, DC supplies/filters, readout, ADC, controller, health tests и conditioner. Если нужны lasers или active vacuum pumps, их footprint и wall-plug power входят в total, включая внечиповые узлы.

\[
P_{\rm total}=P_{\rm RF,wall}+P_{\rm DC}+P_{\rm readout}
+P_{\rm control}+P_{\rm optical}+P_{\rm vacuum}
+\overline P_{\rm loading}.
\]

При finite campaign startup добавляется напрямую:
\[
E_{\rm useful}=
\frac{E_{\rm startup}+P_{\rm run}t_{\rm run}
+\sum E_{\rm reload}}{N_{\rm accepted\ bits}}.
\]

Нельзя одновременно суммировать wall-plug RF power и уже включённые в него losses как независимые потребления.

Проверяем silicon doping/contacts, thermal budget, oxide charging, UV damage, metal contamination от source, outgassing, wafer bonding, hermeticity, alignment между слоями, vibration/shock, CMOS supply domains, RF EMI, substrate coupling и temperatures processor. Первый диапазон 300–350 K — только proposed requirement; более горячий CPU environment рассматривается sensitivity scenario.

Made of silicon ≠ CMOS-compatible. Нужен process-flow audit: lithography/DRIE, doping/annealing, contacts, release, cleaning, bonding, getter activation и source deposition, с допустимыми temperatures/contaminants по этапам. Back-end MEMS/chiplet integration может быть реалистичнее monolithic front-end CMOS.

90-nm CMOS-foundry surface trap уже демонстрировался, используя top-metal electrodes: [Mehta et al., 2014](https://doi.org/10.1063/1.4892061). **90 nm — process node, не ion–electrode distance.** Это не доказательство four-hyperbolic silicon electrode geometry или integration всех RNG components.

### 10.2. Три разных minimum sizes

1. **RF structure:** bounding box четырёх finite electrodes с feeds/anchors; отдельно active facing area и total silicon volume.
2. **Trapping region:** \(2r_0\), usable radius \(\eta r_0\), axial safe extent и capture volume; не lithographic node.
3. **Functional RNG module:** die area, package volume/height и external supporting volume, включая source/readout/cavity.

Для normal-offset geometry characteristic transverse RF span порядка \(2[r_0\cosh u_{\max}+t_{\max}]\), без anchors; axial length L. Эта формула даёт лишь geometry envelope. Total area нельзя оценивать как \(\pi r_0^2\).

Начальный module-size scenario — 1 mm² die area и 1 mm³ package; он **не считается established fabrication capability** для предложенного устройства. Проверенные sealed-ion precedents cm³-scale задают reference comparison. Требуемый размер остаётся результатом co-design.

### 10.3. Semiconductor TRNG baseline

Primary CMOS benchmark: latch TRNG в 130-nm technology, reported area 5343 µm² = 0.005343 mm², energy 0.116 pJ/bit при 0.3 V, включая on-chip von Neumann post-processing. [Wang et al., 2024](https://doi.org/10.1109/TVLSI.2023.3328602); [author institution record](https://waseda.elsevierpure.com/en/publications/a-0116-pjbit-latch-based-true-random-number-generator-featuring-s/).

Это опубликованная метрика конкретной operating point, не independently verified entropy-certified system baseline при всех PVT conditions. Для окончательного сравнения требуются полный paper, common entropy definition, throughput, power boundary и health-test overhead.

Иллюстрация: 1 mW total при 100 kbit/s даёт 10 nJ/useful bit — примерно \(8.6\times10^4\) больше указанной 0.116 pJ/bit. Это budget example, **не prediction** мощности ion RNG. Он показывает, насколько сильным должно быть иное преимущество для конкуренции.

Возможные преимущества для проверки: количественно обоснованная quantum entropy certification, независимость от substrate noise или использование уже существующего trapped-ion subsystem. Для standalone CPU RNG преимуществ energy/density/reliability пока не показано. Если измеряется только Johnson noise, прямой silicon noise sampler должен быть обязательным контролем.

## 11. Feasibility map и неопределённости

Параметры \(\theta\): geometry, species, RF/DC, conductivity/contact model, pressure/composition, temperatures, surface-noise model, source, readout и workload. Неопределённые параметры \(\xi\): roughness, enhancement factor, leak/outgas, cross sections, noise amplitude/correlation, circuit noise и drift.

Для каждого ограничения вводим normalized margin \(M_i(\theta,\xi)\), где \(M_i\ge1\) означает выполнение; для верхних пределов используем limit/value. Например:

\[
M_D=\frac{D_{\rm accessible}}{\chi k_BT_{\rm eff}},
\quad M_E=\frac{E_{\rm allow}}{E_{\max}},
\quad M_R=\frac{R_{\rm out}}{R_*},
\quad M_P=\frac{P_*}{P_{\rm total}}.
\]

\(M_D\) используется только как screening; nonthermal lifetime проверяется отдельно. Общая feasible set:

\[
\mathcal F=\{\theta: M_i(\theta,\xi)\ge1\ \forall i\}.
\]

Показываем optimistic, central и conservative scenarios. «Unknown» — отдельный статус, не «feasible». Robust feasibility означает выполнение для declared uncertainty envelope; вероятностный вывод допустим только при обоснованных distributions, а не при произвольных priors.

Нужны карты: r0–frequency, r0–voltage, r0–pressure, depth–readout bandwidth, module area–energy/bit; дополнительные dimensions фиксируются явно. Каждая boundary содержит uncertainty и dominant violated constraint. Active constraint может меняться при r0; возможны nonmonotonic feasible islands.

Минимальный RF size, minimum aperture и minimum module size ищутся раздельно на одной \(\mathcal F\). Оценка «минимум 100 nm» без specs/temperature/lifetime/readout/uncertainty не считается результатом.

Отрицательный результат классифицируется как: fundamental under stated physics; model-limited; technology-limited; requirement-limited; или unresolved due to missing data. Отсутствие литературы о noise на nm — основание для неопределённости, не доказательство невозможности.

Рабочая boundary register (значения — лишь уже оговорённые сценарии; перед итоговыми картами добавить source/uncertainty interval для каждого входа):

| Observable boundary | Classification | Present bound / uncertainty |
|---|---|---|
| Ideal Mathieu stable region | F | Для $a=0$ первая граница около $|q|=0.908$; размер не входит напрямую после nondimensionalization. Не включает collisions/end-loss/surface. |
| Thermal acceptance vs E-field | R + T + U | При q=.2, $\eta=.5$, 300 K и $\chi=20$ получается V≈82.7 V; при assumed $E_{allow}=10^7$ V/m и $g_E=1$ — r0≈8.27 µm. Порог $\chi$ — workload/lifetime proxy, field limit — неподтверждённая Si технологическая предпосылка; cooling меняет thermal term. |
| Field emission / breakdown | T + U | Для proposed doped-Si/oxide face and rounded/truncated hyperbolas нет применимого validated interval; крайнее поле включает enhancement и edge radius, а breakdown зависит от pressure/surface state. Не указывать единственное $E_{max}$ как known. |
| Image attraction / quantum motion | F under conductor model + U for device | Plane image-energy scales are known analytically, but exact four-surface Green function and material screening/oxide response are absent; continuum image model becomes questionable at nm separation. |
| Readout SNR/size | R + T + U | $A=10$–500 nm, $d_{eff}=0.5$–5 µm and R/B ranges are sensitivity assumptions only. Terminal cancellation and transfer function may shift signal from near zero to a detectable value; no size boundary inferred. |
| Noise/collision lifetime | U + R | Surface $S_E$ spans nine decades in literature under nonmatching conditions; proposed local spectrum unknown. Collision coefficient uses 0.3–3×10$^{-15}$ m³/s scenario, loss probability/event unknown; survival target depends on workload. |
| Complete module minimum | T + R + U | Electrode die, cavity, source/ionizer, vacuum/getter, high-frequency readout and control have no integrated measured floor for this four-hyperbola concept. Any minimum area/volume remains unbounded until architecture and acceptance criteria are fixed. |

Диапазоны здесь являются explicit assumptions/known literature envelope; там, где нет численного measured interval, uncertainty считается unbounded for decision purposes. Следовательно, сейчас можно сообщить conditional scaling and competing constraints, но не численную universal minimum trap/module size.

## 12. Архивный checkpoint незавершённой сессии

**Архивное состояние на момент прерывания, 27 октября 2026.** Тогда пункты ниже ещё не были завершены. Они сохранены как история восстановления; их завершённые версии находятся в разделах 13–18 и в companion documents. Исследовательская реализация не запускалась.

### 12.1. Staged computational research plan — INCOMPLETE

Предполагалась последовательность: literature/specification screening → analytical/Floquet reference → validated 2D hyperbolic boundary-value problem → finite 3D RF/DC fields → stochastic particle dynamics → physical readout/entropy model → system feasibility map. Подробные deliverables, resource limits и gates ещё не оформлены.

Приоритеты следующей разрешённой сессии: сначала analytical depth–field–frequency constraints и readout observability; затем проверка loading/cavity overhead; только после утверждения плана — численная реализация. Большие sweeps и remote workloads не разрешены.

### 12.2. Quantitative go/no-go gates — INCOMPLETE / NOT AGREED

Численные критерии перехода между стадиями ещё не завершены и не согласованы с пользователем. Нужно разделить numerical validation tolerances, physical confinement criteria, readout/conditional entropy и competitiveness criteria.

Обсуждавшиеся при подготовке варианты, **не утверждённые требования**: field error порядка 0.1% для analytical benchmark; physical nonquadrupole field residual порядка 1% в заданной safe region; сравнение EQS/full EM с допустимым расхождением 5%; screening depth 20 kBT с отдельной проверкой nonthermal tails; явные targets для rate, power, startup и availability. Примеры 1 mW, 100 kbit/s, 1 mm² и 1 mm³ в предыдущих разделах — сценарии, не доказанные достижения.

Короткие stable trajectories не доказывают часовое/годовое удержание. Для rare-loss rate нужны confidence bounds и validated first-passage/rare-event methods; при zero losses за total exposure T Poisson upper 95% bound примерно 3/T. Эта методика ещё не разработана для данной модели.

### 12.3. Python project structure — INCOMPLETE

Планируется изолированный namespace для новой hyperbolic RNG задачи: configuration/materials, analytical models, geometry, fields, dynamics, readout, entropy, system budgets, validation tests и reproducible results. Точное дерево, interfaces и dependency lock пока не составлены. Существующие src/config/results не менять и не считать новой архитектурой.

Предварительные candidates: Python, NumPy/SciPy, Gmsh и scikit-fem для лёгкого начала; конкретная совместимость и dependency versions ещё не проверены. Full EM solver следует выбирать отдельно. Ничего не установлено и не реализовано в этой сессии.

### 12.4. Numerical validation/convergence suite — INCOMPLETE

Сохранённый checklist для будущего развёртывания:

- units и differential/peak/RMS voltage conventions;
- точная topology четырёх ветвей, curved-boundary approximation, isolation и physical tags;
- manufactured potential test отдельно от physical equipotential-hyperbola geometry test;
- h/p refinement, внешняя граница, fit radius/order, derivative/Hessian convergence;
- full-domain reference до использования symmetry;
- Laplace trace, multipole symmetry, charge/current conservation и power balance;
- EQS/EM/material dispersion applicability;
- Mathieu/Floquet benchmark, timestep refinement, независимые integrators и RF phase;
- stochastic PSD/FDT, detailed balance и collision kinematics;
- absorbing boundaries, escape saddle/path и rare-loss confidence;
- weighting fields и cancellation в paired pickup;
- readout transfer/noise/aliasing/ADC и backaction;
- entropy estimators на известных synthetic controls, restart tests;
- ion-absent/no-bath/noise-only controls, полный energy/area accounting.

Пороги, test cases и acceptance protocol ещё предстоит определить.

### 12.5. Preliminary novelty assessment — PARTIAL / UNVERIFIED

Уже установлено: silicon traps, CMOS-foundry traps, trapped-ion RNG, geometry/power scaling, micro-ovens и sealed ion packages имеют precedents. Ни одна из этих идей сама по себе не является заявляемой новизной.

Близкий источник, найденный, но ещё не включённый в развёрнутую novelty comparison:
Yang Li et al., *Simulation Study of a Hyperbolic Linear Ion Trap with Asymmetric Curvature Radii*, DOI [10.1021/jasms.5c00307](https://doi.org/10.1021/jasms.5c00307), online 16 декабря 2025, journal issue 2026. Изучены publisher metadata и доступный PMC text/excerpt. Работа посвящена modified hyperbolic LIT и mass-spectrometric ejection, а не complete silicon RNG. Abstract употребляет “experimental results”, хотя заголовок/методика говорят о simulation; независимую аппаратную демонстрацию по этой фразе не утверждать.

Кандидат на вклад: количественная **совместная** feasibility map с uncertainty, loading/restart, physical readout и total energy/footprint; возможен model-based no-go result. Утверждать отсутствие аналогов или “first” пока нельзя. Структура публикации и итоговая comparative novelty table не завершены.

### 12.6. Историческая запись о статусе на момент прерывания

На момент создания checkpoint эти элементы оставались незавершёнными. Они завершены в разделах 13–18 ниже, в COMPUTATIONAL_ROADMAP.md и LITERATURE_AUDIT_RU.md. Полный проект представлен пользователю для рассмотрения; numerical implementation ожидает явного одобрения плана.

Подробное состояние, источники и порядок возобновления: [SESSION_HANDOFF.md](SESSION_HANDOFF.md). Исходные поисковые результаты и прочитанные excerpts сохранены отдельно в LITERATURE_SESSION_NOTES.json; они содержат непроверенные leads и не заменяют curated bibliography.

## 13. План вычислительного исследования и milestones

Этапы выполняются последовательно. Pass одной модели не означает pass полной системы; у каждой стадии есть вариант остановки и публикации отрицательного результата.

| Этап | Модель / вопрос | Минимальный объём | Go condition | No-go |
|---|---|---|---|---|
| M0: specification | Operating envelope и источники параметров определены? | species, T, pressure, geometry, noise priors, package/process bounds, source registry | Каждый вход имеет источник или явный uncertainty range; unknown не маркируется feasible | Нет ограниченного диапазона критического параметра: только conditional feasibility |
| M1: analytical | Существует ли окно q, RF frequency, voltage, depth, field, axial mode? | независимые формулы/units, Floquet, species comparison | Equations и units validated; остаётся хотя бы один physical scenario | Обнаружено противоречие уже в известных физических пределах; FEM не запускается для его маскировки |
| M1.5: early observability/entropy/system screen | Может ли ion-specific observable отделиться от noise, а stochastic/cooling/lifetime/power budgets совместиться? | image-current upper screen, FDT/noise and collision inequalities, explicit ranges and controls | Хотя бы один conditional interval остаётся, unknowns не переименованы в passes | Во всех обоснованных диапазонах сигнал ниже измеримого уровня или heating/loss невозможны для выбранных требований; missing essential data означает “unresolved”, не pass |
| M2: 2D fields | Numerical Laplace solution воспроизводит идеальную задачу? | один geometry, exact reference boundary, 3 mesh levels | M1/M1.5 оставили выбранный bounded scenario; сходятся field/Hessian, ideal coefficients | Нет convergence после geometry/exterior-boundary audit либо M1.5 отозвал физическую мотивацию выбранной точки |
| M3: 3D fields | Finite electrodes и DC ends сохраняют полезный quadrupole и axial stability? | baseline + symmetry control, coarse/medium/fine selected points, без sweep | конечная глубина, modes, multipoles и field peaks с mesh uncertainty | DC defocusing, end escape или multipole errors закрывают M1 window |
| M4: materials/circuit/heat | Doped Si безопасно доставляет RF при допустимой потере/нагреве? | complex conductivity, contacts, dielectrics, thermal solve; один higher-level EM check если нужен | EQS/EM разница в decision observables ≤5%; margins robust | conductivity/emission/power limits нарушены либо параметры не ограничены |
| M5: dynamics/lifetime | Удержание сохраняется при collisions, stray fields, RF switching? | deterministic reference first; затем validated kinetic collision MC или FDT-consistent Langevin | все mode stability; timestep convergence; survival CI meets target | RF heating/loss closes window, либо lifetime не обосновать |
| M6: readout/entropy | Ион добавляет measurable conditional entropy beyond circuit noise? | weighting fields, multiport circuit, front end, ion/no-ion controls, conservative entropy bound | lower confidence entropy rate ≥ useful target; source priors ограничены | apparent entropy от seed, drift, amplifier, unmeasured assumptions |
| M7: loading/module | Можно ли запускать, verify-one-ion, reload и упаковать узел? | capture, heating/ionization, getter/leak/outgas, state machine, energy and footprint | startup, availability, energy, area/volume pass M0 scenario | Нужен внешний macroscopic source/pump вне принятого integrated boundary |
| M8: feasibility map | Есть ли robust intersection и какой предел активен первым? | после валидации модулей — adaptive/multifidelity cases | Robust nonempty feasible set | Пустая область публикуется как отрицательный результат; критерии не ослабляются молча |

Numerical gates: для ideal 2D target refined-to-refined field/Hessian change ≤0.1%; для candidate finite geometry quadrupole residual ≤1% в \(\rho≤0.5r_0\); полевые величины должны сходиться на физических edge radii. Pseudopotential/RF frequencies должны совпадать в пределах 2% только для declared small-q range. EQS/full-EM agreement ≤5% для decision observables. Это proposed numerical tolerances, а не новые результаты или физические законы.

Не запускать FEM до прохождения M1 и первичного M1.5 review. Это исследовательский приоритет, не запрет публиковать analytical no-go. M2 — один bounded 2D benchmark, не большой sweep. M3 — выбранные 3D points. Расширение вычислений требует отдельного рассмотрения ресурсов; remote work не входит в план автоматически.

### 13.1. Предлагаемый benchmark полного устройства

У исходного запроса нет целевых throughput, power, area или startup. Для сравнительного исследования можно сохранить сценарий: 300 K reference с sensitivity до 350 K; useful output ≥100 kbit/s; \(P_{\rm total}≤1\) mW; die ≤1 mm²; package ≤10 mm³; startup/survival/availability targets по явному workload. Дополнительно rate/power combinations (10 kbit/s, 100 kbit/s, 1 Mbit/s; 0.1, 1, 10 mW) переводятся в energy/bit без объявления preferred point. Они нужны для P/E/C сравнения, а не для физического go/no. Каждая output-bit оценка должна указывать entropy accounting и extractor assumptions.

Это **сценарные исследовательские ориентиры**, не требования пользователя, не универсальные пороги и не measurements. Пример 1 mW / 100 kbit/s означает 10 nJ/useful bit — примерно \(8.6×10^4\) раз выше опубликованных 0.116 pJ/bit CMOS baseline для указанной operating point. Это критично для уровня C, но само по себе не опровергает уровень P.

При Poisson survival \(S(1h)=0.99\) соответствует \(\lambda_{\rm loss}≤2.79×10^{-6}\) s\(^{-1}\). Односекундная замена вносит downtime порядка \(\lambda t_{\rm recovery}\) при редких независимых потерях. Collision/technical hazard может быть non-Poissonian; формула тогда только initial screening.

### 13.2. Physical entropy gate

Entropy pass требует:

1. Один и тот же observation model для ion-present/no-ion и bath-only controls.
2. Демонстрации в модели, что процесс зависит от physical stochastic driver после conditioning на доступную боковую информацию.
3. Bounded source PSD/cross sections; если диапазон неизвестен, вывод явно conditional.
4. Связки noise, damping/cooling, RF heating, ion loss и source power.
5. Conditional min-entropy **rate**, достаточной для useful output; SNR и NIST pass сами по себе не проходят gate.
6. Entropy pipeline, который считает loss/health-test rejection/availability и не считает seeded simulation PRNG реальным источником.

Zero observed loss не доказывает long lifetime: для пуассоновской модели upper 95% rate при суммарной экспозиции \(T\) примерно \(3/T\). Метод rare-event/first-passage для предлагаемой RF-модели требует отдельной валидации.

20 \(k_BT\), ранее использованные для screening, не следуют из универсального lifetime law и **не являются hard gate**. Для физического go используются рассчитанные survival estimates. \(10^7\) V/m — иллюстративная field value, не гарантированный предел silicon. До получения применимого эмиссионного/material range это остаётся параметром, а не go/no boundary.

### 13.3. M1.5: ранний аналитический screen observability–entropy–cooling–power

До 2D FEM проверяем только необходимые условия и чувствительность, не заявляя entropy yield устройства. Начальная reproducible point: $^{40}\mathrm{Ca}^+$, $r_0=1\,\mu\mathrm m$, $V_{\mathrm{RF}}=10$ V differential peak, $q=0.2$, $T=300$ K, axial $f_z=1$ MHz. Auxiliary ranges в скриптах маркируются как sensitivity assumptions, пока не получены конкретные process/device data.

**Readout upper screen.** Для image current $i=Q\,\mathbf v\cdot\mathbf E_w$, где weighting field задаётся terminal unit-voltage problem. До расчёта field используем только масштаб $E_w\sim1/d_{\mathrm{eff}}$:
\[
I_{pk}\sim |Q|\omega_r A/d_{eff},\qquad
i_{n,rms}=\sqrt{4k_BT B/R_{in}}.
\]
При $A=100$ nm, $d_{\mathrm{eff}}=1\,\mu\mathrm m$, $f_r\simeq174.8$ MHz получаем $I_{pk}\sim17.6$ pA. Для условных $R_{in}=1$ MΩ, ENBW $B=1$ kHz и 300 K тепловой Norton-current noise равен 4.07 pA rms; номинальное amplitude/noise ratio ≈4.3. Диапазоны sensitivity: $A=10$–500 nm, $d_{\mathrm{eff}}=0.5$–5 µm, $R=10$ kΩ–1 MΩ, $B=1$ kHz–1 MHz; совместное сочетание экстремумов не является design. Не учтены terminal cancellation, cable/feedthrough, amplifier noise, input capacitance, resonant matching, ADC, aliasing, saturation и weighting-field geometry. $B$ — noise-equivalent bandwidth вокруг несущей, а не полоса baseband вместо 175 MHz. Нужен высокочастотный front end/demodulation с измеренной transfer function. Для полной симметрии четырёх RF-пар linear pickup может сокращаться. Следующий результат — FEM weighting field + multiport circuit, не coordinate readout.

**Entropy ceiling.** При $B=1$ kHz идеальный band-limited real process имеет порядок не более $2B=2$ k independent real degrees/s до учёта memory, ADC и rejected samples. Для однобитового event/no-event output это не более 2 kbit/s при таком ENBW. Чтобы один бинарный символ на sample потенциально поддержал 100 kbit/s, требуется как минимум $B=50$ kHz. Тогда в выбранном hypothetical $R_{in}=1$ MΩ Johnson noise возрастает до 28.8 pA rms и nominal 17.6 pA image-current scale / noise падает до ≈0.61 (без amplifier/cancellation/feedthrough). Для sparse independent Poisson indicator максимум достигается при $\Gamma\Delta t=\ln2$: в этом bandwidth scenario $\Gamma=6.93\times10^4$ s$^{-1}$ и $p\simeq0.096$–0.96 Pa для заданного k-range. Multi-bit analog samples могут нести больше одного бита, но entropy зависит от \(\sigma/\mathrm{LSB}\), tails, memory, physical source and side information. Частота инноваций сама не задаёт bits/sample. Conditional min-entropy остаётся unknown unless measured/model-supported; без физической distribution и temporal correlations entropy-rate claim отсутствует.

**Gas-collision branch.** Для идеального газа $\Gamma=nk=p k/(k_BT_g)$. При сценарном Langevin-scale $k=0.3$–$3\times10^{-15}$ m³/s (не измеренное значение для выбранной пары) и $T_g=300$ K при $p=10^{-8}$ Pa получаем $\Gamma\simeq7.2\times10^{-4}$–$7.2\times10^{-3}$ s$^{-1}$. Для независимого sparse Poisson event/no-event process бинарный min-entropy ceiling равен $R_{H,\min}=\Gamma/\ln2\simeq10^{-3}$–$10^{-2}$ bit/s. Event identity/timing может нести больше информации, но нужны timing resolution, event distinguishability, conditional correlations и ion-state reset. Для illustrative 100 kbit/s только event-indicator ceiling требует $\Gamma\ge6.9\times10^4$ s$^{-1}$, $p\simeq0.096$–0.96 Pa и ENBW не менее 50 kHz, если выход — бинарная серия event/no-event. Это не прогноз и не достаточное условие: нужен collision–RF energy exchange, species-specific charge exchange, cooling и survival analysis. При $S(1\,\mathrm h)=0.99$ требуется $\lambda_{loss}\le2.79\times10^{-6}$ s$^{-1}$; при целевой event rate доля теряющих столкновений должна быть <4.1×10$^{-11}$ при одинаковой статистике событий. Это жёсткое условие модели, не доказательство неизбежной потери. RF collision heating зависит от mass ratio и micromotion phase ([Viehland & Goeringer](https://doi.org/10.1063/1.1691405); [DeVoe](https://doi.org/10.1103/PhysRevLett.102.063001)).

**Circuit / surface-field noise and cooling.** Для weak one-sided $S_E(f_r)$, где $f$ задана в Hz, применяется
\[
\dot{\bar n}=\frac{Q^2S_E(f_r)}{4m\hbar\omega_r},\qquad
P_{ion}=\hbar\omega_r\dot{\bar n}=\frac{Q^2S_E(f_r)}{4m}.
\]
Для Ca$^+$ при $f_r\approx174.8$ MHz и $D_{\eta=0.5}\simeq0.0625$ eV formal sensitivity over $S_E=10^{-15}$–$10^{-6}$ V²m⁻²Hz⁻¹ даёт $\dot n\sim8.3\times10^{-4}$–$8.3\times10^5$ quanta/s и no-cooling energy-to-acceptance scale $D_\eta/P_{ion}\sim10^8$–$10^{-1}$ s. Это **подстановка noise range, измеренного в других ловушках**, не прогноз: Brownnutt et al. сводят значения примерно $10^{-15}$–$10^{-6}$ преимущественно при 0.1–20 MHz и 30–3500 µm; здесь ~175 MHz и micron/submicron Si/oxide. Extrapolation неизвестна, а $D_\eta$ не escape saddle. Такая 9-decade sensitivity не позволяет признать surface noise установленным entropy source. Circuit Johnson contribution следует считать из actual $\Re Z(\omega,T)$ и transfer field, не приравнивать к surface PSD. Stationary ion требует cooling/removal не ниже суммы absorbed motional power с учётом RF/collision heating; cooling может подавить полезную флуктуацию и расходует laser/resistive/feedback power ([Brownnutt et al.](https://doi.org/10.1103/RevModPhys.87.1419); [Turchette et al.](https://doi.org/10.1103/PhysRevA.61.063418)).

Для white Johnson noise correlation time определяется полосой и network poles; surface charge noise обычно colored/telegraph, и её correlation time нельзя вывести из одной PSD value. M1.5 не подменяет неизвестный spectrum независимыми Gaussian samples. Energy budget: $P_{avg}=P_{steady}+E_{start}/T_{restart}$ и $E_{useful\,bit}=P_{avg}/R_{useful}$. Например, 1 mW/100 kbit/s = 10 nJ/bit, C-level comparison scenario. Полный $P_{steady}$ неизвестен: RF source+resonator, ionization/oven, cavity/getter, cooling/laser, readout+ADC и controller/health tests. Для условного startup 1 mJ раз в час средняя добавка 0.278 µW; нужно подтвердить, что такие energy и interval относятся к реальной single-ion loading.

Сводное сравнение трёх физических кандидатов:

| Механизм | Coupling и ожидаемый signal | Heating/cooling и lifetime | Correlation / entropy-rate inference | Integration |
|---|---|---|---|---|
| Johnson–Nyquist | Circuit thermal voltage/current noise; ион чувствует только через transfer impedance и electrode field. Номинальный circuit current 4.07 pA rms в выбранном $1$ MΩ/$1$ kHz screen; ion-specific 17.6 pA scale задаётся assumed amplitude/weighting, не выводится из Johnson source. | FDT связывает drive и dissipation; для stationary motion нельзя добавлять noise без matching damping. Cooling/thermalization от сети может уменьшать SNR; heating по $S_E$; loss через excitation при недостаточной глубине. | Correlation time ограничена network bandwidth/poles, если noise model known. Direct electronic sampling может быть дешевле и не хуже по entropy; trap adds no entropy by data processing alone. | Низкая сложность как electronic baseline; high-frequency electrode impedance, parasitic C, phase noise and pair cancellation still require circuit design. |
| Buffer-gas collisions | Discrete momentum transfers; image-current bursts/changes in secular amplitude possible. At assumed $k$ range and $p=10^{-8}$ Pa: $\Gamma=7.2\times10^{-4}$–$7.2\times10^{-3}$ s$^{-1}$; ideal event-indicator ceiling $10^{-3}$–$10^{-2}$ bit/s. | Each collision can cool or heat depending on mass ratio, RF phase and differential cross-section; nonthermal tails/charge exchange can eject the ion. Targeting 100 kbit/s by event indicator implies $p\sim0.1$–1 Pa and <4.1×10$^{-11}$ loss probability/event for 99% hourly survival in the simple Poisson model. | Collision interval is not necessarily innovation time of observed motion; secular state retains memory. Poisson bound is an optimistic model ceiling, not a min-entropy claim for a real trap. | Gas reservoir, pressure control/seal, purity and collision/loss calibration; moderate pressure simplifies vacuum but increases dynamical and contamination burden. |
| Silicon surface-charge / field fluctuations | Direct stochastic force $Q E_s(t)$; secular image current filtered through weighting field. No reliable project-specific amplitude is available. | Same fluctuations heat motion via $\dot n=Q^2S_E/(4m\hbar\omega)$; if strong enough to observe they may require commensurate cooling and reduce lifetime. | PSD amplitude alone gives neither correlation time nor conditional Hmin rate; patch/TLS spectra may be colored and nonstationary. Local GHz/nm–µm extrapolation unresolved. | No added source hardware, but surface preparation, oxide/interface control, drift, radiation/charge management and process reproducibility are difficult. A CMOS surface-noise readout is a required control. |

**Decision:** не назначать primary entropy mechanism. Johnson–Nyquist — primary *control/model-validation bath*, gas collisions — conditional candidate для ion-mediated event entropy, а surface charge — high-uncertainty candidate. Продолжение возможно только при source-specific measurement priors or defensible bounds. У trapped ion есть data-processing constraint: его readout не может иметь больше conditional entropy, чем доступно в physical driver plus independent measurement noise; последние следует разнести controls.

### 13.4. Три уровня критериев и классификация каждой найденной границы

| Уровень | Условие для pass в заявленном workload | Что не считается pass |
|---|---|---|
| P — Physical | Floquet/full-RF stable volume + model-supported stochastic driver + observable ion-mediated contribution + defensible conditional entropy/lifetime interval | Одна stable trajectory, SNR без weighting model, seeded RNG/NIST pass |
| E — Engineering | Integrated source, verified single-ion initialization/recovery, seal/pressure, process flow, thermal/power/reliability budgets in chosen package boundary | Silicon electrodes alone или trap-only power |
| C — Competitive | Common-boundary comparison with CMOS TRNG: useful entropy rate, joule/bit, area, uptime, health/security overhead and lifetime | Любые универсальные пороги 100 kbit/s, 1 mW, 1 mm², 10 mm³ |

Для каждого reported size boundary в таблицах/фигурах указывать метку: F (fundamental under stated assumptions), T (technology), R (requirement), U (model/data unresolved). Например, Mathieu stability curve — F для ideal equation; field-emission onset в doped Si — T; acceptable ion survival/output latency — R; nanometer-scale $S_E$, charge-exchange probability и sealed-cavity leakage — U. У границы приводить interval/sensitivity and source. Если U перекрывает feasible region, вывод остаётся conditional; нельзя рисовать однозначную минимальную ловушку.

## 14. Предлагаемая структура Python проекта

Новая модель должна быть отделена от существующих прямоугольных trap и semiconductor impedance расчётов:

~~~text
research/ion_trap_rng/
  RESEARCH_PLAN_RU.md
  RESEARCH_SUMMARY_RU.md
  COMPUTATIONAL_ROADMAP.md
  SESSION_HANDOFF.md
  LITERATURE_AUDIT_RU.md
  pyproject.toml
  README.md
  configs/ species.yaml, baseline.yaml, uncertainty.yaml
  data/ source_registry.yaml, material_models/
  src/iontrap_rng/
    constants.py
    geometry/ hyperbola.py, finite_electrodes.py
    analytic/ mathieu.py, scaling.py, pseudopotential.py, collisions.py
    fields/ electrostatics.py, eqs.py, weighting.py, multipoles.py
    dynamics/ floquet.py, rf_integrator.py, collision_mc.py, langevin.py
    circuit/ ports.py, readout.py, thermal.py
    entropy/ estimators.py, extractor_budget.py, health_tests.py
    system/ loading.py, vacuum.py, power_area.py, feasibility.py
    io/ manifests.py, units.py
  tests/ unit/, analytic_references/, geometry/, convergence/, stochastic_controls/
  scripts/ run_reference_case.py, run_selected_convergence.py, make_feasibility_map.py
  results/ manifests/, selected_cases/, figures/
~~~

SI units внутренне; входные/выходные values имеют explicit units. YAML configurations immutable per run and content-hashed. Каждый manifest сохраняет source revision, solver, geometry, mesh, random seed, materials и uncertainty case. Большие mesh outputs в Git не сохраняются.

Candidate stack: Python 3.12, NumPy/SciPy, pytest, Gmsh и scikit-fem для переносимой начальной 2D/selected-3D стадии. Это кандидаты, а не выбранные и проверенные зависимости. Перед установкой нужно проверить версии, licence, Windows support и solver capabilities. Full-Maxwell production solver не следует писать с нуля для этого проекта. Существующий root .venv/requirements относится к другим моделям проекта; менять его не требуется.

## 15. Validation и convergence protocol

### 15.1. Analytical / units

- Утвердить differential peak voltage convention; проверять знаки electrode pair и charge reversal.
- Проверить \(\nabla^2\Phi=0\), идеальный E field, trace-free DC Hessian, четыре branch equations, curvature и normal orientation.
- Независимо сверить Mathieu characteristic exponent/Floquet multipliers с первой stability region; recover small-q \(\beta\simeq q/\sqrt2\).
- Проверить коэффициент \(1/4\) pseudopotential для **peak** field сравнением с cycle-average forced micromotion energy; запрещено смешивать peak/RMS.
- Log–log unit checks scaling laws; units image energy, collision rate, PSD, damping and diffusion.

### 15.2. Geometry / fields

- Topology: four isolated conductive branches, central opening, clearances, normals, supports and electrical nets.
- Разделить manufactured harmonic-polynomial solve (численная ошибка) и physical hyperbola equipotential problem (геометрическое и физическое отклонение).
- Full-domain 2D до symmetry reduction; затем сравнить symmetry/full models.
- Независимо refine mesh, curved boundaries, external domain, edge radius, fit radius \(\eta=0.2,0.3,0.5\) и polynomial order.
- Track field residual, null displacement, Hessian, multipoles, surface fields, flux and balance — относительная ошибка в нулевом field point не используется.
- 3D limit \(L/r_0\to\infty\) сравнить с центральным 2D reference; торцевое DC field — с аналитической квадрупольной формой.
- Углы всегда имеют finite radius; sharp-corner FEM peak не равен реальному breakdown.

### 15.3. EM/circuit/thermal/dynamics

- EQS проверить на homogeneous capacitor, analytic slab и RC line; reaction-current, surface flux, stored energy, Joule loss и source power согласовать.
- Проверить \(\Omega\epsilon/\sigma\), skin depth, contacts и equipotential-face approximation; выбранный худший режим сопоставить с higher-level EM, если критерии EQS нарушены.
- Heat solver: аналитическое slab resistance, power conservation, thermal grid refinement и температурная зависимость conductivity.
- Dynamics: harmonic oscillator, Mathieu reference, RF phase, независимые integrators; timestep \(T_{\rm RF}/16,32,64\) на нескольких выбранных точках.
- В static/no damping проверять energy conservation; при RF вести work–energy balance \(\Delta E=\int Q\mathbf E\cdot\mathbf v\,dt\).
- Langevin: thermal covariance/equipartition in harmonic well, correct diffusion/FDT, timestep convergence и colored-noise PSD.
- Collision MC: Maxwell sampling, energy-dependent scattering, rate \(n\langle\sigma v\rangle\), momentum conservation.
- Absorbing surfaces/ends, saddle, chemical loss; survival confidence intervals и независимый first-passage benchmark на synthetic case.

### 15.4. Readout/entropy/system

- Weighting field сверить с dipole/plate analytic tests; проверить reciprocity, terminal signs, current conservation, RF-pair cancellation and imbalance.
- Electrical transfer/noise восстанавливает known RLC/Johnson references; тестировать alias, saturation, ADC resolution, phase-locked feedthrough.
- Synthetic entropy controls: biased IID, Markov, periodic, colored Gaussian, deterministic-plus-noise; known PRNG seed обозначить simulation convenience.
- Сравнить ion present/no ion, source on/off, bath only и electronics-only; не считать amplifier noise энтропией движения иона.
- Restart/false-ready/multiple-ion errors считать отдельно. Условные entropy estimates и extractor budget перепроверить на controls с известным ground truth.
- Energy boundary: wall power уже включает соответствующие losses; не суммировать их повторно. Startup/reload energy отдельно, затем усреднить по workload.

Ни один большой sweep до узких analytic references. Сначала selected points; затем, только при успехе, adaptive parameter sampling. Сходимость подтверждает математическое решение заявленной модели, но не измеряет неизвестный material/noise input.

## 16. Литература, novelty и пределы вывода

| Работа | Verified DOI | Что именно установлено | Не устанавливает |
|---|---|---|---|
| Paul (1990), *Electromagnetic traps for charged and neutral particles* | [10.1103/RevModPhys.62.531](https://doi.org/10.1103/RevModPhys.62.531) | Основы RF confinement | Наноловушку |
| Leibfried et al. (2003), *Quantum dynamics of single trapped ions* | [10.1103/RevModPhys.75.281](https://doi.org/10.1103/RevModPhys.75.281) | Mathieu/Floquet, quantum и motional dynamics | RNG устройство |
| Taylor et al. (2001), silicon MEMS quadrupole mass filter | [10.1116/1.1359172](https://doi.org/10.1116/1.1359172) | Silicon mass-filter technology precedent | Четыре silicon hyperbolas/ion trap |
| Pau et al. (2006), 20-µm cylindrical doped-polysilicon trap array | [10.1103/PhysRevLett.96.120801](https://doi.org/10.1103/PhysRevLett.96.120801) | Silicon microfabricated mass-analysis trap | Atomic-ion motional RNG |
| Britton et al. (2009), *Scalable arrays of rf Paul traps in degenerate Si* | [10.1063/1.3254188](https://doi.org/10.1063/1.3254188) | Doped conductive silicon RF Paul traps and measured heating | Nanometer scaling |
| Mehta et al. (2014), *Ion traps fabricated in a CMOS foundry* | [10.1063/1.4892061](https://doi.org/10.1063/1.4892061) | Surface trap in 90-nm CMOS process | 90-nm electrode distance |
| Turchette et al. (2000), ion motional heating from electric noise | [10.1103/PhysRevA.61.063418](https://doi.org/10.1103/PhysRevA.61.063418) | Field noise heats trapped motion | Noise spectrum of proposed nanosurface |
| Sedlacek et al. (2018), surface-trap noise distance scaling | [10.1103/PhysRevA.97.020302](https://doi.org/10.1103/PhysRevA.97.020302) | Exponent near 4 in a particular device | Universal nanometer law |
| Tian et al. (2014), ion-trap miniaturization scaling | [10.1002/jms.3343](https://doi.org/10.1002/jms.3343) | Voltage/frequency/power/depth constraints for cylindrical traps | Silicon single-ion RNG |
| Viehland & Goeringer (2004), RFQ buffer-gas kinetics | [10.1063/1.1691405](https://doi.org/10.1063/1.1691405) | Collision-derived damping/stochastic forces | Our species’ rates without cross sections |
| DeVoe (2009), RF trap ion–buffer-gas power-law distributions | [10.1103/PhysRevLett.102.063001](https://doi.org/10.1103/PhysRevLett.102.063001) | RF collisions can generate non-Maxwellian tails | Universal effective temperature |
| Li et al. (online 2025), hyperbolic LIT mass-spectrometer simulation | [10.1021/jasms.5c00307](https://doi.org/10.1021/jasms.5c00307) | Hyperbolic geometry and multipole/ejection design | Fabricated silicon RNG |
| Nguyen et al. (2025), trap geometry scaling | [10.3390/e27060576](https://doi.org/10.3390/e27060576) | Geometry affects height, harmonicity, depth, frequency | Complete module feasibility |
| Sterk et al. (2026), low-RF-power ion-trap designs | [10.1116/5.0332336](https://doi.org/10.1116/5.0332336) | Fabricated loss-reduction approaches | General nanoscale limit |
| Pironio et al. (2010), Bell-certified trapped-ion randomness | [10.1038/nature09008](https://doi.org/10.1038/nature09008) | Quantum-certified randomness protocol | Classical motional randomness |
| Um et al. (2013), contextuality RNG with a single trapped ion | [10.1038/srep01627](https://doi.org/10.1038/srep01627) | Qutrit measurement and min-entropy bound | Bare deterministic secular motion |
| Kumar et al. (2025), integrated micro-oven | [10.1063/5.0271048](https://doi.org/10.1063/5.0271048) | 90×90 µm² heater, tens-of-ms neutral fluorescence response, mW scale | Same single-ion capture time; lifetime claim was calculated |
| Schwindt et al. (2016), miniaturized vacuum package | [10.1063/1.4948739](https://doi.org/10.1063/1.4948739) | ~1 cm³ trapped-ion package precedent | mm³ package |
| Thrasher & Schwindt (2025), 3-cm³ sealed package | [10.1063/5.0243517](https://doi.org/10.1063/5.0243517) | Ion lifetime up to 50 days in reported sealed packages | Lifetime for proposed trap |
| Ramo (1939), induced current theorem | [10.1109/JRPROC.1939.228757](https://doi.org/10.1109/JRPROC.1939.228757) | Weighting-field approach to induced current | Adequate SNR in this design |
| Setzer et al. (2021), monolithic SPAD ion-fluorescence detector | [10.1063/5.0055999](https://doi.org/10.1063/5.0055999) | Integrated optical ion detection | Required continuous entropy throughput |
| NIST SP 800-90B (2018) / SP 800-22 Rev.1a (2010) | [10.6028/NIST.SP.800-90B](https://doi.org/10.6028/NIST.SP.800-90B) / [10.6028/NIST.SP.800-22r1a](https://doi.org/10.6028/NIST.SP.800-22r1a) | Entropy-source/health and statistical-test frameworks | Security proof or true-source model |
| Wang et al. (2024), 130-nm CMOS latch TRNG | [10.1109/TVLSI.2023.3328602](https://doi.org/10.1109/TVLSI.2023.3328602) | Reported 5343 µm² and 0.116 pJ/bit incl. von Neumann processing | Independently audited common-boundary comparison |

Records were checked against publisher, official NIST, or institutional/author records within the available source audit. A fuller earlier audit is [research/LITERATURE_AUDIT.md](../LITERATURE_AUDIT.md). Raw search results remain in [LITERATURE_SESSION_NOTES.json](LITERATURE_SESSION_NOTES.json); that archive has duplicate/unverified leads and is not a bibliography.

The focused literature scan used 23 queries and returned 120 URLs, of which 116 were unique. It is not exhaustive. Silicon microtraps, degenerate-Si Paul traps, CMOS surface traps, ion quantum RNG, micro-ovens, hermetic packages, scaling studies and hyperbolic mass-spectrometer simulations are all precedents. No source in this bounded set combines four separate silicon hyperbola-branch electrodes, finite geometry, axial DC, loading/cavity, physically grounded readout, entropy–cooling–lifetime tradeoffs and complete module budgets. This absence is not proof of priority. Repeat a citation-chain review before publication. Defensible novelty candidate: a robust uncertainty-aware joint feasibility map or a negative result naming the first constraint that closes it. Avoid “first-ever” and component-level novelty claims.

## 17. Итоговая физическая проверка

- Differential peak RF convention is consistent in potentials, Mathieu q, frequency and \(D_\eta\).
- The four branches create diagonal apertures; \(2r_0\) is vertex spacing. Finite-domain FEM needs a specified exterior boundary.
- Axial confinement requires extra DC structures and necessarily radial defocusing.
- The infinite ideal pseudopotential has no escape saddle; \(D_\eta\) is explicitly an acceptance metric.
- Fixed-q scaling does not preserve energy in units of \(k_BT\), image forces, noise, circuit or readout.
- The thermal displacement table assumes equilibrium; RF noise/collisions may heat and eject the ion. Cooling/noise obey a joint balance.
- Chaos and seeded simulation are not physical entropy. Classical motional dynamics is deterministic absent a stochastic driver.
- Symmetric paired electrodes can cancel linear image-current pickup; weighting fields, complete readout circuit and ion/no-ion controls are necessary.
- Pressure and lifetime depend on species-specific collisions and leakage. The example \(k_{\rm loss}=10^{-18}\) m³/s is assumed, not measured.
- Silicon nanogap emission limit is unknown here; \(10^7\) V/m is illustrative only.
- CMOS process node is not a trap gap; source-heater area excludes package; cm³ package does not prove mm³ feasibility.
- NIST tests cannot certify cryptographic security or unmeasured model assumptions.

Главное unresolved physical question: найдётся ли совместимый уровень stochastic excitation и cooling, при котором условная entropy rate полезна, ion-specific signal измерим, а lifetime и power остаются приемлемыми? Это первичный falsifiable gate, до 3D optimization.

## 18. Статус и рекомендуемый старт

Analytical model M1 is implemented in src/iontrap_rng and validated by 12 unit tests; the baseline and M1.5 sensitivity results are in results/m1_reference.json, results/m15_screen.json and M1_ANALYTICAL_REPORT_RU.md. The axial DC electric-field sign issue found during review is corrected and regression-tested. The ideal hyperbolic trap is Floquet-stable at selected q=.2, but entropy rate, device-specific noise/lifetime, loading and full-module power remain unproved. M1.5 did not select a primary entropy mechanism. The bounded M2a 2D weighting-field/readout study is now complete and separately reported in M2A_READOUT_REPORT_RU.md; it refines the earlier scalar image-current estimate but does not establish receiver performance or RNG feasibility. The current package passes 17 tests. No 3D FEM, stochastic trajectories, large sweeps or remote computation have been run.

## 19. Stage M2a: ограниченный 2D readout benchmark

M2a solved two distinct problems. The ideal solver validation used a rectangular domain with exact analytical quadrupole boundary data. The finite device used four separate truncated hyperbolic strip electrodes with thickness 0.1*r0, hyperbolic extent 1.5*r0, rounded in-plane caps of radius 0.05*r0, and a grounded square boundary at 6*r0. The finite-device external boundary was not forced to the ideal quadrupole. The four conductors were separate terminals and weighting potentials were solved by one-terminal-at-1 V boundary conditions.

Local triangular P1 h refinement covered h_min/r0 = 1/20, 1/40, 1/80, 1/160, and 1/320; CAD boundary sampling was independently checked at 50, 100, and 200 points per branch. At the finest mesh the ideal-reference center-Hessian relative error was 0.173%, the normalized differential weighting gradient was 0.399741/r0, the differential-port capacitance was 3.74953e-11 F/m, and the per-electrode capacitance was 7.29524e-11 F/m. The 1/160-to-1/320 changes were ~0.060% in differential gradient and ~0.019% in per-electrode capacitance. Global relative electric-field L2 error remained ~1.73% because the far domain is coarsened. A local degree-four multipole fit remains more sensitive and is not treated as a precision coefficient determination.

Symmetrically paired opposite electrodes suppress the linear center weighting gradient to 1.76e-4/r0, about 4.39e-4 of the single-electrode gradient. A normalized differential port restores approximately 0.39974/r0. For a conditional 100 nm Ca+ amplitude at 174.784 MHz and r0=1 µm, that gives 7.03 pA peak. With a 1 pF external parasitic and the stated illustrative amplifier noise, an unmatched 1 MΩ parallel-RC receiver gives only ~0.10 signal-to-rms-noise ratio at 1 kHz ENBW; resonant scenarios improve only under restrictive Q/capacitance assumptions. The former 17.6 pA scalar-gradient estimate is therefore not a readout feasibility result.

The report explicitly distinguishes F/T/R/U limits. At 10 nm, the chosen 1 nm electrode backing makes continuum doping/conductivity and smooth surfaces questionable (T/U); FEM element size is numerical, not an atomistic physical resolution. No universal minimum trap size is inferred. Conventional RF-terminal pickup is a no-go under symmetry; differential sensing merits a conditional readout-modeling step only after parasitic capacitance, loaded Q, amplifier voltage/current noise, CMRR, RF feedthrough, and electrode length are specified. No physical entropy or useful RNG rate is established. Details and figures are in [M2A_READOUT_REPORT_RU.md](M2A_READOUT_REPORT_RU.md), [m2a_readout.json](results/m2a_readout.json), and the two [M2a figures](results/figures/m2a_potential_weighting.png).

**Historical M2a stop point:** at that time 3D FEM was not authorized. This restriction was superseded by the later explicit M2b authorization; see the addendum below.

## M2b — 3D validation preflight (27 October 2026)

### Approved scope completed

- Added `scripts/run_m2b_manufactured.py`, configuration `configs/m2b_validation.yaml`, and `tests/test_m2b_manufactured.py` for the exact harmonic reference Φ=(x²−y²)/2. Across three mesh sizes P1 field error converges first order and P2 reproduces the exact quadratic reference to floating-point/linear-solve precision. P1 analytical element-moment integration independently agrees with FE quadrature. P1's constant element Hessian has 100% RMS error for this case; P2's recovered element Hessian error is <4e−13.
- Because the focused literature audit did not support the requested exact silicon electrode/support/dielectric stack, the physical phase was deliberately a preflight rather than an asserted realistic process. The declared geometry in `configs/m2b_preflight.yaml` has four finite hyperbolic branches, L=10r0, two annular axial DC electrodes, a grounded enclosure, and empty vacuum. Silicon bodies are ideal equipotentials. SiO2, supports, contacts, doping-dependent conductivity and packaging dielectrics are absent.
- `scripts/run_m2b_3d_preflight.py` uses Gmsh tetrahedra and scikit-fem P1; three bounded meshes produce RF/DC potentials, six terminal weighting solutions, differential weighting gradient and an energy-derived Maxwell capacitance matrix. `scripts/plot_m2b_3d_preflight.py` makes the saved potential/field slices and geometry image.

### Validation and limits

Manufactured results: P1 global E-L2 = 8.333%, 4.167%, 2.083% for 12/24/48 divisions per side; local ρ<r0 E-L2 = 56.87%, 28.55%, 14.34%. P2 global errors are (2.5\cdot10^{-15}ldots2.9\cdot10^{-14}) and local errors (6.1\cdot10^{-15}ldots6.9\cdot10^{-14}). These tiny values validate polynomial representation and implementation on the manufactured square; they do not bound device geometry error.

In 3D the two finest mesh levels change differential weighting-gradient magnitude by about 0.21%, differential capacitance by 2.42%, and the Maxwell matrix norm by 3.4%. RF Hessian components change by 9–82%, DC Hessian components by 24–60%; RF Hessian trace residual remains about 10% of its norm. Thus 3D second derivatives/multipoles are not numerically qualified. The coarse center fit is under-resolved and not used for interpretation. The 3D model also does not calculate dielectric-boundary distortion, material loss, RF heating, parasitics, cavity modes or a port transfer function.

### M2b decision

The ideal geometry gives a plausible quadrupole-like center field, axial DC curvature and nonzero differential weighting field. At r0=1 µm the finest conditional weighting-gradient coefficient is 0.4135/r0 for the normalized unit-terminal weighting potential; for assumed Ca+ motion at 174.8 MHz with A=0.1r0 this corresponds to a Ramo current scale about 7.3 pA. This is not measured SNR, entropy or RNG throughput. Differential capacitance is about 0.54 fF relative to the modeled grounded enclosure; it is not total receiver input capacitance.

**Decision: no-go for ion-dynamics or entropy calculations using this 3D field as a validated device model.** First correct the central-fit diagnostics, demonstrate local h/p and curved-boundary convergence for Hessian/multipoles, and select/source-bound a real dielectric/support/terminal stack. The geometry remains an ideal electrostatic preflight only. See [M2b report](M2B_3D_FEASIBILITY_REPORT_RU.md), [JSON](results/m2b_3d_preflight.json), [manufactured validation](results/m2b_manufactured_validation.json), and [figures](results/figures/m2b_3d_potential_slices.png).

Boundary limits remain classified rather than universal: continuum Laplace physics is F within assumptions; fabrication of this geometry/stack is T and unsupported; required SNR/rate/power are R and not set here; dielectric effects, local material response and numerical Hessian error are U. No sub-100-nm material extrapolation was made.







