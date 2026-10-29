# Научное резюме

**Проект:** физические пределы миниатюризации четырёхэлектродной кремниевой RF quadrupole как основы интегрированного RNG. M1–M2b.2, короткий deterministic M3a и DC-force audit M3a.1 выполнены локально. M3a сохранил короткие траектории и локальные Floquet modes; M3a.1 выявил 4–13% межсеточное расхождение слабой DC силы и заблокировал длительную динамику. Тесты: 31 passed. Никакие axial-cycle trajectories, stochastic modelling, entropy или RNG не запускались в M3a.1.

## Гипотеза

Проверить, существует ли robust диапазон геометрии, RF/DC параметров, среды, шума, cooling, readout и ion loading, где один ион долго удерживается, даёт измеримый conditional min-entropy, регулярно перезапускается и помещается в полные budgets. Размеры электродов, trapping aperture и package считаются раздельно.

Основная geometry — четыре отдельные гиперболические ветви \(x^2-y^2=\pm r_0^2\), экструдированные вдоль z, с alternating RF на противоположных парах. Дополнительные DC electrodes создают axial confinement, не заменяя RF quadrupole.

## Что известно

Идеальный потенциал
\[
\Phi=\frac{V_{\rm RF}\cos\Omega t}{2r_0^2}(x^2-y^2)
\]
даёт linear Mathieu equations. При fixed q уменьшение \(r_0\) требует роста frequency как \(1/r_0\), если differential voltage постоянно. Оно не сохраняет energy barrier относительно \(k_BT\). Axial DC обязательно defocuses radial plane. Стабильная идеальная траектория одного иона детерминирована.

Есть precedents silicon microtraps, degenerate-silicon Paul traps, CMOS-foundry surface traps, trapped-ion quantum RNG, MEMS micro-ovens и sealed cm³ ion packages. Novelty может состоять в совместной количественной feasibility всей системы.

## Главные барьеры

- Низкая энергия удержания относительно thermal/collision tails при росте RF frequency, voltage и surface field.
- Electric-field noise даёт stochastic heating и может служить сигналом, но одновременно сокращает ion lifetime; nanometer surface spectrum неизвестен.
- Без cooling шум нагревает ион до потери; охлаждение требует power и добавляет backaction.
- Симметрично объединённые electrodes могут гасить линейный image-current readout.
- Источник, getter, cavity, readout и control могут превосходить ловушку по размеру и энергии.
- Chaos, seeded PRNG и NIST statistical pass сами по себе не доказывают физическую min-entropy.

Аналитическая оценка для \(^{40}\mathrm{Ca}^+\), q=0.2, differential peak 10 V даёт \(f_{\rm RF}\approx2.47\) GHz при \(r_0=1\) µm и 24.7 GHz при 100 nm. При \(\eta=0.5\), \(D_\eta\approx0.0625\) eV \(\approx2.42k_BT\) при 300 K. Это ideal screening, не finite trap depth и не FEM result. M1 Floquet result при q=0.2: \(\beta=0.1425513\), \(f_{sec}=176.18\) MHz без axial DC и \(\det M=0.99999999999985\). С добавлением \(f_z=1\) MHz small-q radial estimate — 174.78 MHz.

## Уровни feasibility

Физическая возможность P (confinement + измеряемый physical entropy), инженерная реализуемость E (автономная загрузка, cavity, recovery, процесс и полный power budget) и конкурентоспособность C (сравнение полной joule/bit/rate/area/reliability с CMOS TRNG) оцениваются отдельно. Сценарии 100 kbit/s и 1 mW не являются физическими go/no thresholds.

## M1.5: ранний observable/entropy screen

Для прежней scalar оценки \(A=100\) nm, \(d_{\mathrm{eff}}=1\) µm и \(^{40}\mathrm{Ca}^+\) получалось около 17.6 pA; M2a показал, что такой эффективный градиент нельзя принимать без weighting-field решения. В выбранной конечной гиперболической геометрии нормированный differential gradient равен \(0.39974/r_0\), или около 7.03 pA peak для тех же 100 nm на 174.8 MHz. Обычная симметричная пара подавляет линейный центр-сигнал примерно в 2 300 раз относительно differential pickup. При 1 pF, 1 MΩ, ENBW 1 kHz и заданных иллюстративных шумовых плотностях усилителя signal/rms-noise около 0.10; это receiver-сценарий, не измеренная производительность. Conditional min-entropy rate остаётся неизвестной.

Подстановка literature sensitivity range $S_E=10^{-15}$–$10^{-6}$ V²m⁻²Hz⁻¹ в weak-noise equation даёт для Ca+ формальный heating $8\times10^{-4}$–$8\times10^5$ quanta/s. Эти данные относятся к другим расстояниям, частотам и поверхностям; при micron Si на 175 MHz это только демонстрация девяти порядков неопределённости. При условном gas collision coefficient $k=0.3$–$3\times10^{-15}$ m³/s и $p=10^{-8}$ Pa идеальный independent-Poisson event-indicator ceiling составляет $0.001$–$0.01$ bit/s. Для идеального 100 kbit/s binary-event ceiling требуется $p\approx0.096$–$0.96$ Pa и ENBW не менее 50 kHz; при M1.5 ENBW 1 kHz бинарный выход ограничен 2 kbit/s. На минимальных 50 kHz assumed 1 MΩ Johnson noise уже 28.8 pA rms против nominal 17.6 pA image-current scale. Простой 99%/hour survival budget тогда требует loss probability $<4.1\times10^{-11}$ на collision. Это условный расчёт, до оценки RF heating/cooling.

Johnson–Nyquist — наиболее удобный baseline/control, но может давать entropy дешевле напрямую в CMOS circuit. Gas collisions — conditional candidate for ion-mediated entropy, ограниченный RF heating/lifetime. Silicon surface fluctuations — физически прямое coupling, но local spectrum/корреляции не ограничены данными; ни один источник пока нельзя обоснованно выбрать как primary RNG mechanism.

## Приоритет моделирования

1. Аналитические/Floquet reference equations и dimensional checks.
2. Один 2D Laplace benchmark четырёх идеальных hyperbolas с exact outer boundary и mesh convergence.
3. Затем finite 3D geometry, axial DC, физические multipoles, finite escape saddle.
4. До широких sweeps отдельно проверить ion-specific readout, noise–cooling–lifetime и loading/package budgets.

## Возможный вклад

Uncertainty-aware feasibility map, которая связывает ideal hyperbola fields с finite silicon, ion lifetime, conditional entropy, measurable readout, loading/restart и full-module power/footprint. Пустое robust-feasible пересечение или первый constraint, который закрывает окно, — ценный отрицательный результат. “First” claim не обоснован без systematic review.

## Историческое резюме после M2b (пересмотрено после M3a)

M1/M1.5 и M2a/M2a.1 прошли ограниченную аналитику и 2D readout screen. M2b P2 manufactured test подтвердил точное представление квадрупольной функции и ожидаемую P1 сходимость. Один 3D vacuum/perfect-conductor preflight построен с четырьмя конечными hyperbolic branches и кольцевыми DC endcaps. Получены differential weighting field порядка 0.4135/r0, differential capacitance около 0.54 fF и шесть terminal basis solutions. RF/DC Hessian не сошёлся: различия между сетками остаются велики, RF trace residual около 10%. Эти результаты пригодны как диагностический preflight, но не для предсказания ионной динамики.

Фокусная литература не даёт полностью совпадающей реализации с объёмными четырьмя гиперболическими Si электродами, конкретным dielectric/support stack и differential terminal routing. Поэтому M2b не является CMOS process model. Условные 7.3 pA при (^{40}\mathrm{Ca}^+\), 174.8 MHz и амплитуде движения 0.1r0 — только image-current scale для заданной траектории. SNR, физическая энтропия, entropy rate, охлаждение, полный энергобаланс, ion loading/restart и competitiveness против CMOS TRNG не установлены.

**Рекомендуемый следующий шаг:** не начинать stochastic trajectories. Сначала улучшить 3D local-fit diagnostic и провести локальную P1/P2/h-refinement валидацию одного электростатического CAD, включая curved boundaries, terminal flux/energy reciprocity и Hessian/multipole estimates; параллельно определить источник-поддержанный dielectric/support/terminal stack. Только после этого переходить к dynamics. Энтропию оценивать при наличии физической модели stochastic source, cooling и receiver transfer/noise.

## Актуальное резюме M3a — 28 октября 2026

## Актуальное резюме M3a.1 — 28 октября 2026

При фиксированном `VDC=0.08697882845 V` на каждом DC кольце три P2 сетки дают локальные осевые частоты 997.220–1,000.001 kHz, однако сами restoring forces на смещениях `0.1–0.2r0` расходятся на 3.96–12.66%; fine-grid осевая сила при `±0.1r0` асимметрична примерно на 13.7%. Target 1 MHz ранее был задан voltage calibration по fine-grid Hessian, а не независимо предсказан. Полиномиальные scalar-potential fits дали 4.02–6.50% RMS gradient mismatch с элементным P2 полем внутри `0.3r0`.

**Вывод:** numerical NO-GO для длинных траекторий и любых последующих RNG/entropy расчётов на текущем поле. Это не доказательство физической невозможности. Расширенные траектории не запускались. Далее требуется независимое уточнение central mesh, force symmetry/gradient recovery при неизменном DC voltage и только после успешного gate — staged deterministic integration. Полный отчёт: `M3A1_EXTENDED_DYNAMICS_VALIDATION_RU.md`.

Выполнена короткая детерминированная динамика (^{40}\mathrm{Ca}^+) по сохранённому M2b.2 P2 FEM полю: `r0=1µm`, `VRF=10V` differential peak, `q=0.2`, `fRF=2.472831GHz`, `VDC≈0.086979V` на каждом DC кольце для локальной осевой цели 1MHz. Полное time-dependent RF force интегрировалось в dimensionless coordinates; сравнивались RK4 и DOP853 с RF work accounting.

Короткие (20 RF-period) траектории для radial/axial/combined смещений и фаз оставались в `|r|<0.154r0`. Linear Floquet modes: две radial около `176.251MHz`, axial `1.0000005MHz`; Mathieu benchmark и unit tests прошли. Mesh comparison выбранных коротких trajectory endpoint дал менее `1.5e−5r0`.

Оговорка: при pointwise сравнении последнего mesh pair относительная неопределённость слабого DC gradient в центральной сфере составляет около 5–15%, хотя абсолютные различия малы; P2 gradients разрывны на element faces. Окно 8.1ns — только около 1.43 radial периода и 0.008 axial периода. Псевдопотенциальная разность на выбранной сфере `0.2r0` около `0.325µeV` — не escape depth. Поэтому результат — conditional deterministic GO, но stochastic entropy/RNG stage пока NO-GO до axial-cycle convergence и контролируемой field recovery. Подробный отчёт: `M3A_DETERMINISTIC_ION_DYNAMICS_RU.md`.

## Актуальное резюме R0 — 29 октября 2026

Аналитический скрипт сравнил 40Ca+, 88Sr+, 171Yb+, 9Be+ и электрон. Для идеального (r_0=1\,\mu m), (q=.2), 10 V сценария частоты RF примерно 2.47, 1.67, 1.20, 5.21 GHz и 0.667 THz. При фиксированных q/r0/V идеальная глубина одинакова 0.25 eV; при ограничении 100 MHz и пересчёте напряжения получаются Ca+ ~0.409 meV, электрон ~5.61 neV. Это аналитическая идеальная геометрия, не новый конечный FEM результат.

Главный риск RNG — конкуренция измеримости и сохранения термодинамического движения. FDT связывает thermal excitation с damping/readout loading; охлаждение, signal bandwidth, lifetime и power нельзя считать независимыми. Поверхностный шум для конкретного Si/oxide неизвестен и может одновременно давать entropy и heating/loss. Опубликован image-current detection электронов в room-temperature Paul trap, но не CMOS-integrated single-electron RNG. Новые отчёты и machine-readable реестр сохранены. R0 remote setup не исполнялся; первый ограниченный `.014` P2 run условен по RAM/disk и approval. M3a.3 numerical NO-GO для дальнейшей динамики остаётся в силе.

Публичный репозиторий создан: [ion-trap-rng-feasibility](https://github.com/snejokzboxik/ion-trap-rng-feasibility). R0.1 сверил public visibility, локальный и удалённый SHA, опубликовал инвентарь из 196 файлов и проверил настоящий свежий clone; 37 тестов прошли в clone. Публикация документационного handoff после этой проверки находится в текущем обновлении R0.1.

