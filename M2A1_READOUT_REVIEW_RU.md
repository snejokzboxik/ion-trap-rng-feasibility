# Stage M2a.1 — обзор точности FEM и архитектуры считывания

**Дата:** 27 октября 2026 г.  
**Область:** ограниченный аудит сохранённой 2D модели, дополнительные лёгкие P1 проверки и аналитическое масштабирование lumped receiver. Полная 3D FEM, стохастическая динамика и опытный дизайн не выполнялись.

## Резюме и решение

M2a численно ответил на узкий вопрос о симметрии weighting field, но не обеспечил глобальную точность поля 0.1%. После исправления интегрирования элементной ошибки получено плато \(\|\mathbf E_h-\mathbf E\|_2/\|\mathbf E\|_2\simeq2.46\%\) для P1 сетки с неизменным максимальным размером дальних элементов \(h_{max}=0.4r_0\). Старые 1.73% получались центроидной квадратурой, пропускавшей внутритреугольную вариацию эталонного линейного градиента, и занижали истинную L2 ошибку. Точечное уменьшение \(h_{min}\) у центра почти не меняет глобальный интеграл. Независимая равномерная серия дала global error 3.319%, 1.664%, 0.833% при \(h/r_0=0.4,0.2,0.1\); trend близок к первому порядку. Это не даёт sub-0.1% гарантии и не оправдывает экстраполяцию на нанометровые материалы.

Физически наиболее обоснованный кандидат считывания — отдельный дифференциальный приёмный мод с частотой secular motion, развязанный от trapping RF mode, вместе с дифференциальной связью к четырём гиперболическим электродам. В 2025 г. image-current signal электронов действительно был зарегистрирован при комнатной температуре в Paul trap с двумя разнесёнными cavity modes; это подтверждает топологический принцип, но не интеграцию Ca+ в кремний. Обычный differential LC с переносом опубликованной добротности напрямую на CMOS был бы необоснованным. Решение: **условный GO только на один будущий 3D electrostatic/port benchmark**, после задания электрического стека и разводки; **NO-GO** на RNG-, энтропийные, мощностные или миниатюризационные утверждения.

## 1. Геометрия и соглашения

Сохраняется исходный поперечный RF профиль: четыре самостоятельных терминала на ветвях

\[
x^2-y^2=+r_0^2\;(R,L),\qquad x^2-y^2=-r_0^2\;(T,B),
\]

с гиперболическим обращённым к вакуумной апертуре лицом, противоположными парами RF и дифференциальным sense-портом. Идеальная опорная функция

\[
\Phi_{RF}=V_d\cos(\Omega t)\frac{x^2-y^2}{2r_0^2},\quad
\mathbf E=-\nabla\Phi_{RF}.
\]

Первичная среда — пустая герметичная вакуумная полость; кремний, SiO2, изоляторы и проводники находятся вне доступного иону объёма. Ambient/buffer gas здесь не моделировались. Существующая конечная 2D схема использует четыре узких гиперболических полосы, thickness \(0.1r_0\), extent \(1.5r_0\), sampled polygon boundary и grounded square. Уравнение weighting potential в пустой области было \(\nabla^2\phi_w=0\). При диэлектрическом окружении его надо заменить на

\[
\nabla\cdot[\epsilon(\mathbf r,\omega)\nabla\phi_w]=0,
\]

с непрерывностью \(\phi\) и нормальной компоненты \(\mathbf D\) на интерфейсе при отсутствии заданного свободного заряда. Случайный свободный поверхностный заряд — отдельный источник возмущений поля, не часть линейного weighting-field решения.

## 2. Аудит точности электростатики

### 2.1 Почему глобальная ошибка застряла около 2.46%

Идеальный manufactured-reference имеет квадрат \([-6r_0,6r_0]^2\), точный квадруполь на его внешней границе и P1 линейные треугольники. Локальная сетка мелкая у центра, но размер растёт до неизменных \(0.4r_0\) вдали. Поскольку точный потенциал квадратичный, ошибка потенциала мала; однако его точный градиент линейный, а P1 градиент внутри каждого элемента постоянен. Глобальный \(L_2\) также весит большую площадь дальнего поля. Поэтому последовательное сгущение около точки не устраняет дальнюю ошибку. Это численное следствие выбранной нормы и сетки, не физическая неоднородность.

| Сетка идеальной задачи | Узлы / треугольники | Глобальная ошибка E | E внутри \(\rho<r_0\) | Ошибка fitted Hessian в центре |
|---|---:|---:|---:|---:|
| adaptive, \(h_{min}=r_0/20\) | 2,141 / 4,120 | 2.4673% | 9.3535% | 0.7689% |
| adaptive, \(h_{min}=r_0/40\) | 2,240 / 4,318 | 2.4644% | 8.7502% | 0.4222% |
| adaptive, \(h_{min}=r_0/80\) | 2,350 / 4,538 | 2.4601% | 8.3718% | 0.3323% |
| adaptive, \(h_{min}=r_0/160\) | 2,475 / 4,788 | 2.4559% | 8.1649% | 0.3387% |
| adaptive, \(h_{min}=r_0/320\) | 2,610 / 5,058 | 2.4567% | 8.1231% | 0.1730% |
| uniform, \(h=0.4r_0\) | 1,124 / 2,126 | 3.3188% | 24.3608% | under-resolved fit |
| uniform, \(h=0.2r_0\) | 4,334 / 8,426 | 1.6643% | 11.6503% | \(1.94\cdot10^{-5}\%\)* |
| uniform, \(h=0.1r_0\) | 16,962 / 33,442 | 0.8326% | 5.7869% | \(6.21\cdot10^{-6}\%\)* |

\*Fit использует специально сбалансированный симметричный набор узлов и потому почти точно возвращает коэффициент параболы; это не общий bound на локальную ошибку поля и не подтверждение столь высокой accuracy конечного устройства. Для adaptive mesh Hessian least-squares fit зависит от распределения узлов и радиуса fit, поэтому наблюдаемая немонотонность на последних уровнях — диагностический сигнал.

Около линейный спад равномерной глобальной ошибки согласуется с \(O(h)\) градиента P1. Простая экстраполяция исправленного trend к 0.1% потребовала бы порядка \(h\sim0.012r_0\) по всему квадрату, то есть нескольких миллионов треугольников в 2D-области; это оценка стоимости, не выполненный расчёт и не error bound. Для локального решения разумнее использовать раздельные нормы по области интереса, адаптивный residual estimator и goal-oriented оценку градиента около иона.

### 2.2 Weighting field, кривизна, Hessian и multipoles

В прежней finite-electrode серии \(h_{min}/r_0=1/20\ldots1/320\): между двумя последними уровнями дифференциальный центр-gradient изменился на 0.060%, per-electrode \(C'\) — на 0.019%, fitted finite-drive quadrupole coefficient — примерно на 0.107%. Коэффициенты fitted higher multipoles на радиусе \(0.4r_0\) поменялись с 0.221% до 0.209% между \(1/160\) и \(1/320\) (около 5.4% относительно последнего значения). Поэтому малый численный multipole fraction пока нельзя выдавать за сходящуюся 0.2%-ную физическую характеристику.

Отдельное CAD sampling \(N=50,100,200\) точек на ветвь при \(h=r_0/160\) изменило differential gradient менее чем на 0.01%, а \(C'_{diff}\) — менее чем на 0.03%. Это ограничивает чувствительность интегральных центральных величин к выбранной polygon sampling density, но **не** заменяет curved-boundary elements: локальное поле у торцов, radius/roughness и curvature derivative могут иметь существенно большую ошибку. P2 isoparametric curved elements should be the next method comparison. В manufactured ideal case \(\Phi\) — полином второй степени, он содержится в triangular P2 space; при точной граничной интерполяции и точной квадратуре P2 должен воспроизвести идеальное поле и Hessian до линейного алгебраического допуска. Это теоретический benchmark-свойство, а не уже выполненный P2 run.

### 2.3 Чувствительность к заземлённому экрану

Для finite-electrode модели при \(h_{min}=r_0/80\), неизменном профиле и grounding на \(R=4,6,10r_0\) получено:

| \(R/r_0\) | drive Hessian coefficient (V/r0² при 10 V) | центр-градиент / квадруполь | \(|\nabla\phi_{w,d}|r_0\) | electrode \(C'\) (F/m) | differential \(C'\) (F/m) |
|---:|---:|---:|---:|---:|---:|
| 4 | 10.062904 | \(3.29\cdot10^{-4}\) | 0.399082 | \(7.6810\cdot10^{-11}\) | \(3.9386\cdot10^{-11}\) |
| 6 | 10.062829 | \(2.84\cdot10^{-4}\) | 0.399077 | \(7.2985\cdot10^{-11}\) | \(3.7511\cdot10^{-11}\) |
| 10 | 10.062916 | \(3.21\cdot10^{-4}\) | 0.399086 | \(7.1222\cdot10^{-11}\) | \(3.6729\cdot10^{-11}\) |

В выбранном идеальном vacuum cross-section центральный differential gradient и drive Hessian почти нечувствительны к расстоянию внешнего ground в этом диапазоне. Ёмкость заметно зависит от него: изменение \(C'\) между 4 и 10 \(r_0\) около 7.8%, потому что внешний ground физически меняет capacitance domain. Это **не** mesh-converged enclosure study и не модель конкретного корпуса. Упаковка, floating/grounded silicon, oxide supports и routing могут дать больший эффект; без реального стека field/capacitance нельзя ограничить.

### 2.4 Предел численного аудита

Проведены manufactured-solution comparison и h-series, но независимого residual-based estimator и P2/curved solve нет. Следовательно, строгое неопределённостное число для finite weighting gradient, Hessian или multipole fit пока не установлено. Значение 0.060% — разность двух сеток, а не доказанная полная относительная ошибка. На последующей FEM итерации требуется: P1/P2 сравнение; локальные и global refinement; fit-radius sweep; boundary/ground relocation; сравнение реакционного заряда, surface-flux integral и energy identity; Zienkiewicz–Zhu/recovery или residual estimator; отдельно error bands от идеализации торцов и шероховатости. Mesh size задавать отношением к \(r_0\), но проверять геометрическую адекватность самого континуума при \(r_0\le100\) nm.

## 3. Опубликованные опорные параметры receiver и применимость

Эти результаты дают **якоря отдельных параметров**, а не статистический интервал одного интегрируемого приёмника. Ни одна найденная публикация не демонстрирует выбранные четыре гиперболических электрода Ca+ + CMOS readout при 175 MHz.

| Источник / первичный результат | Опубликованные параметры | Что нельзя переносить на проект |
|---|---|---|
| Joshi et al., 65-nm CMOS charge buffer, [DOI 10.1109/ISCAS.2016.7539210](https://doi.org/10.1109/ISCAS.2016.7539210) | входная ёмкость 7 fF; площадь 193 µm²; unity gain; авторское обозначение power density 6 µW/MHz | Нет валидированной здесь input-noise density на 175 MHz, RF rejection, вакуумной упаковки или trapped-ion signal. 6 µW/MHz нельзя читать как power устройства без частотной/полосовой дефиниции из полного дизайна. |
| Liu et al., CMOS C–V ASIC для gyroscope, [DOI 10.1063/5.0186178](https://doi.org/10.1063/5.0186178) | 180-nm CMOS, area 0.039 mm²; bandwidth 6–8 kHz; output-noise density <150 nV/√Hz; resolution 840 zF/√Hz | Измеряемая ёмкостная модуляция медленная; эти цифры не задают шум/полосу на сотнях MHz/GHz. |
| Yang et al., released CMOS-MEMS spiral inductor, [DOI 10.3390/s111009798](https://doi.org/10.3390/s111009798) | measured \(L=1.8\) nH at 1 GHz; peak \(Q=15\) at 15 GHz; self-resonance 25.5 GHz; three-and-half-turn layout with 136 µm inner diameter; post-CMOS release/Si removal | Это MEMS-postprocessed inductor, не стандартная CMOS metal stack. Q=15 — максимум у 15 GHz, не значение при 175 MHz или 1.75 GHz. Никакого 0.8-µH on-chip component из этого не следует. |
| Chaddha et al., cryogenic Penning trap detector, [DOI 10.1007/s10909-025-03309-z](https://doi.org/10.1007/s10909-025-03309-z) | detector near 39.1 MHz; ~200 kHz band; Q up to 3500 at 4.2 K; ~1 mW tuned pHEMT amplifier; 40 dB gain; resonant electron-cloud signal observed | Cryogenic, external helical resonator, Penning electron cloud; no monolithic CMOS, Ca+ single ion, or 175 MHz validation. |
| Wang et al., superconducting Penning resonator/amp, [DOI 10.1007/s41365-025-01797-y](https://doi.org/10.1007/s41365-025-01797-y) | superconducting NbTi resonator \(Q=98{,}004\) near 1 MHz; GaAs amplifiers ~0.8 nV/√Hz at 1 MHz, >27 dB gain, ~6 mW, 0.1–10 MHz | Cryogenic/superconducting resonance and low frequency; extreme upper reference only, not plausible CMOS loaded Q at target frequency. |
| Taniguchi & Noguchi, room-temperature Paul electron image current, [DOI 10.1103/6mpr-5s2h](https://doi.org/10.1103/6mpr-5s2h) | separate cavity readout mode 619 MHz (\(Q_i\sim1300,Q_e\sim20{,}000\)); trap mode 3.105 GHz (\(Q_i\sim2300,Q_e\sim3300\)); room-temperature image-current readout; trap-mode 92 V and 1.5 W microwave input in current coupling (estimated 0.8 W at critical coupling) | Electron, 3D gold-plated copper coaxial cavities; mode-specific port placement, not four-hyperbola silicon electrode ASIC. High drive power illustrates system cost rather than a target. |

У Taniguchi & Noguchi loaded-Q оценки из \(Q_L^{-1}=Q_i^{-1}+Q_e^{-1}\) равны около 1,220 для readout mode и 1,356 для trap mode. Это реальная комнатно-температурная подсказка для **архитектуры** разнесённых мод и RF feedthrough rejection, а не гарантированный Q на нашем secular mode. Их эксперимент — наиболее прямое доказательство, что сильная trapping microwave и слабый image-current signal могут сосуществовать при специализированном cavity/port design; одновременно 3D cavity и сотни mW–W не соответствуют chip-level energy ambition.

Для RF rejection есть ещё одна разница полос: при выбранном \(q=0.2\), \(f_{sec}/f_{RF}\approx0.0707\) (около 14:1). Одного частотного разнесения недостаточно для оценки feedthrough: важны coupling matrix, common-mode conversion, фазовая асимметрия, нелинейная перегрузка входа и конечные Q фильтров. Никакой численный CMRR/phase-balance target для предлагаемого layout пока не подтверждён публикацией, поэтому он оставлен параметром будущего port-моделя.

## 4. Размерный screen для \(^{40}\mathrm{Ca}^+\)

Рассматривается **согласованная** идеальная серия при постоянных \(q=0.2\), differential peak \(V_d=10\) V, геометрическом отношении и \(A/r_0=\eta=0.1\). Тогда \(f_{RF},f_{sec}\propto r_0^{-1}\), \(A=0.1r_0\), а 2D capacitance для выбранной длины \(L_z=10r_0\) пропорциональна \(r_0\). Weighting gradient \(G_w\simeq0.39974/r_0\), поэтому

\[
I_{pk}=e(2\pi f_{sec})A G_w\propto r_0^{-1},\quad
E_\eta=\frac{q e V_d\eta^2}{8}.
\]

Это scaling law для идеальной модели, а не обещание, что фиксированные 10 V и q сохраняют физический материал, RF network или полную ловушку при нанометровом масштабе.

| \(r_0\) | \(A=0.1r_0\) | \(f_{sec}\) | \(f_{RF}\) | Ramo \(I_{pk}\) | differential trap \(C\), 2D при \(L_z=10r_0\) | ideal \(E_{0.1}\) / \(E_{0.5}\) |
|---:|---:|---:|---:|---:|---:|---:|
| 10 µm | 1 µm | 17.48 MHz | 0.2472 GHz | 0.703 pA | 3.75 fF | 2.50 / 62.5 meV |
| 3 µm | 300 nm | 58.26 MHz | 0.824 GHz | 2.34 pA | 1.125 fF | 2.50 / 62.5 meV |
| 1 µm | 100 nm | 174.78 MHz | 2.472 GHz | 7.03 pA | 0.375 fF | 2.50 / 62.5 meV |
| 300 nm | 30 nm | 582.61 MHz | 8.24 GHz | 23.45 pA | 0.1125 fF | 2.50 / 62.5 meV |
| 100 nm | 10 nm | 1.748 GHz | 24.72 GHz | 70.34 pA | 0.0375 fF | 2.50 / 62.5 meV |

При 300 K \(k_BT\approx25.85\) meV; радиальная \(\eta=0.5\) acceptance в идеальном pseudopotential — лишь 2.42 \(k_BT\). Это не escape barrier с RF phase, DC, axial end fields и collisions. Амплитуда \(A=0.1r_0\) — условная точка для сопоставления receiver, не равновесная амплитуда и не гарантированная динамика; без охлаждения крупные тепловые колебания и потеря должны оцениваться отдельно.

### 4.1 Частотная полоса, входная ёмкость и индуктивность

Для первого порядка использовано \(C_{tot}=C_{in}+C_{trap}\), \(|Z_C|=(2\pi f_{sec}C_{tot})^{-1}\). Значения \(C_{in}=7\) fF и 1.02 pF — опубликованные отдельными CMOS sensor designs якоря; это не measured parasitic range для данного chip. При \(r_0=1\) µm получаем \(|Z_C|=123\) kΩ и \(V_{pk}\approx0.868\) µV в capacitive limit для 7 fF; при 1.02 pF — 892 Ω и только 6.28 nV. По пяти масштабам это напряжение с 7 fF меняется примерно 0.60–0.91 µV, с 1.02 pF остаётся около 6.26–6.28 nV: рост Ramo current с уменьшением размера компенсируется ростом secular frequency и падением capacitive impedance. Эти числа не включают amplifier noise, actual electrode-network mismatch, RF leakage или ion back-action.

Для параллельной настройки на \(f_{sec}\):

\[
L=\frac{1}{(2\pi f_{sec})^2C_{tot}},\quad
R_p=\frac{Q_L}{2\pi f_{sec}C_{tot}},\quad
B_{3dB}\simeq\frac{f_{sec}}{Q_L},\quad
\tau_{amp}\simeq\frac{Q_L}{\pi f_{sec}}.
\]

При опубликованном CMOS input-cap anchor 7 fF, требуемый L уменьшается с 7.7 mH при 10 µm до 1.18 µH при 100 nm; даже самый малый вариант требует inductance далёкой от опубликованных nanohenry on-chip coils. Для 1.02 pF L меняется с 82.6 µH (10 µm) через 0.826 µH (1 µm) до 8.13 nH (100 nm). Уменьшение размера облегчает passive LC tuning только если сохраняется примерно пикофарадный load, но такой load уменьшает unmatched node impedance и подавляет сигнал.

| \(r_0\) | \(B_{3dB}\), \(Q=15\)* | \(B_{3dB}\), \(Q=1220\)** | \(B_{3dB}\), \(Q=3500\)*** |
|---:|---:|---:|---:|
| 10 µm | 1.165 MHz | 14.3 kHz | 4.99 kHz |
| 3 µm | 3.884 MHz | 47.8 kHz | 16.6 kHz |
| 1 µm | 11.65 MHz | 143 kHz | 49.9 kHz |
| 300 nm | 38.84 MHz | 478 kHz | 166 kHz |
| 100 nm | 116.5 MHz | 1.43 MHz | 499 kHz |

\* Только опубликованный peak Q=15 CMOS-MEMS inductor на 15 GHz; не переносить на эти частоты.  **\(Q_L\sim1220\) из room-temperature electron Paul-cavity readout mode at 619 MHz, также другая структура.  ***Q=3500 измерен в cryogenic Penning receiver около 39 MHz.

На одной и той же ловушке рост Q резко увеличивает амплитуду напряжения узкополосного resonant sensor, но сужает полосу и удлиняет ring-down. Поэтому SNR не даёт entropy rate: быстрые независимые инновации могут быть усреднены resonator memory; медленный процесс может остаться наблюдаемым. ENBW, tank bandwidth, sampling frequency и correlation/innovation rate не тождественны. Полный Lorentzian thermal voltage integral в narrowband parallel-RLC приближении равен \(\sqrt{k_BT/C}\); amplifier contribution и loaded ion-circuit back-action должны добавляться с реальными transfer functions. Текущая таблица не превращает выбранный Q в источник физической энтропии.

## 5. Архитектуры приёмника: сравнение

| Топология | Преимущество | Главные проблемы в этой ловушке | Оценка |
|---|---|---|---|
| Direct differential capacitive buffer / charge amplifier | Широкая полоса; опубликовано 7 fF CMOS input-C; проще сохранить быстрые инновации, чем высоким-Q tank | Нужны input-referred noise и CMRR на 175 MHz–1.75 GHz, не опубликованы для подходящего 7 fF макета; на 1pF feedthrough-limited signal около 6 nV при \(\eta=0.1\); RF trap source способен перегрузить вход через parasitics | Условная исследовательская базовая линия. При 7 fF даёт суб-microvolt узел без resonance, но реальный шум не установлен. |
| Differential parallel LC / resonant transimpedance | Хорошо известная ion-trap методика; чувствительность увеличивается как \(Q_L\); 39–40MHz ion/electron trap circuits уже реально детектировали слабый image signal | Для micron-scale и 7 fF нужны десятки/сотни µH; для 1 pF на 1 µm нужно ~0.83 µH, не подтверждённый integrated; высокий Q узит полосу; resonator связывает измерение с damping/ion back-action; winding/coil area может доминировать над trap | Нет доказательств integrated receiver на требуемом диапазоне; остаётся сравнительным LC benchmark. |
| Разнесённые trapping/readout modes в cavity, differential port | Прямой опыт Taniguchi & Noguchi показывает room-temp electron image current, separate 3.105 GHz trapping / 619 MHz readout mode и подавление утечки за счёт расположения портов | Cavity объемная и не CMOS; их RF drive потребовал ~0.8–1.5 W; Ca+ weighting field и четырёхгиперболическая cavity/electrode связка не тестировались | **Приоритетный кандидат для следующего физического port-моделирования**, потому что он прямо адресует RF feedthrough при image-current sensing. Интеграционная конкурентность не следует из доказанной detectability. |

Дифференциальный порт остаётся важным: в M2a идеальное связывание противоположных правой/левой ветвей подавляло центральный линейный weighting gradient примерно до \(4.4\cdot10^{-4}\) от one-electrode значения; отдельная differential sensing комбинация даёт \(0.39974/r_0\). Реальная фазовая/амплитудная ошибка двух ветвей нарушает подавление common-mode и добавляет trapping tone на приёмник. Электрическая точность split-drive в дискретном усилителе не равна CMRR всего электрод–кабель–амп–package пути.

## 6. Различие между абсолютной наблюдаемостью и entropy readout

Пусть стохастический процесс модулирует motion amplitude/phase с неизвестной скоростью инноваций \(\Gamma_{innov}\). Тогда резонатор с bandwidth \(f/Q\) и ring-down \(Q/(\pi f)\) фильтрует сигнал и способен уменьшить число различимых быстрых изменений. Обратная связь может внести damping и изменить само распределение motion. Уровень несущего image current или высокий coherent SNR не характеризует \(\Gamma_{innov}\), \(H_{min}\) или entropy per accepted sample. Нужны физическая модель bath/collision/surface fluctuations, transfer function до ADC, ion-present/absent subtraction, conditioning и leakage/side-information model. Это остаётся за рамками M2a.1.

## 7. Физические неопределённости и размерные границы

| Ограничение | Класс | Статус для данной работы |
|---|---|---|
| \(f_{RF},f_{sec}\propto1/r_0\), \(I\propto1/r_0\) при фиксированных q,V и A/r0 | F внутри идеальной quasistatic model | Аналитическая scaling relation; теряет силу при меняющейся проводимости, RF losses, геометрии и circuitry. |
| 7fF–1.02pF input-cap anchors | T / U | Параметры разных опубликованных CMOS сенсорных устройств; фактический портовый C, routing и package неопределён. |
| Q от нескольких десятков до ~10³ | T / U | Есть лишь межсистемные демонстраторы; нельзя объявить эту вилку диапазоном одного интегрированного приёмника. Q=3500 cryogenic Penning; 98k superconducting at 1MHz не кандидат для room-temperature CMOS. |
| Grounded shield and dielectric interface | T / U | Дистанционный экран влияет на capacitance около 8% в данном 2D toy enclosure; layout-реальные boundary condition/material properties отсутствуют. |
| Nanometer continuum/conductivity/roughness | T / U | При r0=100nm выбранная толщина кремниевого электродного backing 0.1r0=10nm; при r0=10nm она была бы 1nm. Классическая геометрия не доказывает локальную проводимость, непрерывность материала, стабильность оксида или малую шероховатость. |
| Useful entropy/throughput | R / U | Задаются приложением и физическим stochastic process; из сигналового SNR и полосы не получены. |
| Полный module size/energy | E / U | RF generator, isolation, receiver, cavity, source/loading, packaging отсутствуют в present calculation. |

Ни один из этих результатов не задаёт универсальный минимальный \(r_0\). По мере уменьшения масштаба field grows как \(V/r_0\), secular and trapping frequencies растут, но electric signal при согласованном \(A/r_0\) растёт лишь \(1/r_0\); пассивное tuning и RF feedthrough меняются по-разному в зависимости от параситных ёмкостей. На nanoscales material model становится доминирующей нерешённой неопределённостью.

## 8. Узко определённый будущий 3D расчёт

После отдельной авторизации предлагается **один 3D open/cavity-port benchmark при \(r_0=1\) µm для Ca+**, не sweep:

1. Сохранить четыре гиперболические branches, конечную толщину, определённую конечную длину (стартовый datum \(L_z=10r_0\)) и существующую осевую DC структуру; электродные terminal nets отдельно доступны для differential sense.
2. Центральная область — вакуумный empty cavity. Вокруг задать только один конкретный cross-section стек (Si electrode/support, SiO2 insulation, grounded metal enclosure) с геометрией, толщинами и conductivity/permittivity из выбранного процесса/источника. Не заменять полость сплошным диэлектриком.
3. Рассчитать weighting potential каждого из четырёх терминалов, Maxwell capacitance matrix и differential Ramo coefficient; отдельно RF trap field. Сравнить mode-selective sense port с простым differential lumped портом и вычислить coupling/feedthrough между портами при 175 MHz readout / 2.472 GHz drive.
4. Замкнуть только эквивалентный lumped network: измеряемые/обоснованные \(C_{in}\), \(Q_L\), amplifier noise PSD, source impedance, RF-isolation transfer. Пока конкретный CMOS macro не выбран, оформить explicit parameter intervals вместо выдуманной номинальной схемы.
5. Mesh convergence: global and local h, P1/P2, curved boundaries; independent charge/energy checks; material-boundary refinement; enclosure relocation. Проверить open-circuit response and termination sensitivity. Не вычислять stochastic entropy или битовый поток.

**Gate:** считать model valuable только если differential signal transfer, C matrix и RF leakage converged в выбранной неопределённости, и существует хотя бы один физически реализуемый parameter interval из измеренных схемных источников. Даже успешный gate позволит лишь перейти к анализу coupled ion damping/noise; он не докажет RNG.

## 9. Воспроизводимость и файлы

Результаты limited screens сохраняются в `results/m2a1_numeric_audit.json` и `results/m2a1_scale_screen.json`. Их расчётные сценарии воспроизводятся через `scripts/run_m2a1_bounded_audit.py` и `scripts/run_m2a1_scale_screen.py`. M2a FEM source remains `src/iontrap_rng/readout_fem.py`; входные значения прежних расчётов — `configs/m2a_readout.yaml`. Числа отчёта относятся только к эталонной аналитической Ca+ scaling и простой lumped circuit, не к измеренному модулю.

## Источники

1. Taniguchi & Noguchi, “Image current detection of electrons in a room-temperature Paul trap,” *Phys. Rev. A* 112, 022420 (2025), [DOI 10.1103/6mpr-5s2h](https://doi.org/10.1103/6mpr-5s2h). Прямое image-current/cavity-mode доказательство, но другая частица и геометрия.
2. Joshi, Kim & Cauwenberghs, “A 6 μW/MHz charge buffer with 7 fF input capacitance in 65nm CMOS for non-contact electropotential sensing,” ISCAS 2016, [DOI 10.1109/ISCAS.2016.7539210](https://doi.org/10.1109/ISCAS.2016.7539210). CMOS input-cap/area anchor.
3. Liu et al., “A low noise capacitance to voltage readout circuit based on active feedback charge sensitive amplifier for micro-gyro,” *Rev. Sci. Instrum.* 95, 065001 (2024), [DOI 10.1063/5.0186178](https://doi.org/10.1063/5.0186178). Низкочастотный CMOS sensor readout, не RF-ion noise specification.
4. Yang, Dai & Hong, “Manufacture and Characterization of High Q-Factor Inductors Based on CMOS-MEMS Techniques,” *Sensors* 11, 9798–9806 (2011), [DOI 10.3390/s111009798](https://doi.org/10.3390/s111009798). Measured inductor/Q/SRF, requiring MEMS post-process.
5. Chaddha et al., “Development and Characterization of a Cryogenic Tuned Amplifier for Penning Trap Applications,” *J. Low Temp. Phys.* (2025), [DOI 10.1007/s10909-025-03309-z](https://doi.org/10.1007/s10909-025-03309-z). Cryogenic comparator for measured ion-trap resonant receiver.
6. Wang et al., “Development and characteristics of a superconducting resonator and cryogenic low-noise amplifiers for Shanghai Penning Trap,” *Nucl. Sci. Tech.* 36, 217 (2025), [DOI 10.1007/s41365-025-01797-y](https://doi.org/10.1007/s41365-025-01797-y). Extreme cryogenic Q/noise upper reference.
