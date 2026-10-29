# Сравнение кандидатов частиц: ограниченный аналитический скрининг

Дата: 29 октября 2026. Расчёты воспроизводятся `scripts/compare_species.py`; здесь нет новых FEM-полей, траекторий или экспериментальных параметров. Использованы нейтральные изотопные атомные массы NIST; вычитание массы электрона для положительного иона изменяет значения менее чем на несколько (10^{-5}) относительно и не меняет округлённые выводы. Массы не следует путать с экспериментально измеренной массой конкретного ионного состояния.

## Допущения и два сценария

Общий идеальный quadrupole имеет (r_0=1\,\mu m), (|q|=0.2), одинарный заряд и small-(q) формулы. RF дифференциальное напряжение — peak. Для аналитической модели (Omega=\sqrt{2|Q|V_{RF}/(m r_0^2|q|)}), (f_{sec}\simeq |q|f_{RF}/(2\sqrt2)), а радиальная pseudopotential на (r=r_0) равна (|q|eV_{RF}/8). Это геометрическая оценка идеального поля; физическая поверхность выхода и DC ось имеют другую структуру.

| Частица | Сценарий A: (f_{RF}), фиксированы (V=10\,V, |q|,r_0) | (f_{sec}) | Идеальная (U_{ps}(r_0)) | Сценарий B (V) при 100 MHz | (U_{ps}(r_0)) при 100 MHz | Сценарий B (V) при 1 GHz | (U_{ps}) при 1 GHz |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40Ca+ | 2.473 GHz | 174.9 MHz | 0.250 eV | 16.35 mV | 0.409 meV | 1.635 V | 40.9 meV |
| 88Sr+ | 1.667 GHz | 117.9 MHz | 0.250 eV | 35.97 mV | 0.899 meV | 3.597 V | 89.9 meV |
| 171Yb+ | 1.196 GHz | 84.55 MHz | 0.250 eV | 69.94 mV | 1.749 meV | 6.994 V | 174.9 meV |
| 9Be+ | 5.208 GHz | 368.2 MHz | 0.250 eV | 3.687 mV | 0.0922 meV | 0.369 V | 9.22 meV |
| electron | 667.5 GHz | 47.20 GHz | 0.250 eV | 0.2245 µV | 5.61 neV | 22.45 µV | 0.561 µeV |

Округление чисел скрипта; сценарий B на 1 GHz — engineering comparison point, а не универсальный предел источника. Сценарий B фиксирует (q,r_0,f_{RF}) и пересчитывает напряжение; он показывает, что при ограничении RF частоты глубина масштабируется как (m f_{RF}^2). Радиальная глубина электрона при таких частотах чрезвычайно мала. Сценарий A, наоборот, сохраняет идеальную глубину, но требует для лёгких частиц очень высокую RF частоту и усиливает требования к паразитным модам, изоляции, RF-потерям и электронике. Номинальное поле (V/r_0=10\,MV/m) в сценарии A одинаково для всех, но не учитывает усиление на краях.

При неизменной DC кривизне, соответствующей (f_z=1\,MHz) для Ca+, идеальные осевые частоты составляют 1.000, 0.674, 0.484, 2.106 и 269.9 MHz соответственно. Это масштабирование по массе, не пересчёт напряжения по конечному FEM полю. Для Ca+ axial rms тепловая амплитуда при 300 K на 1 MHz равна около 40.7 µm, то есть в десятки раз больше (r_0); чтобы классическая тепловая амплитуда была порядка (0.1r_0), требуется милли-кельвиновая температура движения, около 1.9 mK. При этом нагрев/охлаждение, частотный диапазон и image-current readout связаны, а не регулируются независимо.

Для условного радиального single-mode readout при (g_w=|\nabla\phi_w|=0.4/r_0), взятом как оценка чувствительности из прежней ограниченной модели, (I_{rms}=|Q|g_w\sqrt{k_BT/m}). В 1mK это примерно 29 fA для Ca+, 20 fA для Sr+, 14 fA для Yb+, 62 fA для Be+ и 7.9 pA для электрона. Это не реальный сигнал 3D электрода: weighting gradient, mode, receiver transfer и взаимодействие с резонатором не проверены для этих частиц. Электронный большой ток — не автоматическое преимущество: trapping RF в сценарии A около 0.67 THz; сценарий B даёт крайне малую глубину, а image interaction становится существенной.

## Практический вывод по кандидатам

- **40Ca+** остаётся главным сравнительным reference из-за доступности лазерного охлаждения/флуоресцентного обнаружения и сохранённого M3a конфигурационного выбора. Не означает преимущества для электрического термошумового RNG: лазерная система и охлаждение создают значительные габариты/энергетические затраты.
- **88Sr+ и 171Yb+** уменьшают требуемую частоту при фиксированном (q,V,r_0), но тяжелее и снижают условный тепловой image-current при одинаковых (T) и weighting gradient. Сохраняются лазерные и вакуумные требования. Пакетные/часовые архитектуры с Yb+ продемонстрированы значительно крупнее chip module.
- **9Be+** повышает secular frequency и условный тепловой ток относительно Ca+, но для того же (q,V,r_0) требует около 5.2 GHz RF. У охлаждения/лазерного источника нет очевидного интеграционного преимущества.
- **Электрон** — интересный контрастный объект для microwave image-current detection, но не очевидный интегрируемый кандидат. При высоком RF частота попадает в sub-THz, а при инженерно ограниченной RF частоте псевдопотенциальная глубина опускается ниже neV–µeV. Плоская оценка image charge у идеальной плоскости при (d=1\,\mu m) около −0.36 meV (при 100 nm — −3.6 meV), что намного выше сценария B 1GHz глубины; это приближение и не полный барьер, но сигнализирует о критичности поверхности. Опубликована microwave image-current detection электронов в комнатном Paul trap (Taniguchi & Noguchi, PRA 112, 022420, 2025, DOI 10.1103/6mpr-5s2h); публикация говорит о trapped electrons и потенциальной возможности single-electron detection/cooling, но не доказывает CMOS-integrated, устойчивый single-electron RNG.

**Рекомендация:** сохранить Ca+ как траекторный/проверочный reference, а Sr+/Yb+ как варианты снижения RF; Be+ только если высокочастотная электроника и охлаждение подтверждены. Электрон не выбирать без отдельного анализа image-potential + finite RF capture + 100-GHz-class readout. Ни один кандидат пока не имеет положительного решения по полному RNG модулю.

## Загрузка, охлаждение и считывание — системный контекст

| Кандидат | Установленные способы захвата/инициализации (общая научная база) | Охлаждение/обнаружение | Что мешает автономному chip restart |
|---|---|---|---|
| Ca+ | Нейтральный Ca из oven/атомного источника, photoionization; электронная бомбардировка также применялась | Doppler laser cooling и fluorescence; возможен image-current, но сигнал/ион в малой ловушке остаётся неподтверждённым | Laser wavelengths, optics/detector, atom source, sealed vacuum/getter; атомный micro-oven опубликован как component, не в данной монолитной ловушке |
| Sr+ | Oven + photoionization | Laser cooling/fluorescence; electronic readout needs mode-sensitive resonator | несколько лазеров/оптический путь и вакуумная система |
| Yb+ | Atom source + photoionization | Laser cooling/fluorescence; sealed-package research uses Yb+ | isotope/laser/source and package; demonstrated packages are cm³ scale, not CPU die |
| Be+ | Atomic source + ionization | Laser cooling; established ion species but no inherent simplified detector | species-specific atom handling and laser subsystem; no chip-level source evidence in cited scope |
| electron | Electron gun/field emission/photoionization-type sources in trap literature | Microwave cavity image current; recent room-temperature Paul trap work demonstrates electron signal detection, while single-electron operation/cooling is discussed as future direction | sub-THz RF for high depth under scenario A, surface/image forces, source capture, microwave isolation and vacuum; not integrated module evidence |

These are architecture classes, not claims that the full listed loading chain is fabricated in a CMOS-compatible package. Before an RNG system claim, autonomous ion-loss detection, reloading, charge-state verification, startup time, source lifetime and recovery energy need an explicit budget. A trapped particle that disappears on power loss or shutdown makes operation non-autonomous unless a retained source and capture sequence exist.

Readout is also species/mode-dependent. Fluorescence may be strong but requires optical cycling and laser infrastructure; image current avoids photons but is often fA-scale for ions and needs a narrow-band resonator/amplifier. Narrow bandwidth can improve voltage sensitivity while extending correlation time and coupling the resonator thermally to the mode. An electron signal can be larger at equal temperature and gradient but its trapping frequency and surface interaction are severe constraints. Every candidate must be compared to a direct receiver Johnson-noise sample without the ion.

For practical comparison, a published CMOS latch-based TRNG reports 0.116 pJ/bit including its stated von Neumann conditioning (Wang et al. 2024, DOI 10.1109/TVLSI.2023.3328602; see literature audit). This is not a matched benchmark: process, throughput, health tests, energy boundary and entropy model must be aligned before numerical comparison. No ion module energy/bit is presently available. It must include RF generation and resonant loss, vacuum package/getter, particle source/reloading, cooling, analog front-end, conditioning and idle power. The ion proposal has no demonstrated energy, area, restart, or reliability advantage over direct CMOS thermal noise.

## Источники

- NIST isotopic composition tables, Ca/Sr/Yb; [NIST isotope masses](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii2). Масса 9Be: NIST Atomic Weights and Isotopic Compositions. Массы в этом скрипте — neutral isotope reference.
- K. Taniguchi & A. Noguchi, “Image current detection of electrons in a room-temperature Paul trap,” *Phys. Rev. A* 112, 022420 (2025), [DOI 10.1103/6mpr-5s2h](https://doi.org/10.1103/6mpr-5s2h).
- Для image current / weighting field и ion heating: см. проверенную подборку DOI в [LITERATURE_AUDIT_RU.md](LITERATURE_AUDIT_RU.md), в частности Ramo (1939), DOI 10.1109/JRPROC.1939.228757, Brownnutt et al. (2015), DOI 10.1103/RevModPhys.87.1419.

### Валидность

Известная физика: идеальный Mathieu scaling, тепловой гармонический осциллятор, поверхностная image-charge модель. Условные приближения: квадруполь бесконечной длины, small-q secular approximation, идеальный проводник, заданная weighting field. Гипотезы для будущего: какой вид частицы имеет лучший полезный entropy/readout баланс в реальном модуле. Неизвестные параметры: вакуумная упаковка и загрузка, охлаждение/нагрев, поверхность, 3D detector coupling, receiver noise/bandwidth, RF isolation и system power.
