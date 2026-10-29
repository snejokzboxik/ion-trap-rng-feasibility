# Реестр физических механизмов и модельных пробелов

Машиночитаемый реестр — `PHYSICAL_EFFECTS_REGISTER.json`. Величины для reference (r_0=1\,\mu m), Ca+, 10 V RF peak и идеального проводника являются аналитическими порядками величины, а не результатами измерений/новой FEM. В частности, материалы и проводимость silicon/oxide не задавались. Каждая «можно опустить» формулировка условна: незнание величины не является доказательством её малости.

## Приоритеты

1. **Johnson–Nyquist / receiver noise — обязательный baseline**, включая прямое измерение шума схемы без trapped particle. Классическая односторонняя PSD (S_V=4k_BT\operatorname{Re}Z), когда (hf\ll k_BT); на высоких частотах/низкой T необходимо согласованное квантовое FDT выражение и конкретное соглашение PSD. Эта энтропия генерируется схемой и может быть эффективнее оцифрована напрямую; ион имеет преимущество только если докажется дополнительный measurable component.
2. **Тепловое движение иона с согласованным damping/cooling — условная главная гипотеза RNG.** Уравнение Ланжевена должно подчиняться fluctuation–dissipation theorem: сила возбуждения, сопротивление/охлаждение, receiver loading, мотионная температура и lifetime связаны. Слишком малое демпфирование снижает передачу/установление статистики, слишком большое демпфирование охлаждает/гасит сигнал; инженерное считывание имеет конечную полосу. Стохастический seeded integration не даёт физической энтропии без экспериментально ограниченной модели шума.
3. **Поверхностное поле, adsorbate и зарядовые флуктуации — потенциально значимый, но плохо переносимый на новую поверхность кандидат.** Их измеренная спектральная плотность зависит от материала, подготовки, температуры, частоты и расстояния. Это может давать heating и entropy, но одновременно разрушать lifetime/стабильность; для unbuilt nanoscale Si нельзя выбрать число по универсальному (d^{-4}) закону.
4. **Ион-газ столкновения — не основной vacuum RNG кандидат**, а вторичный сравнимый режим. Пуассоновское событие пригодно как физический процесс только при реальных partial pressures, cross-sections, потерях заряда и восстановлении. Столкновения могут вызывать потерю/charge exchange и non-Maxwellian tails, а не просто слабый белый шум.

## FDT и signal baseline

Для независимого электрического порта (S_V^{(1)}=4k_BT\operatorname{Re}Z) даёт при 300 K и 1 kΩ около 4.07 nV/√Hz до фильтра/передаточной функции. Это не равно напряжению на входе реального CMOS receiver. Ion image current по Shockley–Ramo: (i=Q\mathbf v\cdot\mathbf E_w), (\mathbf E_w=-\nabla\phi_w) при единичном bias соответствующей sensing terminal и остальных заземлённых. Предыдущая условная оценка (g_w=0.4/r_0) даёт для Ca+ при 1mK rms ток около 29 fA; 3D weighting field, портовое Z, bandwidth и обратное воздействие ещё должны быть рассчитаны. Один и тот же резистор/резонатор, который задаёт измеримое напряжение, демпфирует частицу и вносит флуктуации по FDT.

Шум поля может дать heating rate (в принятой Brownnutt односторонней PSD convention, малая micromotion): (\dot{\bar n}\simeq Q^2 S_E(\omega)/(4m\hbar\omega)). Для RF traps с micromotion учитываются sidebands/Floquet response; формула служит лишь первичным бюджетом. Min-entropy rate нельзя вывести из secular frequency, аналогового SNR или равновесной температуры без совместной модели тракта, дискретизации, корреляции, adversarial/technical component и утверждённой conditioning методики.

## Эффекты, охват и критерий исключения

JSON содержит для каждого эффекта origin/model, безразмерные группы, reference scale, evidence, uncertainty, условия допустимого упрощения и следующий расчёт/измерение. Приоритетные unresolved: (i) FEM force convergence (текущий numerical NO-GO), (ii) 3D differential weighting gradient/receiver transfer, (iii) реальные RF/DC losses и semiconductor complex conductivity, (iv) surface noise/charge drift и image potential, (v) ion cooling/loading/restart, (vi) sealed cavity pressure evolution and total module power. Менее первичные, но не «нулевые»: поле магнитов и токов, thermal strain, vibration и fabrication tolerances; включать только после получения диапазонов, но проверять критерии неравенств до их исключения.

## Источники с проверяемыми DOI

- Brownnutt et al., *Rev. Mod. Phys.* 87, 1419 (2015), [10.1103/RevModPhys.87.1419](https://doi.org/10.1103/RevModPhys.87.1419): обзор heating/surface-noise данных; не универсальный закон для нового Si.
- Turchette et al., *Phys. Rev. A* 61, 063418 (2000), [10.1103/PhysRevA.61.063418](https://doi.org/10.1103/PhysRevA.61.063418).
- Sedlacek et al., *Phys. Rev. A* 97, 020302 (2018), [10.1103/PhysRevA.97.020302](https://doi.org/10.1103/PhysRevA.97.020302): конкретный surface trap и расстояние-зависимость.
- DeVoe, *Phys. Rev. Lett.* 102, 063001 (2009), [10.1103/PhysRevLett.102.063001](https://doi.org/10.1103/PhysRevLett.102.063001): collision-driven non-Maxwellian regime.
- Ramo, *Proc. IRE* 27 (1939), [10.1109/JRPROC.1939.228757](https://doi.org/10.1109/JRPROC.1939.228757): induced-current theorem.
- Flowers-Jacobs et al., *Metrologia* 54, 730 (2017), [10.1088/1681-7575/aa7b3f](https://doi.org/10.1088/1681-7575/aa7b3f): Johnson noise thermometry.
- Электронная room-temperature Paul trap image current: Taniguchi & Noguchi (2025), [10.1103/6mpr-5s2h](https://doi.org/10.1103/6mpr-5s2h).
- Sealed miniature packages: Schwindt et al. (2016), [10.1063/1.4948739](https://doi.org/10.1063/1.4948739); Thrasher & Schwindt (2025), [10.1063/5.0243517](https://doi.org/10.1063/5.0243517). Эти достижения не означают CMOS-monolithic ion-RNG.

Проверенная проектная подборка также перечислена в `LITERATURE_AUDIT_RU.md`. Сырые snippets поисковой сессии не использовать как источники; соответствующий tool dump исключается из публикационного набора до правовой проверки.
