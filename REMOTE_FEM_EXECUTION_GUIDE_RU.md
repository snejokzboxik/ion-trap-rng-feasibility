# Подготовка Windows-хоста для будущего bounded FEM запуска

Статус: подготовка локальных скриптов, удалённый запуск не выполнялся. Целевая машина Ryzen 5 7430U / около 16 GiB RAM может быть занята другой задачей; свободная память измеряется непосредственно в момент preflight. Все команды ниже выполняются только после отдельного разрешения на вычисления соответствующего этапа. Режим dry-run безопасен и не запускает FEM.

## Чистая установка

1. После публикации клонировать публичный репозиторий в отдельную папку проекта. До фактического GitHub push URL ещё не существует; не использовать родительский репозиторий `impedance_semiconductor_trap` как clone target.
2. Установить 64-bit Python 3.12 и Git для Windows. В PowerShell из корня проекта выполнить `.emotesetup_windows.ps1 -DryRun`, затем при готовности `.emotesetup_windows.ps1`. Установка выполняется в локальную `.venv`, lock содержит проверенные версии NumPy 2.5.3, SciPy 1.18.1, PyYAML 6.0.3, scikit-fem 12.0.2, Gmsh 4.15.2, Matplotlib 3.11.2 и pytest 9.1.1.
3. Проверить конфигурацию и тесты до mesh/solver запуска. Python lock фиксирует top-level пакеты, но не платформенные транзитивные wheels; для межкомпьютерной научной воспроизводимости фиксировать Python minor, `pip freeze`, OS/build и Gmsh бинарную версию в manifest после установки.

## Preflight и рекомендуемый первый FEM run

`configs/r0_remote_fem.json` задаёт один кандидат: существующая октантная mesh `plateau014`, DC и RF-basis P2 Laplace solve, один поток. Не строится новый более плотный mesh и не запускаются geometry sweep/траектории. Уже измеренные reference DOF: 45,824 для решённого `.018`, mesh-only 63,759 для `.014`. `.018` LU имел 32.33 млн ненулевых записей L+U; оценка для `.014` 50–80 млн записей была экстраполяцией, а не профилированием. При локальном RAM 1.28 GiB `.014` solve не запускался.

Конфигурация требует на машине минимум 13 GiB свободной RAM и 10 GiB диска, оставляя nominal 5 GiB на ОС/конкурирующую работу поверх modeled 8 GiB limit. Это консервативный операционный guard, не гарантия отсутствия пиков памяти; администратор может поднять thresholds. Если другой workload занят/память ниже guard, запуск отказать и перенести. **Первый overnight job:** только один DC/RF basis call на уже сгенерированной `.014` сетке, если измеренные ресурсы проходят preflight и имеется явное одобрение. Число повторов — ноль; после решения проверить manifest, CPU time, peak working set, LU nnz, force comparison и сохранность старых файлов. Затем остановиться для оценки, а не автоматически уточнять сетку.

Текущий P2 матричный оператор Laplace с зафиксированным grounded enclosure положительно определён после наложения Dirichlet DOF, так что CG с SPD-предобуславливателем физически/математически подходит. Наиболее перспективный sparse-memory кандидат — smoothed-aggregation algebraic multigrid (AMG) как preconditioner; SciPy CG API прямо требует Hermitian positive definite A/M. Однако **никакой iterative/AMG solve не валидирован для данной P2 mesh** и AMG может также хранить значительную иерархию. Перед применением к `.014` нужен небольшой контрольный benchmark на `.018`: сравнить CG+AMG с SuperLU по true residual (||b-Ax||/||b||), центральной DC силе/Hessian/flux, времени и пиковому RSS; назначить tolerances по observable sensitivity, а не только solver residual. Если амг-настройка не держит нужную силу, оставаться на direct solve/увеличить ресурс либо прекратить. Источники API: [SciPy CG](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html), [PyAMG smoothed aggregation](https://pyamg.readthedocs.io/en/latest/generated/pyamg.aggregation.html).

Команды:

```powershell
.\remote\setup_windows.ps1 -DryRun
.\remote\preflight.ps1 -DryRun
.\remote\run_fem.ps1 -DryRun
# Только после явного разрешения на запуск и пройденного ресурсного preflight:
.\remote\run_fem.ps1 -Execute
.\remote_collect_results.ps1 -Destination D:\iontrap-results -DryRun
```

Solver не поддерживает безопасный checkpoint внутри SuperLU factorization; прерывание собственного процесса может оставить неполный output, который verifier должен отклонить. Run wrapper наблюдает за рабочим набором собственного Python процесса и общей свободной RAM, прекращает только этот FEM-процесс при пересечении настроенных границ, не трогает фоновые задачи, не меняет их приоритет и не повторяет задачу. Это polling guard, не абсолютная защита от очень быстрого скачка памяти. Успех проверки файлов — проверка комплектности, не физической точности/сходимости.

## Известная граница

M3a.3 остаётся numerical NO-GO для extended dynamics: даже два решённых mesh уровня показывают до 2.51% DC-force изменения; полный и октантный fine results отличаются до 5.24%. `.014` решается только как потенциальный дополнительный уровень для диагностики. Адаптивный estimator должен информировать локальное уточнение по DC force functional, а не служить основанием для объявления результата достоверным без независимого cross-check.

Windows setup устанавливает зависимости из top-level version lock в локальную `.venv`. Сам src-layout пакет не ставится editable: тесты импортируют `src/` через `tests/conftest.py`, что избегает скачивания незафиксированной build-backend версии. Перед будущим чистым environment run следует повторно валидировать lock для Windows wheel availability.
