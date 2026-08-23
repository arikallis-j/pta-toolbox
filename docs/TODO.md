# TODO

Да. Я бы здесь довольно существенно перестроила исходный план. Главная идея: **не пытаться сначала довести до конца весь задуманный экосистемный проект**, а выделить минимальное научное ядро, которое уже само по себе является новой методологией и может стать первой статьёй.

Ваша исходная точка хорошо ложится на существующую литературу: Gair et al. показали разложение PTA-отклика на gradient/curl-моды и принципиальную нечувствительность PTA к curl-компоненте; при этом их формализм работает непосредственно с гармоническими компонентами фона. ([arXiv][1]) Ваше потенциальное описание, судя по изложению, добавляет другой уровень интерпретации: **данные PTA интерпретируются как локальное/сглаженное измерение скалярного градиентного потенциала**, а не просто как набор измеренных ($a_{\ell m}$). Именно это я бы сделала центром первой работы.

Есть ещё один важный момент: современные исследования уже рассматривают восстановление анизотропий через Fisher/covariance formalism и оптимизацию чувствительности PTA, поэтому оптимизационная часть вашей программы имеет самостоятельную научную ценность, но её лучше отделить от собственно новой геометрической интерпретации. ([APS Journals][2])

---

# 1. Как я бы устроила весь проект

Я вижу здесь **три уровня**:

$$
\boxed{\text{новый формализм}}
\quad\longrightarrow\quad
\boxed{\text{проверка и восстановление}}
\quad\longrightarrow\quad
\boxed{\text{оптимизация PTA}}.
$$

При этом реальная работа может быть разбита примерно на **11 научных этапов**, плюс отдельный подготовительный этап.

---

# Этап 0. Подготовительный

Это то, что я сознательно **не считала бы научными результатами**.

### Что сделать

1. Привести существующие теоретические выкладки в единую систему обозначений.
2. Перепроверить все коэффициенты, нормировки и фазовые соглашения.
3. Разобраться с `s2fft`, переходами

$$
\text{map}\leftrightarrow a_{\ell m},
$$

включая spin-weighted transforms.
4. Реализовать базовые операции:

* scalar SHT;
* spin-2 SHT;
* gradient/curl decomposition;
* переходы между поляризациями и потенциалами;
* SVD матрицы отклика;
* построение собственных мод PTA.

5. Создать маленький воспроизводимый numerical sandbox.
6. Привести существующий текст к состоянию, в котором каждое утверждение можно однозначно связать с формулой или численным тестом.
7. Написать минимальную версию введения и literature review.

### Результат

Не статья и не новый научный результат, а **рабочая исследовательская инфраструктура**:

```text
potential.py / Julia equivalent
polarization.py
spherical_harmonics.py
pta_response.py
svd_basis.py
tests/
examples/
```

и единая математическая конвенция.

### Важная оговорка

Сама разработка кода — не научный результат. Но **результаты, которые этот код позволяет получить**, уже являются научной работой.

То есть:

> «Я освоила s2fft» — подготовка.

но:

> «На $10^4$ случайных конфигурациях я численно проверила тождество между polarization и potential representations с относительной ошибкой $<10^{-X}$» —

уже результат статьи.

---

# Этап 1. Зафиксировать математическое ядро формализма

Это, на мой взгляд, **самый первый настоящий научный этап**.

### Задача

Полностью вывести цепочку

$$
h_{+,\times}
\longleftrightarrow
a_{\ell m}^{G,C}
\longleftrightarrow
\Phi_G,\Phi_C
\longleftrightarrow
r_a
$$

и явно показать, какие именно операторы связывают эти величины.

Особенно важно разделить три утверждения:

1. разложение GW поля;
2. определение спиновых потенциалов;
3. выражение PTA response через потенциал.

И отдельно аккуратно установить:

$$
r_a^{\rm E} = \mathcal{R}_a[\Phi_E]
$$

и в идеальном случае

$$
r_a^{\rm ideal}\propto
\Phi_E(\hat{\mathbf p}_a).
$$

### Результат

Получается **замкнутая аналитическая теория**, где любую формулу из дальнейшей статьи можно получить из нескольких базовых уравнений.

Это фактически skeleton первой статьи.

---

# Этап 2. Численная верификация эквивалентности представлений

До сложных источников я бы обязательно сделала очень скучный, но исключительно важный тест.

### Задача

Сгенерировать произвольный набор

$$
a_{\ell m}^{E,B}
$$

и получить два независимых описания:

$$
(E/B)\rightarrow h_{+,\times}\rightarrow r_a
$$

и

$$
(E/B)\rightarrow
(\Phi_E,\Phi_B)\rightarrow r_a.
$$

Для одной и той же конфигурации пульсаров сравнить:

$$
r_a^{(1)}-r_a^{(2)}.
$$

### Результат

График вроде

$$
\frac{|r^{(1)}-r^{(2)}|}
{|r^{(1)}|}
$$

как функция ($L_{\max}$), числа пульсаров и случайных реализаций.

Нужна практически машинная проверка:

$$
\boxed{
\text{polarization formalism}
\equiv
\text{potential formalism}
}
$$

в пределах численной ошибки.

### Почему это раньше всего

Потому что иначе дальше вы одновременно отлаживаете:

* математику;
* (s2fft);
* PTA response;
* spin harmonics;
* SVD.

А тогда совершенно невозможно понять источник ошибки.

---

# Этап 3. Полностью исследовать PTA beam / диаграмму направленности

Это, на мой взгляд, **главный отдельный теоретический результат первой статьи**.

Вы уже фактически нашли очень красивую физическую интерпретацию:

$$
r_a
=
\int d\Omega',
K_a(\hat{\mathbf n},\hat{\mathbf n}')
\Phi_G(\hat{\mathbf n}').
$$

Причём

$$
K_a\rightarrow
\delta(\hat{\mathbf n}-\hat{\mathbf p}_a)
$$

в некотором пределе.

### Задача

Нужно:

1. вывести аналитический вид ($K$);
2. разложить

$$
K(\gamma)=
\sum_\ell K_\ell
P_\ell(\cos\gamma);
$$

3. определить коэффициенты ($K_\ell$);
4. сравнить их с коэффициентами сферической дельта-функции;
5. найти характерный ($L_\mathrm{break}$);
6. исследовать зависимость от:

   * ($L_{\max}$);
   * расстояния до пульсара;
   * длины GW;
   * количества используемых мультиполей.

### Результат

Это уже отличный набор figures:

**Figure 1:** ($K(\theta)$) для нескольких ($L_{\max}$).

**Figure 2:** ($K_\ell$) против ($\ell$) и соответствующие коэффициенты ($\delta$)-функции.

**Figure 3:** переход

$$
\text{quasi-delta}
\rightarrow
\text{oscillatory kernel}
\rightarrow
\text{butterfly}.
$$

**Figure 4:** ошибка локального приближения как функция параметров PTA.

Это очень хорошая физическая история.

---

# Этап 4. Показать предел малости пульсарного члена

Это следует непосредственно из предыдущего этапа.

### Задача

Не просто написать:

> pulsar term is negligible.

а сформулировать **критерий**, при котором это верно.

Например, некоторый безразмерный параметр вида

$$
\epsilon_{\rm PT} =
F(
L_p,\lambda_{\rm GW},\theta,
L_{\max},\ldots).
$$

И затем определить область:

$$
\epsilon_{\rm PT}\ll1.
$$

### Результат

У вас появляется:

$$
\boxed{
\text{criterion for the local-response approximation}
}
$$

а не просто демонстрационный график.

И это особенно ценно, потому что тогда фраза «для текущего PTA можно пренебречь pulsar term» превращается в проверяемое утверждение.

---

# Этап 5. Обобщение детерминированного источника

Теперь можно переходить к физическим источникам.

Я бы **первым рассматривала детерминированный источник**, а не stochastic background.

Вы уже имеете L-polarization. Следующий шаг:

$$
h_+,\ h_\times
$$

или эквивалентные потенциальные компоненты с произвольными амплитудами и фазами.

### Задача

Получить:

$$
\Phi_G(\hat n;\theta_s,\mathcal A,\psi,\phi_0),
$$

а затем PTA response.

Проверить:

$$
\text{source parameters}
\rightarrow
\Phi_G
\rightarrow
r_a.
$$

### Результат

Аналитическая формула для arbitrary-polarization point source.

После этого можно показать очень красивую вещь:

$$
\text{точечный источник}
\rightarrow
\text{потенциальная карта}
\rightarrow
\text{локальный максимум}.
$$

---

# Этап 6. Восстановление потенциала как обратная задача

Вот здесь появляется вторая половина вашей идеи.

Вы уже предлагаете:

$$
\mathbf r
\rightarrow
\text{PTA eigenmodes}
\rightarrow
\text{mode amplitudes}
\rightarrow
a_{\ell m}
\rightarrow
\Phi_G.
$$

Это надо сделать отдельным полноценным numerical experiment.

### Задача

Построить матрицу отклика

$$
R_{a,\ell m}.
$$

Сделать SVD:

$$
R=U\Sigma V^\dagger.
$$

Определить измеряемые и плохо измеряемые направления в пространстве функций.

Затем реконструировать

$$
\hat\Phi_G.
$$

### Результат

Нужно показать:

$$
\Phi_G^{\rm true}
\rightarrow
{r_a}
\rightarrow
\hat\Phi_G
$$

и измерить качество:

$$
\mathrm{MSE},
\qquad
C_\ell,
\qquad
\text{angular localization error}.
$$

Очень полезно показать, что **число реально восстанавливаемых мод определяется не просто ($N_\mathrm{pulsar}$), а спектром singular values**.

Это может оказаться одним из наиболее сильных результатов всей работы.

---

# Этап 7. Стохастический фон

Теперь ваш stochastic case естественно появляется как статистическое обобщение предыдущей задачи.

### Задача

Переписать stochastic GW background через потенциальные коэффициенты:

$$
\langle a_{\ell m}a_{\ell' m'}^*\rangle =
C_{\ell m,\ell'm'}.
$$

Для анизотропного фона получить связь с covariance PTA response.

Здесь как раз появляется machinery с 3j-symbols, используемая в направлении Gair et al. ([arXiv][1])

### Результат

Должна получиться схема:

$$
\boxed{
\text{GW covariance}
\leftrightarrow
\text{potential covariance}
\leftrightarrow
\text{PTA covariance}
}
$$

и демонстрация на:

1. isotropic background;
2. dipole/quadrupole anisotropy;
3. локальный hotspot.

### Очень важное замечание

Я бы **не пыталась сразу строить полный realistic stochastic Bayesian pipeline**.

Сначала:

$$
\text{Gaussian signal + Gaussian noise}
$$

и только затем добавлять:

* white noise;
* red noise;
* timing model;
* realistic frequencies.

Иначе статистическая часть моментально станет отдельной огромной статьёй.

---

# Этап 8. End-to-end восстановление

Теперь собрать всё вместе.

### Задача

Сделать synthetic PTA:

$$
\Phi_G
\rightarrow
r_a
\rightarrow
r_a+n_a
\rightarrow
\hat\Phi_G.
$$

Рассмотреть два случая:

### A. Deterministic source

$$
\hat\Phi_G
\rightarrow
\hat{\theta}_s,\hat{\mathcal A},\ldots
$$

### B. Stochastic background

$$
\hat C_{\ell m}
\rightarrow
\text{anisotropy reconstruction}.
$$

### Результат

End-to-end демонстрация работоспособности метода.

Здесь уже нужны realistic noise models:

$$
N=N_{\rm white}+N_{\rm red}.
$$

Это значительно усиливает статью, потому что показывает не просто математическое тождество, а реально работающий inference method.

---

# Этап 9. Совместный поиск deterministic source + stochastic background

Это уже **следующий уровень сложности**.

### Задача

Проверить, можно ли одновременно оценивать

$$
\theta_{\rm CW}
$$

и

$$
C_{\ell m}
$$

без существенного bias.

Например,

$$
p(d|\theta_{\rm CW},C_{\ell m},N).
$$

Сделать injection/recovery.

### Результат

Матрица:

$$
\text{true parameters}
\rightarrow
\text{posterior}
$$

и исследование degeneracies.

### Мой совет

Этот этап **не должен блокировать первую статью**.

Это отличная часть второй статьи.

---

# Этап 10. Оптимальная конфигурация PTA

Вот здесь начинается уже самостоятельная исследовательская программа.

Причём формулировка должна быть немного иной.

Не:

> «найти красивое расположение пульсаров».

А:

> Для заданного класса сигналов найти конфигурацию PTA, максимизирующую определённую информационную/статистическую меру.

Например:

$$
\mathcal I =
F(\theta)
$$

или

$$
\det F,\qquad
\operatorname{Tr}F,\qquad
\lambda_{\min}(F),
$$

где ($F$) — Fisher matrix.

### Задача

Для одного и того же алгоритма сравнить:

$$
\begin{array}{c|c}
\text{source} & \text{optimal PTA} \\
\hline
\text{isotropic} & \text{uniform} \\
\text{point source} & \text{concentrated} \\
\text{cosmic string} & \text{extended patch}
\end{array}
$$

### Результат

Карта:

$$
\boxed{
\text{source morphology}
\rightarrow
\text{optimal pulsar geometry}
}
$$

И вот это уже очень сильный отдельный результат.

При этом это направление имеет заметную связь с существующими задачами оптимизации PTA и выбора пульсаров, поэтому нужно будет особенно чётко показать отличие вашей objective function и физической постановки. ([OUP Academic][3])

---

# Этап 11. Реальное распределение пульсаров и «галактический базис»

Я бы сознательно поставила это **последним**.

### Задача

Вместо идеального распределения

$$
{\hat p_a}_{a=1}^{N}
$$

взять реальную галактическую выборку.

Рассмотреть SVD:

$$
R_{\rm uniform}=U\Sigma V^\dagger
$$

против

$$
R_{\rm Galactic}=U_G\Sigma_GV_G^\dagger.
$$

Исследовать:

$$
\Delta U,\qquad
\Delta\Sigma.
$$

### Результат

Ответ на вопрос:

> насколько собственные моды PTA универсальны, а насколько они являются свойством конкретной геометрии нашей Галактики?

А затем можно определить:

$$
\boxed{
\text{Galactic PTA basis}
}
$$

как адаптированный basis функций.

Но это **совершенно необязательно для первой статьи**.

---

# 2. Что из этого подготовка, а что наука?

Здесь я бы провела довольно жёсткую границу.

| Работа                                    | Статус                    |
| ----------------------------------------- | ------------------------- |
| Исправление английского                   | подготовка                |
| Приведение обозначений в порядок          | подготовка                |
| Литературный обзор                        | подготовка                |
| Освоение `s2fft`                          | подготовка                |
| Освоение SVD/HealPIX tooling              | подготовка                |
| Написание библиотечного кода              | подготовка                |
| Численный unit test известных формул      | подготовка / verification |
| Доказательство нового potential formalism | **наука**                 |
| Вывод response через потенциал            | **наука**                 |
| Аналитический beam                        | **наука**                 |
| Delta-function limit                      | **наука**                 |
| Критерий малости pulsar term              | **наука**                 |
| Arbitrary-polarization source             | **наука**                 |
| Potential reconstruction                  | **наука**                 |
| Stochastic covariance formulation         | **наука**                 |
| End-to-end injection/recovery             | **наука**                 |
| Joint deterministic + stochastic search   | **наука**                 |
| Оптимизация PTA                           | **наука**                 |
| Galactic eigenbasis                       | **наука**                 |

Но есть тонкая граница:

> написание кода — подготовка; разработка **алгоритма**, который этим кодом реализуется, — научный результат.

Например, SVD — сам по себе не научная новизна.

Но:

$$
R=U\Sigma V^\dagger
$$

и обнаружение того, что **из ($N$) пульсаров практически измеримы только ($N_{\rm eff}$) потенциальных мод, причём структура ($U$) имеет физическую интерпретацию**, — уже вполне научный результат.

---

# 3. Как я бы разделила это на статьи

Вот здесь я бы **не делала три равноправные статьи с самого начала**.

Я вижу:

$$
\boxed{\text{Article I}}
$$

$$
\boxed{\text{Article II}}
$$

$$
\boxed{\text{Article III — опционально}}
$$

---

# Статья I — основная методологическая

### Возможное название

**“A potential-based description of pulsar timing array responses to gravitational waves”**

или чуть более физически:

**“Spin potentials and the spatial response of pulsar timing arrays to gravitational waves”**

Это должна быть **главная работа**, на которой стоит весь проект.

### В неё я включила бы

### 1. Motivation

Что известно:

$$
\text{GW field}
\rightarrow
E/B
$$

и PTA insensitive to tensor curl. Это принципиальный результат Gair et al. ([arXiv][1])

Но существующий description рассматривает прежде всего harmonic coefficients.

Вы предлагаете:

$$
\boxed{
E/B
\rightarrow
\text{spin potentials}
\rightarrow
\text{local PTA response}
}
$$

### 2. Formalism

Полный вывод ваших потенциалов.

### 3. PTA response

$$
r_a=\mathcal R_a[\Phi_E].
$$

### 4. Ideal antenna

$$
r_a\propto\Phi_E(\hat p_a).
$$

### 5. Real antenna

$$
r_a=
\int K_a\Phi_E,d\Omega.
$$

### 6. Analytical beam

Это обязательный результат.

### 7. Delta limit

Показать, почему response можно интерпретировать как локальную выборку потенциала.

### 8. Numerical validation

Показать:

$$
r_{\rm polarization}=r_{\rm potential}.
$$

### 9. Finite PTA

Показать влияние конечного числа мультиполей и pulsar term.

### 10. Reconstruction

Минимальная демонстрация:

$$
r_a\rightarrow\hat\Phi_E.
$$

Я бы даже **не стала в первой версии статьи включать весь realistic red-noise Bayesian machinery**.

---

# Что станет главным claim статьи I

Не:

> «мы разработали ещё один способ разложения GW background».

А:

> **PTA response admits a potential formulation in which the timing residual associated with the observable gradient sector can be interpreted as a local measurement of a spin potential, with the finite pulsar response represented as a spatial kernel; this provides a natural map-making and reconstruction framework.**

Вот это действительно выглядит как самостоятельная методологическая идея.

---

# Статья II — восстановление и поиск источников

После первой статьи уже не надо заново доказывать ваш formalism.

### Возможное название

**“Reconstructing gravitational-wave potentials and localized sources with pulsar timing arrays”**

Сюда идут:

### 1. Eigenmode reconstruction

$$
R=U\Sigma V^\dagger.
$$

### 2. Number of measurable modes

$$
N_{\rm eff}.
$$

### 3. Potential map reconstruction

$$
d\rightarrow\hat\Phi.
$$

### 4. Deterministic source

$$
\hat\Phi
\rightarrow
(\theta_s,\phi_s).
$$

### 5. Arbitrary polarization

### 6. Stochastic background

$$
C_{\ell m,\ell'm'}
$$

и 3j formalism.

### 7. Noise

$$
N_{\rm white}+N_{\rm red}.
$$

### 8. Joint stochastic + deterministic search

### 9. Injection/recovery

То есть логика второй статьи:

$$
\boxed{
\text{formalism}
\rightarrow
\text{inverse problem}
\rightarrow
\text{actual inference}
}
$$

Это очень естественная вторая работа.

---

# Статья III — оптимальный PTA

Вот здесь я бы уже дала отдельное название, например:

**“Optimizing pulsar timing array geometry for anisotropic gravitational-wave sources”**

или:

**“Optimal pulsar configurations for gravitational-wave mapping with pulsar timing arrays”**

### Содержание

$$
\text{Fisher / likelihood}
\rightarrow
\text{optimization}
\rightarrow
\text{source morphology}
\rightarrow
\text{optimal PTA}.
$$

Три основных случая:

### Isotropic

$$
I(\hat n)=\mathrm{const}
$$

→ uniform configuration.

### Point source

$$
I(\hat n)=\delta(\hat n-\hat n_s)
$$

→ localization around source.

### Cosmic-string-like source

$$
I(\hat n)=I_{\rm CS}(\hat n)
$$

→ extended optimized configuration.

А затем:

$$
\text{idealized}
\rightarrow
\text{real Galactic population}.
$$

И именно здесь я бы поместила **Galactic basis**.

---

# 4. Что я бы вообще пока убрала

Есть несколько вещей, которые выглядят очень интересно, но опасны как источники scope creep.

## 4.1. «Давайте сразу сделаем идеальный universal Galactic basis»

Нет.

Сначала нужно показать:

$$
\text{SVD of PTA response}
$$

для простой геометрии.

После этого спросить:

> Какие собственные функции возникают у реального PTA?

Только тогда имеет смысл Galactic basis.

---

## 4.2. Космическую струну не надо сейчас делать центральным объектом

Это особенно важно.

Ваш исходный проект начинался с cosmic strings, но нынешняя идея намного фундаментальнее:

$$
\boxed{
\text{What does a PTA actually measure?}
}
$$

Если сделать cosmic string центром первой статьи, вы искусственно сузите фундаментальный результат до одного класса источников.

Лучше наоборот:

$$
\text{general formalism}
$$

→

$$
\text{generic localized source}
$$

→

$$
\text{cosmic string as application}.
$$

---

# 5. Самые трудные этапы

Я бы оценивала трудоёмкость примерно так.

| Этап                     | Трудоёмкость |  Риск |
| ------------------------ | -----------: | ----: |
| Формализм                |         ★★★★ |   ★★★ |
| Численная verification   |           ★★ |    ★★ |
| Аналитический beam       |        ★★★★★ |  ★★★★ |
| Pulsar-term criterion    |          ★★★ |   ★★★ |
| Arbitrary source         |           ★★ |    ★★ |
| Potential reconstruction |         ★★★★ |  ★★★★ |
| Stochastic background    |        ★★★★★ |  ★★★★ |
| Realistic noise          |         ★★★★ |   ★★★ |
| Joint search             |        ★★★★★ |  ★★★★ |
| PTA optimization         |        ★★★★★ | ★★★★★ |
| Galactic basis           |        ★★★★★ | ★★★★★ |

И отсюда вытекает важный практический вывод:

$$
\boxed{
\text{первая статья не должна зависеть от последних четырёх строк}
}
$$

---

# 6. Что является критическим для публикации первой статьи

Я бы поставила четыре **publication gates**.

### Gate A — mathematical closure

Нужно иметь строгую цепочку

$$
h_{+,\times}
\leftrightarrow
\Phi_{G,C}
\leftrightarrow
r_a.
$$

Без неё статья не готова.

### Gate B — numerical equivalence

Нужно доказать численно:

$$
r_{\rm conventional}
\simeq
r_{\rm potential}.
$$

Без этого новый formalism остаётся формальной перестановкой обозначений.

### Gate C — physical interpretation

Нужен beam:

$$
K(\theta)
$$

и доказательство его квазилокального характера.

**Вот это, на мой взгляд, главный потенциально оригинальный результат.**

### Gate D — useful consequence

Нужно показать, что новый язык **что-то позволяет сделать лучше**:

$$
r_a
\rightarrow
\Phi_G
\rightarrow
\text{map}.
$$

То есть хотя бы один end-to-end reconstruction experiment.

---

# 7. Очень важная концептуальная проверка

Я бы ещё до написания статьи сформулировала буквально один вопрос:

> **Что является новым объектом, которого не было у Gair et al.?**

Потому что Gair et al. уже дали:

$$
\text{GW angular field}
\rightarrow
\text{gradient/curl harmonics}
\rightarrow
\text{PTA response}.
$$

И они уже показали, что curl sector PTA не видит. ([arXiv][1])

Поэтому ваша новизна должна быть не сформулирована как:

> «мы используем gradient modes».

А примерно как:

$$
\boxed{
\text{gradient harmonic coefficients}
\rightarrow
\text{potential field}
}
$$

и далее

$$
\boxed{
\text{PTA response}=
\text{sampling/convolution of this potential}
}
$$

с физической интерпретацией соответствующего kernel.

И если аналитический ($K$) действительно позволяет установить переход

$$
K\rightarrow\delta
$$

и физически интерпретировать eigenmodes как пространственные измерительные моды потенциала, то именно здесь, на мой взгляд, находится **сердце новизны**.

---

# 8. Ещё один момент, который я бы обязательно добавила

Я бы не ограничивалась демонстрацией того, что

$$
r_a \propto \Phi(\hat p_a).
$$

Нужно выяснить, **насколько это утверждение фундаментально и насколько оно зависит от выбранного truncation по ($\ell$)**.

Очень интересный эксперимент:

$$
L_{\max}=2,4,6,\ldots,100
$$

и для каждого (L_{\max}):

$$
K_{L_{\max}}(\theta).
$$

После чего определить angular width:

$$
\Delta\theta(L_{\max})
$$

и боковые лепестки.

Тогда можно получить настоящий результат масштаба:

$$
\boxed{
\Delta\theta
\sim
f(L_{\max})
}
$$

и, возможно, связать пространственное разрешение потенциальной карты с числом доступных мультиполей.

А дальше уже появляется естественная связь:

$$
\boxed{
N_{\rm pulsar}
\rightarrow
L_{\rm measurable}
\rightarrow
\Delta\theta_{\rm PTA}.
}
$$

Вот это может стать очень сильным мостом от теории к оптимизации.

---

# 9. Итоговая последовательность, которую я бы реально поставила себе в рабочий план

Если превратить всё выше в буквально линейный roadmap:

### Подготовка

**0. Инфраструктура**
$$
\text{notation + literature + s2fft + SVD + tests}
$$

↓

### Первая статья

**1. Закрыть аналитический potential formalism**

↓

**2. Проверить polarization ↔ potential численно**

↓

**3. Получить аналитический PTA kernel**

↓

**4. Исследовать ($K\rightarrow\delta$)**

↓

**5. Получить критерий малости pulsar term**

↓

**6. Проверить arbitrary-polarization source**

↓

**7. Сделать минимальную potential reconstruction**

↓

$$
\boxed{\text{ARTICLE I}}
$$

↓

### Вторая статья

**8. Полноценное eigenmode/SVD reconstruction**

↓

**9. Deterministic source recovery**

↓

**10. Stochastic anisotropic background**

↓

**11. White + red noise**

↓

**12. Joint deterministic + stochastic inference**

↓

$$
\boxed{\text{ARTICLE II}}
$$

↓

### Третья статья

**13. Определить objective function для optimal PTA**

↓

**14. Изотропный фон**

↓

**15. Точечный источник**

↓

**16. Cosmic-string morphology**

↓

**17. Реальное распределение Galactic pulsars**

↓

**18. Galactic eigenbasis**

↓

$$
\boxed{\text{ARTICLE III}}
$$

---

# И главное: что я бы делала прямо сейчас

Не начинала бы с cosmic strings, realistic noise или оптимизации.

**Первый научный milestone я бы сформулировала всего в одной фразе:**

$$
\boxed{
\text{Derive and numerically validate the PTA response kernel for the gradient potential, and establish its localized/delta-function limit.}
}
$$

То есть последовательность ближайшей работы:

$$
\boxed{
\text{potential formalism}
\rightarrow
\text{polarization/potential equivalence}
\rightarrow
K(\theta)
\rightarrow
K_\ell
\rightarrow
\delta\text{-limit}
\rightarrow
\text{finite pulsar term}.
}
$$

После этого у вас уже будет не просто «набор вычислений к будущей статье», а **готовое научное ядро первой публикации**.

И я бы именно так объясняла структуру проекта руководительнице: не «у меня ещё 18 пунктов», а **«есть основная методологическая работа, а затем две независимые линии приложений — reconstruction и PTA optimization»**. Это существенно проще защищать и организационно, и научно.

[1]: https://arxiv.org/abs/1406.4664?utm_source=chatgpt.com "Mapping gravitational-wave backgrounds using methods from CMB analysis: Application to pulsar timing arrays"
[2]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.111.083039?utm_source=chatgpt.com "Pulsar timing array sensitivity to anisotropies in the gravitational wave background | Phys. Rev. D"
[3]: https://academic.oup.com/mnras/article/423/3/2642/2460717?utm_source=chatgpt.com "optimal schedule for pulsar timing array observations | Monthly Notices of the Royal Astronomical Society | Oxford Academic"


Да. С учётом именно этого уточнения я **меняю свою предыдущую рекомендацию**: я уже не считаю, что eigenmodes нужно обязательно выносить во вторую статью. Более того, теперь мне кажется, что есть сильный аргумент сделать **собственные моды PTA центральной частью первой статьи**, а beam — их естественным следствием и физической иллюстрацией.

Я дополнительно проверила свежую работу Andrianov & Kopeikin: она действительно вышла на arXiv 21 июля 2026 года и прямо формулирует полную Earth–pulsar response как линейную inverse problem; в abstract они уже выделяют четыре режима, cutoff (l_{\rm cut}\sim\omega L), предел разрешения и SVD/Fisher analysis. То есть приписывать вашей работе самостоятельное открытие этих режимов или fundamental resolution limits было бы неправильно. ([arXiv][1])

И вот именно поэтому ваша поправка очень существенна.

---

# 1. Я бы теперь разделила три вещи

Здесь есть очень полезное различие:

### A. Что делает физическая антенна?

Это

[
\Phi_G(\hat{\mathbf n})
\longrightarrow
r_a.
]

Полная антенна имеет некоторую angular response.

### B. Какие режимы антенна физически пропускает?

Это свойства

[
F_\ell(\omega L),
]

переходы Earth-term → transition → pulsar-term → cutoff и т. д.

**Это в значительной степени уже изучено Andrianov & Kopeikin.** Их текущая работа прямо посвящена этим режимам и фундаментальному angular resolution. ([arXiv][1])

### C. Какие **собственные функции измерения** имеет PTA?

Вот здесь появляется совершенно другой вопрос:

[
\boxed{
R:\mathcal H_{\rm sky}\rightarrow\mathbb R^{N_{\rm pulsar}}
}
]

и

[
R=U\Sigma V^\dagger.
]

Что означают (V_i(\hat n))?

В обычном harmonic language это просто линейные комбинации (Y_{\ell m}).

А в вашем formalism это можно интерпретировать как:

[
\boxed{
V_i(\hat n)
===========

\text{eigenmode of the gradient potential seen by the PTA}.
}
]

Вот это уже действительно выглядит как **новый концептуальный объект**.

И это существенно интереснее, чем просто ещё раз показать, что (F_\ell) переходит от одного режима к другому.

---

# 2. Более того: я бы даже переименовала смысл SVD

Раньше я говорила примерно:

> «Сделайте SVD и посмотрите, сколько гармоник PTA измеряет».

Теперь я бы не так это формулировала.

Ваш центральный вопрос можно поставить гораздо глубже:

> **What are the natural modes in which a pulsar timing array measures the observable gravitational-wave field?**

И ответ:

[
\boxed{
\text{They are eigenmodes of the PTA response operator acting on the gradient potential.}
}
]

То есть SVD становится не просто численной техникой реконструкции.

Он становится **частью физического формализма**.

Это принципиальная разница.

---

# 3. Почему здесь потенциал действительно меняет интерпретацию

У Gair et al. основная observable decomposition имеет вид

[
a_{\ell m}^{G}.
]

И Earth-term response можно записать как

[
R_{\ell m}^{G}\propto
N_\ell Y_{\ell m}(\hat p).
]

Это известно: PTA чувствителен к gradient modes и не чувствителен к curl modes. ([Caltech Authors][2])

В potential formulation же

[
a_{\ell m}^{G}=N_\ell \Phi_{\ell m}
]

(с точностью до ваших соглашений о нормировке), и тогда Earth-term response превращается в

[
r_a
\propto
\sum_{\ell m}
N_\ell \Phi_{\ell m}
Y_{\ell m}(\hat p_a).
]

А если выбрать именно такое определение потенциала, чтобы (N_\ell) был поглощён в оператор связи, возникает буквально

[
\boxed{
r_a\propto\Phi(\hat p_a).
}
]

И вот это уже не просто другое обозначение (a_{\ell m}^{G}).

Физическая интерпретация меняется:

> PTA не просто «измеряет gradient spherical-harmonic coefficients».

Он в соответствующем пределе **измеряет значения наблюдаемого поля потенциала в направлениях пульсаров**.

Это гораздо более наглядная картина.

---

# 4. Тогда beam становится не главным результатом, а демонстрацией этой интерпретации

Вот здесь я бы теперь полностью изменила акцент.

Не:

> «Мы подробно изучаем beam и обнаруживаем четыре режима».

Потому что это действительно уже занято Andrianov & Kopeikin. Их статья сама говорит, что они исследуют Earth-term dominated, transition, pulsar-term dominated и exponential cutoff regimes, а также angular resolution. ([arXiv][1])

А:

> «В нашем representation PTA response имеет вид spatial convolution of the gradient potential with an explicit angular kernel.»

И тогда:

[
r_a =
\int d\Omega,
K_a(\hat n)\Phi_G(\hat n).
]

В Earth-term approximation:

[
K_a(\hat n)
\approx
\delta(\hat n-\hat p_a).
]

То есть **ваша работа даёт физическую картину того, что означают уже известные (F_\ell)**.

Это очень важное различие:

[
\text{Andrianov:}
\quad
F_\ell
\longrightarrow
\text{regimes and resolution}
]

против

[
\text{ваша работа:}
\quad
F_\ell
\longrightarrow
K(\hat n,\hat p)
\longrightarrow
\text{spatial measurement of }\Phi_G.
]

---

# 5. Более того, я бы сделала очень красивый мост между beam и eigenmodes

Это, возможно, вообще станет самым сильным местом первой статьи.

Рассмотрим оператор

[
R:\Phi(\hat n)\mapsto
{r_a}.
]

Его kernel в общем случае:

[
R_a(\hat n)=K_a(\hat n).
]

Тогда SVD:

[
R=U\Sigma V^\dagger
]

даёт

[
\Phi(\hat n)
============

\sum_i c_i V_i(\hat n).
]

И теперь можно сказать:

[
\boxed{
V_i(\hat n)
\text{ — естественные пространственные моды градиентного потенциала, измеряемые PTA.}
}
]

Причём corresponding singular value

[
\sigma_i
]

определяет чувствительность PTA к данной моде.

Получается очень естественная физическая картина:

[
\boxed{
\Phi_G(\hat n)
==============

\sum_i c_i V_i(\hat n)
}
]

и

[
r_a
===

\sum_i
\sigma_i U_{ai}c_i.
]

Это уже не просто reconstruction algorithm.

Это **modal description самого измерения**.

---

# 6. И тут появляется ещё одна интересная связь с Gair et al.

У Gair есть maximum-likelihood reconstruction sky map, и современный обзор также описывает PTA reconstruction как задачу, в которой insensitive modes фактически отбрасываются. ([SpringerLink][3])

Но у вас можно сделать следующий conceptual step:

обычная карта:

[
\hat\Phi(\hat n)
]

не является первичным объектом измерения.

Первичными являются:

[
\boxed{
c_i=
\langle \Phi,V_i\rangle
}
]

для observable eigenmodes (V_i).

То есть вместо:

> «Мы пытаемся восстановить бесконечномерное sky field с помощью конечного PTA»

получаем:

> «PTA имеет конечномерное пространство измеримых мод потенциального поля; reconstruction происходит внутри этого observable subspace.»

Это, на мой взгляд, **очень сильная формулировка**.

---

# 7. Поэтому я бы теперь перестроила первую статью

Не так:

> Formalism → beam → beam regimes → reconstruction.

А так:

## I. GW field and gradient potential

[
h_{+,\times}
\rightarrow
a_{\ell m}^{G,C}
\rightarrow
\Phi_G,\Phi_C.
]

Здесь устанавливается ваш formalism.

---

## II. PTA response in potential representation

Получаем

[
r_a=\mathcal R_a[\Phi_G].
]

Здесь же показываем:

[
\boxed{\text{Earth term}}
\Rightarrow
r_a\propto\Phi_G(\hat p_a).
]

Это главный physical interpretation result.

---

## III. Full Earth–pulsar response

Получаем

[
r_a=
\int K_a(\hat n)\Phi_G(\hat n)d\Omega.
]

И только здесь появляется beam.

Причём нужно **очень явно** сказать:

> The spectral properties of the full response and its resolution regimes have been studied previously by Andrianov & Kopeikin. Here we reinterpret the same response as an angular kernel acting on the gradient potential.

Таким образом вы сразу снимаете вопрос с novelty.

---

## IV. Potential-space eigenmodes

И вот здесь я бы сделала **центральную новую секцию**.

[
R=U\Sigma V^\dagger.
]

Показать:

[
V_i(\hat n)
]

для нескольких (i).

И самое главное — сравнить их с harmonic modes.

То есть:

[
V_i
===

\sum_{\ell m}
v_{i,\ell m}Y_{\ell m}.
]

Тогда можно посмотреть:

[
|v_{i,\ell m}|^2
]

и понять, какие (\ell)-режимы входят в каждую собственную функцию.

Это непосредственно связывает:

[
\text{Gair harmonic picture}
\leftrightarrow
\text{your potential picture}
\leftrightarrow
\text{actual PTA eigenmodes}.
]

---

# 8. И только потом beam analysis

Beam здесь становится буквально **sanity check / physical interpretation** eigenmode structure.

Можно показать:

### Earth-term limit

[
K(\theta)
\simeq
\delta(\theta).
]

### Full response

[
K(\theta;\omega L)
]

получает конечную ширину и oscillatory structure.

### Harmonic representation

[
K(\theta)
=========

\sum_\ell
K_\ell P_\ell(\cos\theta)
]

и

[
K_\ell
\leftrightarrow F_\ell N_\ell.
]

И здесь уже дать ссылку на Andrianov & Kopeikin:

> наши (K_\ell) воспроизводят известную спектральную структуру полной PTA response.

Но **не делать из воспроизведения этих переходов отдельный claim**.

---

# 9. А что тогда является вашей потенциальной новой физикой?

Я вижу минимум три результата.

### Результат №1 — potential representation

[
\boxed{
r_a
===

\mathcal R_a[\Phi_G].
}
]

И в Earth-term limit:

[
\boxed{
r_a\propto\Phi_G(\hat p_a).
}
]

### Результат №2 — spatial kernel interpretation

[
\boxed{
\mathcal R_a
\equiv
\text{spatial convolution operator on }\Phi_G.
}
]

Это связывает harmonic response с реальным пространственным измерением.

### Результат №3 — PTA eigenmodes in potential space

[
\boxed{
R=U\Sigma V^\dagger,
\qquad
V_i(\hat n)=\text{observable potential modes}.
}
]

И вот **№3 я теперь считаю главным новым результатом**.

---

# 10. А восстановление карты тогда не надо откладывать во вторую статью

Это следствие №3.

Если

[
\Phi_G
======

\sum_i c_iV_i
]

и

[
d_i
===

\sigma_i c_i+n_i,
]

то reconstruction становится естественной:

[
\hat c_i
========

\frac{d_i}{\sigma_i}
]

в простом noiseless/regularized случае.

То есть:

[
\boxed{
\text{eigenmode decomposition}
\rightarrow
\text{potential reconstruction}
}
]

— фактически **одна научная история**.

Разрывать её между двумя статьями теперь действительно не стоит.

---

# 11. Тогда я бы существенно изменила границу статей

## Статья I — фундаментальная методологическая

### Рабочее название

**“A potential-space description of gravitational-wave measurements with pulsar timing arrays”**

или:

**“Gravitational-wave potential eigenmodes of pulsar timing arrays”**

Вторая формулировка мне даже нравится больше как title, если eigenmodes действительно получатся сильным результатом.

### Содержание

[
\boxed{
\begin{array}{c}
GW\ polarization\
\downarrow\
gradient/curl\
\downarrow\
spin\ potentials\
\downarrow\
PTA\ response\
\downarrow\
potential\ kernel\
\downarrow\
PTA\ eigenmodes\
\downarrow\
potential\ reconstruction
\end{array}
}
]

То есть статьи I **должно хватить** на:

* новый formalism;
* analytic response;
* conventional ↔ potential equivalence;
* Earth-term/local limit;
* full response kernel;
* связь kernel с известными (F_\ell);
* SVD/eigenmodes;
* reconstruction;
* numerical validation.

Именно здесь beam находится **не в центре, а как bridge between old harmonic description and your new potential-space interpretation**.

---

# 12. Статья II тогда становится существенно более физической/statistical

### “Inference of gravitational-wave sources in the potential eigenmode representation of pulsar timing arrays”

И сюда:

[
\Phi_G
\rightarrow
{c_i}
\rightarrow
\text{inference}.
]

Содержание:

* arbitrary-polarization deterministic source;
* localized source;
* stochastic background;
* anisotropic background;
* 3j formalism;
* white/red noise;
* injection/recovery;
* joint CW + SGWB search.

То есть первая статья отвечает:

> **What does a PTA measure?**

Вторая:

> **What can we infer once we formulate the measurement in those natural modes?**

Это очень чистое разделение.

---

# 13. Статья III — оптимизация

А уже здесь:

[
\text{source morphology}
\rightarrow
\text{optimal eigenmode spectrum}
\rightarrow
\text{optimal pulsar geometry}.
]

То есть оптимизация становится естественным развитием первой статьи.

Особенно интересная постановка:

для конфигурации PTA (P)

[
R(P)=U(P)\Sigma(P)V^\dagger(P),
]

и оптимизировать не просто какой-нибудь SNR, а свойства observable subspace:

[
{\sigma_i}.
]

Например:

[
\max \sigma_{\min},
\qquad
\max\sum_i\sigma_i^2,
]

или quantity, связанный с конкретным source model.

Это уже выглядит гораздо органичнее, чем отдельная абстрактная задача размещения пульсаров.

---

# 14. Что делать с beam-пунктом практически

Я бы **не убирала его вообще**. Наоборот, я бы оставила его достаточно подробно — но поменяла статус.

Не:

> «Исследуем новые режимы PTA beam».

А:

> «Derive the spatial response kernel corresponding to the previously known harmonic response and use it to interpret the PTA measurement in potential space.»

И тогда нужны всего четыре ключевых результата:

[
\boxed{
K_\ell \leftrightarrow F_\ell N_\ell
}
]

[
\boxed{
K(\theta)
}
]

[
\boxed{
K\rightarrow\delta
\quad\text{in Earth-term regime}
}
]

[
\boxed{
K\text{ becomes nonlocal when pulsar term matters}.
}
]

А уже графики бабочек, critical (l), переходы и т. п. — **коротко и с явной атрибуцией к предыдущей работе**.

То есть вы не говорите:

> «Мы открыли бабочку».

Вы говорите:

> «Известная в harmonic representation структура response приобретает следующую пространственную интерпретацию в potential representation».

И это совершенно другой scientific claim.

---

# 15. И есть ещё одна вещь, которую теперь я бы поставила в самый центр

Я бы специально сделала **сравнение трёх базисов**:

[
\boxed{
Y_{\ell m}
\quad\leftrightarrow\quad
\text{GW gradient modes}
\quad\leftrightarrow\quad
V_i^{\rm PTA}.
}
]

То есть показать:

### Harmonic basis

[
Y_{\ell m}.
]

### Response-weighted harmonic basis

[
N_\ell Y_{\ell m}.
]

### PTA eigenbasis

[
V_i.
]

И спросить:

> Чем отличается «разложение поля по гармоникам» от «разложения измеряемой части поля по собственным модам прибора»?

Вот это уже очень глубокий методологический вопрос.

Потому что

[
Y_{\ell m}
]

определяются геометрией сферы и **не зависят от конкретного PTA**, тогда как

[
V_i^{\rm PTA}
]

определяются оператором

[
R[{\hat p_a,L_a,\sigma_a,\ldots}].
]

То есть собственные моды естественно становятся **detector-adapted basis**.

И это очень хорошо связывается с вашей дальнейшей идеей «галактического базиса»:

[
V_i^{\rm ideal}
\quad\longrightarrow\quad
V_i^{\rm real\ PTA}
\quad\longrightarrow\quad
V_i^{\rm Galactic}.
]

Тогда Galactic basis уже не выглядит внезапной дополнительной темой. Он становится естественным продолжением общей концепции.

---

# 16. Поэтому мой пересмотренный план выглядит так

### Первая статья

[
\boxed{
\textbf{What does a PTA measure?}
}
]

[
\text{GW}
\rightarrow
\text{gradient potential}
\rightarrow
\text{response operator}
\rightarrow
\text{spatial kernel}
\rightarrow
\textbf{potential eigenmodes}
\rightarrow
\text{reconstruction}.
]

### Вторая

[
\boxed{
\textbf{What can we infer from those measurements?}
}
]

[
\text{eigenmode amplitudes}
\rightarrow
\text{deterministic sources}
+
\text{stochastic backgrounds}
+
\text{noise}.
]

### Третья

[
\boxed{
\textbf{How should we build the PTA to measure them?}
}
]

[
\text{source}
\rightarrow
\text{desired modes}
\rightarrow
\text{optimal geometry}.
]

---

## И я теперь действительно согласна с вашим ощущением

Ваше интуитивное «**собственные моды — это, кажется, и есть фундаментально новый результат**» мне после этого уточнения кажется гораздо более убедительным, чем моя первоначальная схема.

Причём я бы даже усилила его:

> **Новизна не просто в том, что собственные моды вычисляются через SVD. Новизна в том, что вы переопределяете объект реконструкции: вместо абстрактных (a_{\ell m}^{G}) вводите gradient potential как физическое поле, а SVD отклика превращается в построение естественного, detector-adapted базиса этого поля.**

Это вполне естественно соседствует с Gair et al., которые дают gradient/curl harmonic description и показывают blind spot PTA к curl, и с Andrianov & Kopeikin, которые исследуют spectral response полной Earth–pulsar антенны, но не делает вид, будто вы заново открываете их результаты. ([Caltech Authors][2])

И я бы поэтому **не делала beam крупным самостоятельным разделом первой статьи**. Я бы сделала его примерно одним из трёх элементов вокруг главного результата:

[
\boxed{
\text{potential}
\quad-\quad
\text{response kernel}
\quad-\quad
\text{eigenmodes}.
}
]

Причём центральным из них были бы именно **eigenmodes**.

И, кстати, наличие именно сейчас свежей статьи Andrianov & Kopeikin делает такую расстановку акцентов ещё важнее: поскольку их работа от 21 июля 2026 года уже прямо заявляет closed-form response, four regimes, angular cutoff, Fisher и SVD, вашу статью нужно позиционировать рядом с ней не как ещё один analysis of PTA resolution, а как **другой representation и новую интерпретацию измеряемого пространства**. ([arXiv][1])

Это, на мой взгляд, существенно более безопасная и одновременно более сильная позиция для публикации.

[1]: https://arxiv.org/abs/2607.19329 "[2607.19329] Gravitational-Wave Sky Mapping with Pulsar Timing Arrays: The Full Earth-Pulsar Response and Fundamental Resolution Limits"
[2]: https://authors.library.caltech.edu/records/d86bj-e7t49?utm_source=chatgpt.com "Mapping gravitational-wave backgrounds using methods from CMB analysis: Application to pulsar timing arrays | CaltechAUTHORS"
[3]: https://springerlink.fh-diploma.de/article/10.1007/s41114-017-0004-1?utm_source=chatgpt.com "Detection methods for stochastic gravitational-wave backgrounds: a unified treatment | Living Reviews in Relativity | Springer Nature Link"
