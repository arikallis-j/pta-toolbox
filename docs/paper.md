# paper

Я продолжаю исследовать гравитационно-волновой интерферометр (PTA) в новом формализме. Краткая сводка того, что я нашeл:

## Потенциалы гравитационного поля

Мы можем ввести spin-2 потенциалы grad и curl части нашего возмущения гравитационного поля как:

$$
\Psi_G(\Omega) = \sum_{lm} \psi^G_{lm} \cdot Y_{lm}(\Omega)
$$

$$
\Psi_C(\Omega) = \sum_{lm} \psi^C_{lm} \cdot Y_{lm}(\Omega)
$$

где коэффициенты равны:

$$
\psi^G_{lm} = \frac{N_l}{2} \cdot a^G_{lm}
$$

$$
\psi^C_{lm} = \frac{N_l}{2} \cdot a^C_{lm}
$$

где:

$$
N_l = \sqrt{\frac{2(l - 2)!}{(l+2)!}}
$$

$a^G_{lm}$, $a^C_{lm}$ - коэффициенты разложения тензора $h_{ab}(f, \Omega)$ (по частоте), которые можно непосредственно связать с $+ / \times$ или $L/R$ поляризациями тензора $h_{ab}$. В частности, для $L/R$ поляризаций выражения выглядят так:

$$
h_L = \frac{1}{2} \cdot \sum_{lm} (a^G_{lm} + i \cdot a^C_{lm}) \cdot {{}_{+2}Y_{lm}}
$$

$$
h_R = \frac{1}{2} \cdot \sum_{lm} (a^G_{lm} - i \cdot a^C_{lm})  \cdot  {{}_{-2}Y_{lm}}
$$

где ${{}_{+2}Y_{lm}}$ и ${{}_{-2}Y_{lm}}$ - spin-2 сферические гармоники. А $+/\times$ поляризации можно выразить через $L/R$ как:

$$
h_+ = \frac{1}{\sqrt{2}} (h_L + h_R)
$$

$$
h_{\times} = \frac{1}{i\sqrt{2}} (h_L - h_R)
$$

В таком случае, используя связь ${{}_{+2}Y_{lm}}$ и ${{}_{-2}Y_{lm}}$ c $Y_{lm}$ через операторы $\eth$ и $\bar\eth$, мы можем получить выражения для $h_L$ и $h_R$ поляризаций через $\Psi_G$ и $\Psi_C$ как:

$$
h_L = \frac{1}{\sqrt{2}} \cdot \eth^2 [ \Psi_G + i \cdot \Psi_C]
$$

$$
h_R = \frac{1}{\sqrt{2}} \cdot \bar\eth^2 [ \Psi_G - i \cdot \Psi_C]
$$

И через связь $L/R$ и $+/\times$ поляризаций:

$$
h_+ = \frac{1}{2} \cdot \left[ (\eth^2 + \bar\eth^2)[\Psi_G] + i \cdot (\eth^2 - \bar\eth^2)[\Psi_C] \right]
$$

$$
h_{\times} = \frac{1}{2i} \cdot \left[ (\eth^2 - \bar\eth^2)[\Psi_G] + i \cdot (\eth^2 + \bar\eth^2)[\Psi_C] \right]
$$

## Сигнал идеального PTA

Сигнал, принимаемый PTA (спектр остаточных уклонений) $r(f)$ на тензор $h_{ab}(f, \Omega)$ можно выразить как:

$$
r(f) = \sum_{lm} \sum_{P} R^P_{lm}(f) \cdot a^P_{lm}(f)
$$

где $P = \{G, C\}$, а отклик антенны $R^P_{lm}(f)$ равен:

$$
R^P_{lm}(f) = (i 2 \pi f)^{-1} \cdot F^P_{lm}(f)
$$

Обозначим первое слагаемое, связанное только с тем, что мы измеряем не красное смещение, а именно остаточные уклонения, как:

$$
\mathcal{A}(f) = \frac{1}{i 2 \pi f}
$$

И рассмотрим чисто геометрический отклик антенны $F^P_{lm}(f)$. В short-wave appoximation, он будет равен (Gair et al. 2014):

$$
F^G_{lm}(f) \approx 4 \pi \cdot \frac{N_l}{2} \cdot (-1)^l \cdot Y_{lm}(\hat p)
$$

$$
F^C_{lm}(f) \approx 0
$$

где $\hat p$ - единичный вектор направления на пульсар. Как мы видим, отклик на curl часть поля равен нулю, поэтому остановимся только на grad части. Подставим в наше исходное выражение для $r(f)$ и получаем:

$$
r(f) = \mathcal{A}(f) \cdot \sum_{lm} F^G_{lm}(f) \cdot a^G_{lm}(f)
$$

Подставим выражение для $F^G_{lm}(f)$ и сгруппируем слагаемые:

$$
r(f) = 4 \pi \cdot  \mathcal{A}(f) \cdot \sum_{lm} \left[\frac{N_l}{2} \cdot a^G_{lm}(f) \right] \cdot \left[(-1)^l \cdot Y_{lm}(\hat p) \right]
$$

Как можно видеть, второе слагаемое в сумме сворачивается как:

$$
(-1)^l \cdot Y_{lm}(\hat p) = Y_{lm}(- \hat p)
$$

а первое по нашему определению равно $\psi^G_{lm}$. Итого, получаем:

$$
r(f) = 4 \pi \cdot \mathcal{A}(f) \cdot \sum_{lm} \psi^G_{lm} \cdot Y_{lm}(- \hat p) = \mathcal{A}(f) \cdot 4 \pi  \Psi_G(- \hat p)
$$

То есть, идеальный PTA измеряет наш spin-2 потенциал grad части в одной точке - противоположной направлению на пульсар.

## Сигнал реального PTA

Уточним наши рассуждения вне short wave approximation, или же с учетом пульсарного члена. Согласно (Gair et al. 2015), формулы для $F^P_{lm}(f)$ будут равны:

$$
F^G_{lm}(f) = 4 \pi \cdot \frac{N_l}{2} \cdot (-1)^l \cdot Y_{lm}(\hat p) \cdot \mathcal{F}^G_l(2\pi f L/c)
$$

$$
F^C_{lm}(f) = 0
$$

где $L$ - расстояние до пульсара, $c$ - скорость света. То есть, отклик на curl часть поля всё ещё равен нулю, в то время как отклик на grad часть уточняется некоторым спектральным фильтром $\mathcal{F}^G_l$, который равен:

$$
\mathcal{F}^G_l(y) = \frac{i^l e^{-iy}}{2} \cdot \left[ (2 - 2iy + y^2) j_l(y) - i (6 + 4iy + y^2) \frac{\mathrm{d} j_l}{\mathrm{d} y} - (6iy - y^2) \frac{\mathrm{d}^2 j_l}{\mathrm{d} y^2} - i y^2  \frac{\mathrm{d}^3 j_l}{\mathrm{d} y^3} \right]
$$

где $j_l(y)$ - сферические функции Бесселя, а $y = 2\pi f L/c$.

С учетом $\mathcal{F}^G_l$, формула для сигнала PTA перепишется как:

$$
r(f) = 4 \pi \cdot \mathcal{A}(f) \cdot \sum_{lm} \mathcal{F}^G_l \cdot \psi^G_{lm} \cdot Y_{lm}(- \hat p) = 4 \pi \cdot \mathcal{A}(f) \cdot \sum_{l} \mathcal{F}^G_l \cdot \Psi^G_{l}(- \hat p)
$$

где:
$$
\Psi^G_{l}(\Omega) = \sum_m \psi^G_{lm} \cdot Y_{lm}(\Omega)
$$

Введем оператор проекции как:

$$
\mathcal{P}_l(\Omega, \Omega') = \sum_m Y_{lm}(\Omega) \cdot Y_{lm}^* (\Omega')
$$

где сумма по $m$ равна:

$$
\mathcal{P}_l(\Omega, \Omega') = \frac{2l + 1}{4\pi} P_l(\Omega \cdot \Omega')
$$

где $P_l(\Omega \cdot \Omega')$ - полиномы Лежандра. Тогда наша функция $\Psi^G_l$ перепишется в интегральном виде как:

$$
\Psi^G_l(\Omega) = \int_{S^2} \mathrm{d \Omega} \cdot \mathcal{P}_l(\Omega, \Omega') \cdot \Psi_G(\Omega')
$$

А наша сумма по $l$ будет выглядеть как:

$$
\sum_{l} \mathcal{F}^G_l \cdot \Psi^G_{l}(\Omega) = \int_{S^2} \mathrm{d \Omega} \cdot \left[ \sum_l \mathcal{F}^G_l \cdot \mathcal{P}_l(\Omega, \Omega') \right] \cdot \Psi_G(\Omega')
$$

То есть, получается обычная свертка функции $\Psi_G$ c интегральным ядром $K^G(\Omega, \Omega')$:

$$
\sum_{l} \mathcal{F}^G_l \cdot \Psi^G_{l}(\Omega) = \int_{S^2} \mathrm{d \Omega} \cdot K^G(\Omega, \Omega') \cdot \Psi_G(\Omega')
$$

где интегральное ядро равно:

$$
K^G(\Omega, \Omega') = \sum_l \frac{2l + 1}{4\pi} \cdot \mathcal{F}^G_l \cdot  P_l(\Omega \cdot \Omega')
$$

Обозначим операцию свертки функции с нашим спектральным фильтром как:

$$
\mathcal{F}^G[\Psi_G](\Omega) = \int_{S^2} \mathrm{d \Omega} \cdot K^G(\Omega, \Omega') \cdot \Psi_G(\Omega')
$$

Тогда наш итоговый отклик $r(f)$ будет равен:

$$
r(f) = 4 \pi \cdot \mathcal{A}(f) \cdot \mathcal{F}^G[\Psi_G](- \hat p)
$$

То есть, реальный PTA измеряет наш spin-2 свертку нашего потенциала grad части в одной точке, со спектральным фильтром, зависящем от комбинации $\omega L/c$, где $\omega = 2\pi f$ - циклическая частота.

```markdown
Привет! Я продолжаю исследовать гравитационно-волновой интерферометр (PTA) в новом формализме. Краткая сводка того, что я нашла.


## Потенциалы гравитационного поля
Мы можем ввести spin-2 потенциалы grad и curl части нашего возмущения гравитационного поля как:

$$ \Psi_G(\Omega) = \sum_{lm} \psi^G_{lm} \cdot Y_{lm}(\Omega) $$

$$ \Psi_C(\Omega) = \sum_{lm} \psi^C_{lm} \cdot Y_{lm}(\Omega) $$

где коэффициенты равны:

$$ \psi^G_{lm} = \frac{N_l}{2} \cdot a^G_{lm} $$

$$ \psi^C_{lm} = \frac{N_l}{2} \cdot a^C_{lm} $$

где:

$$ N_l = \sqrt{\frac{2(l - 2)!}{(l+2)!}} $$

$a^G_{lm}$, $a^C_{lm}$ - коэффициенты разложения тензора $h_{ab}(f, \Omega)$ (по частоте), которые можно непосредственно связать с $+ / \times$ или $L/R$ поляризациями тензора $h_{ab}$. В частности, для $L/R$ поляризаций выражения выглядят так:

$$ h_L = \frac{1}{2} \cdot \sum_{lm} (a^G_{lm} + i \cdot a^C_{lm}) \cdot {{}_{+2}Y_{lm}} $$

$$ h_R = \frac{1}{2} \cdot \sum_{lm} (a^G_{lm} - i \cdot a^C_{lm})  \cdot  {{}_{-2}Y_{lm}} $$

где ${{}_{+2}Y_{lm}}$ и ${{}_{-2}Y_{lm}}$ - spin-2 сферические гармоники. А $+/\times$ поляризации можно выразить через $L/R$ как:

$$ h_+ = \frac{1}{\sqrt{2}} (h_L + h_R) $$

$$ h_{\times} = \frac{1}{i\sqrt{2}} (h_L - h_R) $$

В таком случае, используя связь ${{}_{+2}Y_{lm}}$ и ${{}_{-2}Y_{lm}}$ c $Y_{lm}$ через операторы $\eth$ и $\bar\eth$, мы можем получить выражения для $h_L$ и $h_R$ поляризаций через $\Psi_G$ и $\Psi_C$ как:

$$ h_L = \frac{1}{\sqrt{2}} \cdot \eth^2 [ \Psi_G + i \cdot \Psi_C] $$

$$ h_R = \frac{1}{\sqrt{2}} \cdot \bar\eth^2 [ \Psi_G - i \cdot \Psi_C] $$

И через связь $L/R$ и $+/\times$ поляризаций:

$$ h_+ = \frac{1}{2} \cdot \left[ (\eth^2 + \bar\eth^2)[\Psi_G] + i \cdot (\eth^2 - \bar\eth^2)[\Psi_C] \right] $$

$$ h_{\times} = \frac{1}{2i} \cdot \left[ (\eth^2 - \bar\eth^2)[\Psi_G] + i \cdot (\eth^2 + \bar\eth^2)[\Psi_C] \right] $$


## Сигнал идеального PTA

Сигнал, принимаемый PTA (спектр остаточных уклонений) $r(f)$ на тензор $h_{ab}(f, \Omega)$ можно выразить как:

$$ r(f) = \sum_{lm} \sum_{P} R^P_{lm}(f) \cdot a^P_{lm}(f) $$

где $P = \{G, C\}$, а отклик антенны $R^P_{lm}(f)$ равен:

$$ R^P_{lm}(f) = (i 2 \pi f)^{-1} \cdot F^P_{lm}(f) $$

Обозначим первое слагаемое, связанное только с тем, что мы измеряем не красное смещение, а именно остаточные уклонения, как:

$$ \mathcal{A}(f) = \frac{1}{i 2 \pi f} $$ 

И рассмотрим чисто геометрический отклик антенны $F^P_{lm}(f)$. В short-wave appoximation, он будет равен (Gair et al. 2014):

$$ F^G_{lm}(f) \approx 4 \pi \cdot \frac{N_l}{2} \cdot (-1)^l \cdot Y_{lm}(\hat p) $$

$$ F^C_{lm}(f) \approx 0 $$ 

где $\hat p$ - единичный вектор направления на пульсар. Как мы видим, отклик на curl часть поля равен нулю, поэтому остановимся только на grad части. Подставим в наше исходное выражение для $r(f)$ и получаем:

$$ r(f) = \mathcal{A}(f) \cdot \sum_{lm} F^G_{lm}(f) \cdot a^G_{lm}(f) $$ 

Подставим выражение для $F^G_{lm}(f)$ и сгруппируем слагаемые: 

$$ r(f) = 4 \pi \cdot  \mathcal{A}(f) \cdot \sum_{lm} \left[\frac{N_l}{2} \cdot a^G_{lm}(f) \right] \cdot \left[(-1)^l \cdot Y_{lm}(\hat p) \right] $$

Как можно видеть, второе слагаемое в сумме сворачивается как:

$$ (-1)^l \cdot Y_{lm}(\hat p) = Y_{lm}(- \hat p) $$

а первое по нашему определению равно $\psi^G_{lm}$. Итого, получаем:

$$ r(f) = 4 \pi \cdot \mathcal{A}(f) \cdot \sum_{lm} \psi^G_{lm} \cdot Y_{lm}(- \hat p) = \mathcal{A}(f) \cdot 4 \pi  \Psi_G(- \hat p) $$

То есть, идеальный PTA измеряет наш spin-2 потенциал grad части в одной точке - противоположной направлению на пульсар.
```
