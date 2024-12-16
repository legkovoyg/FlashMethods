import math

import numpy as np
from base_eos_flash import EOSFlash

class PRFlash(EOSFlash):
    
    def _init_eos_parameters(self):
        """
        Инициализация параметров для уравнения состояния Пенг-Робинсона (PR).

        Возвращает:
            tuple:
                - a_i (ndarray): Параметр "a" для каждого компонента.
                - b_i (ndarray): Параметр "b" для каждого компонента.
                - BIPs (ndarray): Бинарные параметры взаимодействия (Binary Interaction Parameters).
                - psi_i (ndarray): Коэффициенты для расчёта температурной поправки альфа.
                - ac_i (ndarray): Значения параметра "a" при критических условиях.
                - tuple: Набор параметров (b_i, c_a_i) для дальнейших расчётов.

        Описание:
            - Параметры PR рассчитываются на основе критической температуры (Tkr),
            критического давления (Pkr) и ацентрического фактора (w).
            - Температурная поправка alpha_i корректирует параметр "a" в зависимости от температуры.
            - Матрица взаимодействий c_a_i учитывает бинарные взаимодействия между компонентами.
        """
        ac_i = 0.457235 * self._R**2 * self._Tkr**2.0 / self._Pkr
        psi_i = 0.37464 + 1.54226 * self._w - 0.26992 * self._w**2
        alpha_i = (1 + psi_i * (1 - (self._t / self._Tkr)**0.5))**2
        a_i = ac_i * alpha_i
        b_i = 0.077796 * self._R * self._Tkr / self._Pkr
        BIPs = self._c
        c_a_i = (1 - BIPs) * np.sqrt(a_i[:, None] * a_i[None, :]) 
        # для PR не нужен c_i, поэтому c_params = (b_i, c_a_i)

        return a_i, b_i, BIPs, psi_i, ac_i, (b_i, c_a_i)

    def _calculate_fugacity(self, P, T, molar_frac, b_i, c_a_i, isMax=True):
        """
        Рассчитывает фугитивность компонента в смеси на основе уравнения состояния PR (EOS).
        Параметры:
            P (float): Давление системы.
            T (float): Температура системы.
            molar_frac (ndarray): Мольные доли компонентов в смеси.
            b_i (ndarray): Параметры "b" для каждого компонента.
            c_a_i (ndarray): Кросс-параметры "a" для взаимодействия между компонентами.
            isMax (bool): Флаг для выбора корня уравнения состояния. 
                        True - используется максимальный корень, False - минимальный.
        Возвращает:
            tuple:
                - fz_i (ndarray): Фугитивность компонентов в смеси.
                - Z_v (float): Коэффициент сжимаемости смеси (Z-фактор).

        Описание:
            - Используется кубическое уравнение состояния для расчета Z-фактора.
            - Фугитивность определяется через логарифмическую зависимость от Z-фактора и параметров смеси.
            - Метод позволяет учитывать межмолекулярные взаимодействия и выбор корня для фазового состояния.
        """
        
        R = self._R
        aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
        bw = np.dot(molar_frac, b_i)
        Aw = aw * P / (R**2 * T**2)
        Bw = bw * P / (R * T)

        # SRK - Peneloux cubic equation coefficients
        a = 1
        b = -(1 - Bw)
        c = Aw - 3 * Bw**2 - 2 * Bw
        d = -(Aw * Bw - Bw**2 - Bw**3)

        roots = self._solve_cubic(a, b, c, d)
        Z_v = np.max(roots) if isMax else np.min(roots)

        avvv = np.dot(molar_frac, c_a_i)

        log_coeff = (
            np.log(molar_frac * P)
            - np.log(Z_v - Bw)
            + (b_i / bw) * (Z_v - 1)
            - (Aw / (2 * math.sqrt(2) * Bw))
            * ((2 * avvv / aw) - (b_i / bw))
            * np.log((Z_v + (1 + math.sqrt(2)) * Bw) / (Z_v + (1 - math.sqrt(2)) * Bw))
        )
        # SRK - Peneloux
        fz_i = np.exp(log_coeff)

        return fz_i, Z_v

    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        # Либо реализуем так же, как в базовом коде, либо вынесем в базовый класс EOSFlash
        # Если у всех одинаковая логика расчёта энтальпии, можно оставить общую реализацию в EOSFlash
        pass
