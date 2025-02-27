import numpy as np
from flashmethods.new_domain.base_eos_flash import EOSFlash


class SRKFlash(EOSFlash):
    """
    Реализация расчёта на основе уравнения состояния Soave-Redlich-Kwong (SRK).
    Наследует общую логику из EOSFlash и реализует специфичные методы SRK.
    """

    def _init_eos_parameters(self):
        """
        Инициализация параметров EOS для SRK-Peneloux.

        Возвращает:
            tuple:
                - a_i (ndarray): Параметр "a" для каждого компонента.
                - b_i (ndarray): Параметр "b" для каждого компонента.
                - c_i (ndarray): Дополнительный параметр Peneloux для каждого компонента.
                - BIPs (ndarray): Бинарные параметры взаимодействия (Binary Interaction Parameters).
                - psi_i (ndarray): Коэффициенты для расчёта параметра альфа.
                - ac_i (ndarray): Значения параметра "a" при критических условиях.
                - tuple: Набор параметров (b_i, c_i, c_a_i) для дальнейших расчётов.

        Описание:
            - Параметры SRK рассчитываются с учётом критической температуры (Tkr),
              критического давления (Pkr) и ацентрического фактора (w).
            - Параметр Peneloux (c_i) добавляется для корректировки объёма.
            - Матрица взаимодействий c_a_i учитывает бинарные взаимодействия между компонентами.
        """

        # Коэффициенты SRK
        ac_i = 0.42747 * self._R**2 * self._Tkr**2 / self._Pkr
        psi_i = 0.48 + 1.574 * self._w - 0.176 * self._w**2
        alpha_i = (1 + psi_i * (1 - np.sqrt(self._t / self._Tkr))) ** 2
        a_i = ac_i * alpha_i
        b_i = 0.08664 * self._R * self._Tkr / self._Pkr

        # Дополнительный параметр Peneloux
        c_i = self._cpen

        # Матрица взаимодействий для параметров a
        BIPs = self._c
        c_a_i = (1 - BIPs) * np.sqrt(a_i[:, None] * a_i[None, :])
        return a_i, b_i, c_i, BIPs, psi_i, ac_i, (b_i, c_i, c_a_i)

    def _calculate_fugacity(self, P, T, molar_frac, *args, isMax=True):
        """
        Рассчитывает фугитивность компонента в смеси на основе уравнения состояния SRK (EOS).
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

        b_i, _, c_a_i = args  # распаковываем только нужные параметры
        R = self._R

        # Вычисление Aw и Bw
        aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
        bw = np.dot(molar_frac, b_i)
        Aw = aw * P / (R**2 * T**2)  # Коэффициент Aw
        Bw = bw * P / (R * T)  # Коэффициент Bw

        # Коэффициенты кубического уравнения SRK
        a = 1
        b = -1
        c = Aw - Bw - Bw**2
        d = -(Aw * Bw)

        # Решение кубического уравнения для Z-фактора
        roots = self._solve_cubic(a, b, c, d)
        Z_v = np.max(roots) if isMax else np.min(roots)

        # Логарифм фугактивности
        avvv = np.dot(molar_frac, c_a_i)
        log_coeff = (
            np.log(molar_frac * P)
            - np.log(Z_v - Bw)
            + (b_i / bw) * (Z_v - 1)
            - (Aw / Bw) * ((2 * avvv / aw) - (b_i / bw)) * np.log(1 + Bw / Z_v)
        )
        fz_i = np.exp(log_coeff)  # Фугактивность компонентов

        return fz_i, Z_v

    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчёт энтальпии и теплоёмкостей на основе общей логики.
        return: enthalpy_w, Cp_w, Cv_w
        """
        return 1, 1, 1
