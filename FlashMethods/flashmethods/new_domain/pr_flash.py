import numpy as np
from flashmethods.new_domain.base_eos_flash import EOSFlash


# Peng-Robinson EOS
class PRFlash(EOSFlash):
    def _init_eos_parameters(self):
        """Инициализация параметров для уравнения Пенг-Робинсона"""
        ac_i = 0.457235 * self._R**2 * self._Tkr**2.0 / self._Pkr
        psi_i = 0.37464 + 1.54226 * self._w - 0.26992 * self._w**2
        alpha_i = (1 + psi_i * (1 - (self._t / self._Tkr) ** 0.5)) ** 2
        a_i = ac_i * alpha_i
        b_i = 0.077796 * self._R * self._Tkr / self._Pkr
        BIPs = self._c
        c_a_i = (1 - BIPs) * np.sqrt(a_i[:, None] * a_i[None, :])

        c_i = None
        return a_i, b_i, c_i, BIPs, psi_i, ac_i, (b_i, c_i, c_a_i)

    def _calculate_fugacity(self, P, T, molar_frac, *args, isMax=True):
        b_i, _, c_a_i = args  # распаковываем только нужные параметры
        R = self._R

        # Вычисление параметров смеси
        aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
        bw = np.dot(molar_frac, b_i)

        # Приведенные параметры
        Aw = aw * P / (R**2 * T**2)
        Bw = bw * P / (R * T)

        # Коэффициенты кубического уравнения PR
        a = 1
        b = -(1 - Bw)
        c = Aw - 3 * Bw**2 - 2 * Bw
        d = -(Aw * Bw - Bw**2 - Bw**3)

        # Z-фактор
        roots = self._solve_cubic(a, b, c, d)
        roots = roots[np.logical_and(np.imag(roots) == 0, roots > 0)]
        Z_v = np.max(roots) if isMax else np.min(roots)

        # Параметры для фугитивности
        avvv = np.dot(molar_frac, c_a_i)
        sqrt_2 = np.sqrt(2)

        # Расчет логарифма коэффициента фугитивности PR
        log_coeff = (
            np.log(molar_frac * P)
            - np.log(Z_v - Bw)
            + (b_i / bw) * (Z_v - 1)
            - (Aw / (2 * sqrt_2 * Bw))
            * ((2 * avvv / aw) - (b_i / bw))
            * np.log((Z_v + (1 + sqrt_2) * Bw) / (Z_v + (1 - sqrt_2) * Bw))
        )

        return np.exp(log_coeff), Z_v

    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчёт энтальпии и теплоёмкостей на основе общей логики.
        return: enthalpy_w, Cp_w, Cv_w
        """
        return 1, 1, 1
