import numpy as np
from flashmethods.new_domain.base_eos_flash import EOSFlash


class SRKPenelouxFlash(EOSFlash):
    """
    Класс для расчёта с использованием уравнения состояния SRK-Peneloux.
    Реализует специфику SRK-Peneloux, включая дополнительный параметр c_i.
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

    def _calculate_fugacity(self, P, T, molar_frac, b_i, c_i, c_a_i, isMax=True):
        """
        Рассчитывает фугитивности компонентов и коэффициент сжимаемости Z для уравнения состояния SRK-Peneloux.
        """
        R = self._R
        # Проверка корректности входных данных
        if np.any(molar_frac < 0) or not np.isclose(np.sum(molar_frac), 1.0, rtol=1e-3):
            print(f"Некорректные мольные доли: {molar_frac}, сумма = {sum(molar_frac)}")
            return None, None

        # Вычисление Aw, Bw и Cw с защитой от некорректных значений
        try:
            aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
            bw = np.dot(molar_frac, b_i)
            cw = np.dot(self._cpen, molar_frac)

            Aw = aw * P / (R**2 * T**2)
            Bw = bw * P / (R * T)
            Cw = cw * P / (R * T)

            Biw = b_i * P / (R * T)
            Ciw = c_i * P / (R * T)

            # Коэффициенты кубического уравнения
            a = 1
            b = 3 * Cw - 1
            c = 3 * Cw**2 - Bw**2 - 2 * Cw - Bw + Aw
            d = Cw**3 - Bw**2 * Cw - Cw**2 - Bw * Cw + Aw * Cw - Aw * Bw

            # Решение кубического уравнения
            roots = self._solve_cubic(a, b, c, d)

            Z_v = np.max(roots) if isMax else np.min(roots)

            # Расчет фугактивности
            avvv = np.dot(molar_frac, c_a_i)

            log_fugacity = (
                np.log(molar_frac * P)
                - np.log(Z_v + Cw - Bw)
                + (Biw - Ciw) / (Z_v + Cw - Bw)
                - (Aw / Bw)
                * ((2 * avvv / aw) - (b_i / bw))
                * np.log((Z_v + Bw + Cw) / (Z_v + Cw))
                - (Aw / Bw) * (Biw + Ciw) / (Z_v + Bw + Cw)
                + (Aw / Bw) * Ciw / (Z_v + Cw)
            )

            fz_i = np.exp(log_fugacity)

            return fz_i, Z_v

        except Exception as e:
            print(f"Ошибка при расчете фугактивности: {str(e)}")
            return None, None

    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчёт энтальпии и теплоёмкостей на основе общей логики.
        return: enthalpy_w, Cp_w, Cv_w
        """
        Tref = 273.15

        # Вычисление разностей температур для степеней от 1 до 4 с помощью цикла
        powers = [1, 2, 3, 4]
        dT = [T ** p - Tref ** p for p in powers]

        Hid = (
                self._Cp1_id * dT[0]
                + 0.5 * self._Cp2_id * dT[1]
                + self._Cp3_id * dT[2] / 3
                + 0.25 * self._Cp4_id * dT[3]
        )

        # Производные для Yi и ai
        dYi_dT = -psi_i * np.power(self._Tkr * T, -0.5) * (1 + psi_i * (1 - np.sqrt(T / self._Tkr)))
        d2Yi_dT2 = 0.5 * psi_i * np.power(self._Tkr, -0.5) * (1 + psi_i) * np.power(T, -1.5)

        dai_dT = ac_i * dYi_dT
        d2ai_dT2 = ac_i * d2Yi_dT2

        R = self._R

        # Вычисление am
        am = (1 - BIPs) * np.sqrt(a_i[:, None] * a_i[None, :]) * molar_frac[None, :] * molar_frac[:, None]
        am = np.sum(am)

        # Производная am по T
        dam_dT = (1 - BIPs) * np.power(a_i[:, None] * a_i[None, :], -0.5) * molar_frac[None, :] * molar_frac[:,
                                                                                                  None] * (
                         dai_dT[None, :] * a_i[:, None] + dai_dT[:, None] * a_i[None, :]
                 )
        dam_dT = 0.5 * np.sum(dam_dT)

        dAm_dT = dam_dT * P / ((R ** 2) * (T ** 2))

        bm = np.dot(molar_frac, b_i)
        cm = np.dot(cpen, molar_frac)

        Am = am * P / (R ** 2 * T ** 2)
        Bm = bm * P / (R * T)
        Cm = cm * P / (R * T)

        # Расчёт остаточной энтальпии
        H_residual = 1000 * R * T * (Z - 1 - ((Am - T * dAm_dT) / Bm) * np.log((Z + Cm + Bm) / (Z + Cm)))

        enthalpy = H_residual + np.dot(molar_frac, Hid)

        return enthalpy

    def _calculate_Cp_Cv(
            self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчёт энтальпии и теплоёмкостей на основе общей логики.
        return: Cp_w, Cv_w
        """
        return 1,1



