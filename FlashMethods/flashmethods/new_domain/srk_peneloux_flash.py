import numpy as np
from flashmethods.new_domain.base_eos_flash import EOSFlash, log_function

import logging


class SRKPenelouxFlash(EOSFlash):
    """
    Класс для расчёта с использованием уравнения состояния SRK-Peneloux.
    Реализует специфику SRK-Peneloux, включая дополнительный параметр c_i.
    """

    logger = logging.getLogger(__name__)

    @log_function
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

    @log_function
    def _calculate_fugacity(self, P, T, molar_frac, b_i, c_i, c_a_i, isMax=True):
        """
        Рассчитывает фугитивности компонентов и коэффициент сжимаемости Z для уравнения состояния SRK-Peneloux.
        """
        R = self._R

        # Проверка корректности входных данных
        if np.any(molar_frac < 0) or not np.isclose(np.sum(molar_frac), 1.0, rtol=1e-3):
            self.logger.error(
                f"Некорректные мольные доли: {molar_frac}, сумма = {sum(molar_frac)}"
            )
            return None, None

        # Вычисление Aw, Bw и Cw с защитой от некорректных значений
        try:
            aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
            bw = np.dot(molar_frac, b_i)
            cw = np.dot(self._cpen, molar_frac)

            if any(np.isnan([aw, bw, cw])):
                logger.error("Получены nan значения при расчете параметров")
                return None, None

            Aw = aw * P / (R**2 * T**2)
            Bw = bw * P / (R * T)
            Cw = cw * P / (R * T)

            # Коэффициенты кубического уравнения
            a = 1
            b = 3 * Cw - 1
            c = 3 * Cw**2 - Bw**2 - 2 * Cw - Bw + Aw
            d = Cw**3 - Bw**2 * Cw - Cw**2 - Bw * Cw + Aw * Cw - Aw * Bw

            # Решение кубического уравнения
            roots = self._solve_cubic(a, b, c, d)
            if len(roots) == 0:
                logger.error("Не найдены корни кубического уравнения")
                return None, None

            Z_v = np.max(roots) if isMax else np.min(roots)

            # Расчет фугактивности
            avvv = np.dot(molar_frac, c_a_i)

            # Проверяем знак аргумента логарифма
            log_term1 = molar_frac * P
            log_term2 = Z_v + Cw - Bw
            log_term3 = (Z_v + Bw + Cw) / (Z_v + Cw)

            if np.any(log_term1 <= 0) or log_term2 <= 0 or log_term3 <= 0:
                logger.error("Обнаружены отрицательные значения в аргументах логарифма")
                return None, None

            log_coeff = (
                np.log(log_term1)
                - np.log(log_term2)
                + (b_i * P / (R * T) - c_i * P / (R * T)) / log_term2
                - (Aw / Bw) * ((2 * avvv / aw) - (b_i / bw)) * np.log(log_term3)
                - (Aw / Bw) * (b_i + c_i) / (Z_v + Bw + Cw)
                + (Aw / Bw) * c_i / (Z_v + Cw)
            )

            fz_i = np.exp(log_coeff)

            # Проверка результатов
            if np.any(np.isnan(fz_i)) or np.any(np.isinf(fz_i)):
                logger.error("Получены некорректные значения фугактивности")
                return None, None

            return fz_i, Z_v

        except Exception as e:
            logger.error(f"Ошибка при расчете фугактивности: {str(e)}")
            return None, None

    @log_function
    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчёт энтальпии и теплоёмкостей на основе общей логики.
        return: enthalpy_w, Cp_w, Cv_w
        """
        return 1, 1, 1
