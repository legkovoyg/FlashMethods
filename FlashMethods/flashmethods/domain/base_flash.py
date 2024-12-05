import abc

import numpy as np
from interfaces import base_equation
from dtos import equation_dto, fluid_dto


class Flash(base_equation.Equation):
    """
    Базовый класс для уравнения состояния
    """

    def __init__(self, fluid: fluid_dto.FluidDTO) -> None:
        """
        Инициализировать переменные
        :param fluid: модель со свойствами флюида
        """

        self._t = fluid.t  # температура
        self._p = fluid.p  # давление
        self._N = len(fluid.components)  # количество компонент N
        self._R = 0.00831472998267014  # универсальная газовая МПа*м3/(кмоль*К)
        self._MaxIterationStable = 30  # Максимальное количество итераций для проверки 1
        self._MaxIterationFlash = 400  # Максимальное количество итераций для Flash-расчета
        self._component_name = np.empty(self._N, dtype="S10")
        self._z = np.zeros(self._N, dtype=np.float64)
        self._mass = np.zeros(self._N, dtype=np.float64)
        self._Pkr = np.zeros(self._N, dtype=np.float64)
        self._Tkr = np.zeros(self._N, dtype=np.float64)
        self._Vkr = np.zeros(self._N, dtype=np.float64)
        self._w = np.zeros(self._N, dtype=np.float64)
        self._cpen = np.zeros(self._N, dtype=np.float64)
        self._x_i = np.zeros(self._N, dtype=np.float64)
        self._y_i = np.zeros(self._N, dtype=np.float64)
        # Новое
        self._T_boil = np.zeros(self._N, dtype=np.float64)
        self._density_liq_phase = np.zeros(self._N, dtype=np.float64)
        self._Cp1_id = np.zeros(self._N, dtype=np.float64)
        self._Cp2_id = np.zeros(self._N, dtype=np.float64)
        self._Cp3_id = np.zeros(self._N, dtype=np.float64)
        self._Cp4_id = np.zeros(self._N, dtype=np.float64)
        self._Kw = np.zeros(self._N, dtype=np.float64)
        self._Cf = np.zeros(self._N, dtype=np.float64)
        self._Volume = np.zeros(1, dtype=np.float64)
        self._density = np.zeros(1, dtype=np.float64)

        if fluid.BIPs is None:
            self._c = np.zeros((self._N, self._N), dtype=np.float64)
        else:
            self._c = np.array(fluid.BIPs, dtype=np.float64)

        for i, component in enumerate(fluid.components):
            self._component_name[i] = component.component_name
            self._z[i] = component.z
            self._mass[i] = component.mass
            self._Pkr[i] = component.Pkr
            self._Tkr[i] = component.Tkr
            self._Vkr[i] = component.Vkr
            self._w[i] = component.w
            self._cpen[i] = component.cpen / 1000.0
            # Новое
            self._T_boil[i] = component.T_boil
            self._density_liq_phase[i] = component.density_liq_phase

        self._z /= sum(self._z)  # Нормализация состава
        self._V_pkr = np.dot(self._z, self._Vkr)
        self._T_pkr = np.dot(self._z, self._Tkr)
        self._Control_phase = self._V_pkr * self._T_pkr * self._T_pkr

    def _findroot(self, k, eps=0.000001, max_iter=1000):
        fv_min = 1 / (1 - np.max(k))
        fv_max = 1 / (1 - np.min(k))
        a = fv_min + 0.00001
        b = fv_max - 0.00001
        z_k = self._z * (k - 1)
        fa = np.sum(z_k / (1 + a * (k - 1)))
        iter_count = 0

        while abs(a - b) > eps and iter_count < max_iter:
            x = (a + b) / 2
            fx = np.sum(z_k / (1 + x * (k - 1)))
            if fa * fx < 0:
                b = x
            else:
                a = x
                fa = fx
            iter_count += 1

        if iter_count >= max_iter:
            print(f"Warning: _findroot did not converge after {max_iter} iterations")
            return (a + b) / 2

        return (a + b) / 2

    def _solve_cubic(self, a, b, c, d):
        inv_a = 1.0 / a
        b_a = inv_a * b
        b_a2 = b_a * b_a
        c_a = inv_a * c
        d_a = inv_a * d

        Q = (3 * c_a - b_a2) / 9
        R = (9 * b_a * c_a - 27 * d_a - 2 * b_a * b_a2) / 54
        Q3 = Q * Q * Q
        D = Q3 + R * R
        b_a_3 = (1.0 / 3.0) * b_a

        if Q == 0:
            if R == 0:
                x0 = x1 = x2 = -b_a_3
                return np.array([x0, x1, x2])
            else:
                cube_root = (2 * R) ** (1.0 / 3.0)
                x0 = cube_root - b_a_3
                return np.array([x0])

        if D <= 0:
            theta = np.arccos(R / np.sqrt(-Q3))
            sqrt_Q = np.sqrt(-Q)
            x0 = 2 * sqrt_Q * np.cos(theta / 3.0) - b_a_3
            x1 = 2 * sqrt_Q * np.cos((theta + 2 * np.pi) / 3.0) - b_a_3
            x2 = 2 * sqrt_Q * np.cos((theta + 4 * np.pi) / 3.0) - b_a_3
            return np.array([x0, x1, x2])

        AD = 0.0
        BD = 0.0
        R_abs = np.fabs(R)
        if R_abs > 2.2204460492503131e-16:
            AD = (R_abs + np.sqrt(D)) ** (1.0 / 3.0)
            AD = AD if R >= 0 else -AD
            BD = -Q / AD

        x0 = AD + BD - b_a_3
        return np.array([x0])

    @abc.abstractmethod
    def calculate(self) -> equation_dto.EquationResultDTO:
        """
        Выполнить расчет
        :return: рассчитанные параметры
        """

        raise NotImplementedError
