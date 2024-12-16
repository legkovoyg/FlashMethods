import abc

import numpy as np
from interfaces import base_equation
import equation_dto, fluid_dto


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


    def _findroot(self, k: np.ndarray, eps: float = 0.000001, max_iter: int = 1000) -> float:
        """
        Численное нахождение корня методом бисекции.

        Параметры:
        ----------
        k : np.ndarray
            Массив коэффициентов, связанных с расчетом фазового равновесия.
        eps : float, optional
            Точность, с которой нужно найти корень (по умолчанию 1e-6).
        max_iter : int, optional
            Максимальное количество итераций для алгоритма (по умолчанию 1000).

        Возвращаемое значение:
        ----------------------
        float
            Найденное значение корня, которое удовлетворяет заданной точности eps.
        """
        # Вычисляем диапазон возможных значений корня:
        # fv_min и fv_max определяют минимальные и максимальные значения функции,
        # основанные на входном массиве k.
        fv_min: float = 1 / (1 - np.max(k))
        fv_max: float = 1 / (1 - np.min(k))

        # Устанавливаем начальные границы интервала для метода бисекции
        a: float = fv_min + 0.00001
        b: float = fv_max - 0.00001

        # Вектор z_k, представляющий модифицированные значения состава (self._z)
        # и множителя (k - 1).
        z_k: np.ndarray = self._z * (k - 1)

        # Вычисляем значение функции в точке a
        fa: float = float(np.sum(z_k / (1 + a * (k - 1))))


        # Счетчик итераций
        iter_count: int = 0

        # Метод бисекции:
        # Пока разница между границами интервала (a и b) больше заданной точности
        # и число итераций меньше максимального, продолжаем уточнение корня.
        while abs(a - b) > eps and iter_count < max_iter:
            # Средняя точка интервала
            x: float = (a + b) / 2

            # Значение функции в точке x
            fx: float = float(np.sum(z_k / (1 + x * (k - 1))))


            # Если функция меняет знак между a и x, корень лежит в этом интервале.
            # Перемещаем верхнюю границу b в точку x.
            if fa * fx < 0:
                b = x
            else:
                # Иначе корень находится между x и b. Обновляем a и значение fa.
                a = x
                fa = fx

            # Увеличиваем счетчик итераций
            iter_count += 1

        # Если достигнуто максимальное количество итераций, предупреждаем пользователя.
        if iter_count >= max_iter:
            print(f"Warning: _findroot did not converge after {max_iter} iterations")
            # Возвращаем среднее значение границ как приближённое решение.
            return (a + b) / 2

        # Возвращаем среднюю точку интервала, как найденный корень
        return (a + b) / 2


    def _solve_cubic(self, a, b, c, d) -> float:
        """
        Решает кубическое уравнение вида ax^3 + bx^2 + cx + d = 0.

        Параметры:
            a (float): Коэффициент при x^3 (не должен быть равен 0).
            b (float): Коэффициент при x^2.
            c (float): Коэффициент при x.
            d (float): Свободный член.

        Возвращает:
            np.array: Массив из одного или нескольких действительных корней уравнения.
        
        Описание:
            - Используется метод Кардано для решения кубического уравнения.
            - В зависимости от дискриминанта (D), функция возвращает один корень или три действительных корня.
            - Если Q == 0 и R == 0, возвращает троекратно совпадающий корень.
        """
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
