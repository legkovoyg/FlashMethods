import abc
import logging
import numpy as np
from dtos import equation_dto
from base_flash import Flash

logger = logging.getLogger(__name__)

class EOSFlash(Flash):
    """
    Базовый класс для всех уравнений состояния типа ПР, СРК, ПР с Пенеле, СРК с Пенеле.
    Содержит общую логику метода calculate, разбитую на этапы.
    Наследники переопределяют только методы _init_eos_parameters и _calculate_fugacity.
    """

    @abc.abstractmethod
    def _init_eos_parameters(self):
        """
        Инициализирует параметры уравнения состояния, такие как ac_i, psi_i, alpha_i, a_i, b_i, c_i (если нужно).
        Должен быть переопределён в потомках, так как PR, SRK и их Peneloux-версии рассчитывают эти параметры по-разному.
        """
        pass

    @abc.abstractmethod
    def _calculate_fugacity(self, P, T, molar_frac, *args, isMax=True):
        """
        Рассчитать фугактивности и Z-фактор для текущего EOS.
        Параметры (например, b_i, c_a_i) будут зависеть от того, что вернёт _init_eos_parameters.
        Должен быть переопределён в потомках.
        """
        pass

    def _perform_stability_check(self, fz_i, Z_init, K_i, b_i, c_params):
        """
        Выполняет проверку термодинамической стабильности газовой и жидкой фаз 
        для уравнений состояния (EOS).

        Параметры:
            fz_i (ndarray): Фугитивности компонентов в исходной фазе.
            Z_init (float): Начальное значение коэффициента сжимаемости (Z-фактора).
            K_i (ndarray): Начальные коэффициенты распределения компонентов между фазами.
            b_i (ndarray): Параметры "b" для каждого компонента.
            c_params (tuple): Дополнительные параметры уравнения состояния, 
                            такие как матрица взаимодействий (c_a_i) и параметры "b_i" и "c_i".

        Возвращает:
            tuple:
                - Stable (bool): Флаг стабильности системы (True - стабильна, False - нестабильна).
                - TestPTF (bool): Флаг результата проверки стабильности (True - устойчива, False - неустойчива).
                - K_iv (ndarray): Коэффициенты распределения для газовой фазы.
                - K_il (ndarray): Коэффициенты распределения для жидкой фазы.
                - Sv (float): Коэффициент нормы для газовой фазы.
                - Sl (float): Коэффициент нормы для жидкой фазы.

        Описание:
            - Алгоритм проверяет термодинамическую стабильность газовой и жидкой фаз.
            - Для газовой фазы проводится итерационная проверка на основе коэффициентов распределения K_i.
            - Аналогично для жидкой фазы определяется стабильность.
            - Условия стабилизации фаз определяются по суммам квадратичных отклонений (Ri) и 
            функции стабильности (TS).

        Примечания:
            - Максимальное число итераций задаётся параметром `_MaxIterationStable`.
            - Логгируются промежуточные результаты каждой итерации для отслеживания процесса.
            - Итерации прекращаются, если суммарное отклонение Ri падает ниже заданного порога 
            или достигается условие стабилизации (TS < 1e-4).
        """
        logger.debug("Начало проверки стабильности газовой фазы.")
        m = 0
        TS_v_flag = 0
        # Проверка стабильности для газовой фазы
        while m < self._MaxIterationStable:
            Yi_v = self._z * K_i
            Sv = np.sum(Yi_v)
            self._y_i = Yi_v / Sv

            fw_i, Z_v = self._calculate_fugacity(self._p, self._t, self._y_i, *c_params, isMax=True)
            Ri = fz_i / (Sv * fw_i)
            Ri_v = np.sum((Ri - 1)**2)

            logger.debug(f"Итерация стаб. газа {m}: Ri_v={Ri_v}")
            if Ri_v < 1e-12:
                break

            K_i *= Ri
            TS_v = np.sum(np.log(K_i)**2)
            if TS_v < 1e-4:
                TS_v_flag = 1
                break
            m += 1

        K_iv = K_i.copy()

        logger.debug("Проверка стабильности жидкой фазы.")
        # Проверка стабильности для жидкой фазы
        K_i = self._initial_guess_K()
        fz_i, Z_v = self._calculate_fugacity(self._p, self._t, self._z, *c_params, isMax=True)

        ml = 0
        TS_l_flag = 0
        while ml < self._MaxIterationStable:
            Yi_l = self._z / K_i
            Sl = np.sum(Yi_l)
            self._x_i = Yi_l / Sl

            fl_i, Z_l = self._calculate_fugacity(self._p, self._t, self._x_i, *c_params, isMax=False)
            Ri = Sl * fl_i / fz_i
            Ri_l = np.sum((Ri - 1)**2)
            logger.debug(f"Итерация стаб. жидкости {ml}: Ri_l={Ri_l}")

            if Ri_l < 1e-12:
                break

            K_i *= Ri
            TS = np.sum(np.log(K_i)**2)
            if TS < 1e-4:
                TS_l_flag = 1
                break
            ml += 1

        K_il = K_i.copy()

        # Определение общей стабильности:
        Sv = np.sum(self._y_i)
        Sl = np.sum(self._x_i)
        if ((TS_l_flag == 1 and TS_v_flag == 1)
            or (Sv <= 1 and TS_l_flag == 1)
            or (Sl <= 1 and TS_v_flag == 1)
            or (Sv < 1 and Sl <= 1)):
            Stable = True
            TestPTF = True
        else:
            Stable = False
            TestPTF = False

        return Stable, TestPTF, K_iv, K_il, Sv, Sl

    def _initial_guess_K(self):
        """
        Возвращает начальные оценки K-факторов.
        Общая логика для всех EOS, можно вынести в отдельный метод.
        """
        return (
            np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / self._t))
            * self._Pkr / self._p
        )**1.0

    def _perform_flash(self, Stable, TestPTF, K_iv, K_il, a_i, b_i, c_params):
        """
        Выполняет Flash-расчёт для определения фазового равновесия системы.

        Параметры:
            Stable (bool): Флаг стабильности системы (True - стабильна, False - нестабильна).
            TestPTF (bool): Результат теста на фазовую стабильность (True - устойчивость фаз проверена).
            K_iv (ndarray): Коэффициенты распределения для газовой фазы.
            K_il (ndarray): Коэффициенты распределения для жидкой фазы.
            a_i (ndarray): Параметр "a" для каждого компонента.
            b_i (ndarray): Параметр "b" для каждого компонента.
            c_params (tuple): Дополнительные параметры уравнения состояния (например, c_a_i).

        Возвращает:
            tuple:
                - W (float): Качество смеси (мольная доля газа в двухфазной системе).
                - Z_v (float): Коэффициент сжимаемости для газовой фазы.
                - Z_l (float): Коэффициент сжимаемости для жидкой фазы.
                - isStable (bool): Флаг стабильности системы после расчёта (True - стабильна, False - нестабильна).

        Описание:
            - Если система нестабильна или требует проверки (TestPTF=True), выполняется Flash-расчёт.
            - Итерационный процесс корректирует коэффициенты распределения (K_i) для достижения равновесия фаз.
            - Используются метод нахождения корня (`_findroot`) и расчёт фугитивностей для фазовых состояний.
            - Если система стабильна, выбирается фаза на основе текущего состояния (одна или двухфазная система).

        Алгоритм:
            - Итерации продолжаются до тех пор, пока максимальное отклонение (eps_f) не станет меньше 1e-5 
            или не достигнуто максимальное число итераций (`_MaxIterationFlash`).
            - В случае стабильной системы производится выбор фазового состояния с использованием метода `_choose_phase_if_stable`.

        Логирование:
            - В процессе итераций Flash логируются значения отклонений (eps_f), коэффициентов распределения и доли газа (W).
        """
        if not Stable or TestPTF:
            logger.debug("Система нестабильна или TestPTF=True, выполняем Flash-расчет.")

            K_i = self._choose_best_K(K_iv, K_il)
            m = 0
            eps_f = 1
            # Итерации Flash
            while eps_f > 1e-5 and m < self._MaxIterationFlash:
                W = self._findroot(K_i)
                self._x_i = self._z / (1 + W * (K_i - 1))
                self._y_i = K_i * self._x_i

                fw_i, Z_v = self._calculate_fugacity(self._p, self._t, self._y_i, *c_params, isMax=True)
                fl_i, Z_l = self._calculate_fugacity(self._p, self._t, self._x_i, *c_params, isMax=False)

                Rr = fl_i / fw_i
                K_i *= Rr
                eps_f = np.max(np.abs(Rr - 1))
                logger.debug(f"Итерация Flash {m}: eps_f={eps_f}, W={W}")
                m += 1

            return W, Z_v, Z_l, False
        else:
            # Если стабильна или TestPTF
            # Аналогично логике из оригинала, выбор W, фаз.
            W, Z_v, Z_l = self._choose_phase_if_stable()
            return W, Z_v, Z_l, True

    def _choose_best_K(self, K_iv, K_il):
        """
        Выбирает лучшие коэффициенты распределения (K-факторы) для старта Flash-расчёта.

        Параметры:
            K_iv (ndarray): К-факторы, рассчитанные для газовой фазы.
            K_il (ndarray): К-факторы, рассчитанные для жидкой фазы.

        Возвращает:
            ndarray: Оптимальные К-факторы для начала итераций Flash.

        Описание:
            - Расчёт основан на оценке квадратичного отклонения К-факторов от единицы.
            - Если отклонение для жидкой фазы (Kst_l) больше, выбираются K_il, иначе K_iv.
            - Используется для оптимального выбора начальных условий в зависимости от фазовой стабильности.
        """
        Kst_v = np.sum((K_iv - 1)**2)
        Kst_l = np.sum((K_il - 1)**2)
        return K_il if Kst_l > Kst_v else K_iv

    def _choose_phase_if_stable(self):
        """
        Определяет фазу системы, если она признана стабильной, 
        на основе объёма и температуры.

        Возвращает:
            tuple:
                - W (float): Качество смеси (мольная доля газа в фазовой системе).
                - Z_v (float): Коэффициент сжимаемости для газовой фазы.
                - Z_l (float): Коэффициент сжимаемости для жидкой фазы.

        Описание:
            - Решение принимается на основе объёма системы и контрольного параметра температуры.
            - Если температура ниже заданного порога (`_T_pkr < 260`), выбирается газовая фаза.
            - При объёме, превышающем контрольное значение, также выбирается газовая фаза.
            - Если оба условия не выполняются, выбирается жидкая фаза.
            - Обновляются мольные доли компонентов в газовой (`_y_i`) и жидкой (`_x_i`) фазах.
        """
        Volume = 1000.0 * (self._Z_init * self._R * self._t / self._p)
        if self._T_pkr < 260:
            W = 1
            self._x_i = np.zeros(self._N)
            self._y_i = self._z
            Z_v = self._Z_init
            Z_l = 0.0
        else:
            if Volume * self._t**2 > self._Control_phase:
                W = 1
                self._x_i = np.zeros(self._N)
                self._y_i = self._z
                Z_v = self._Z_init
                Z_l = 0.0
            else:
                W = 0
                self._x_i = self._z
                self._y_i = np.zeros(self._N)
                Z_v = 0
                Z_l = self._Z_init
        return W, Z_v, Z_l

    def _calculate_thermo_properties(self, W, Z_v, Z_l, a_i, BIPs, psi_i, ac_i, b_i, cpen):
        """
        Расчёт термодинамических свойств: энтальпии, теплоёмкостей (Cp, Cv), объёма и плотности.

        Параметры:
            W (float): Качество смеси (мольная доля газа в двухфазной системе).
            Z_v (float): Коэффициент сжимаемости для газовой фазы.
            Z_l (float): Коэффициент сжимаемости для жидкой фазы.
            a_i, BIPs, psi_i, ac_i, b_i, cpen: Параметры уравнения состояния.

        Возвращает:
            tuple: Энтальпия, теплоёмкости (Cp, Cv), объём и плотности (общие и фазовые).

        Примечание:
            - Метод в текущем виде может работать некорректно, требуется дополнительная проверка и доработка.
        """
        constants = "1"  # Или другой режим по необходимости
        enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(self._p, self._t, self._y_i, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z_v, constants)
        enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(self._p, self._t, self._x_i, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z_l, constants)

        if W == 1:
            enthalpy = enthalpy_w
            Cp = Cp_w
            Cv = Cv_w
        elif W == 0:
            enthalpy = enthalpy_l
            Cp = Cp_l
            Cv = Cv_l
        else:
            enthalpy = enthalpy_w * W + enthalpy_l * (1 - W)
            Cp = Cp_w * W + Cp_l * (1 - W)
            Cv = Cv_w * W + Cv_l * (1 - W)

        cpen_mix_y = np.sum(self._cpen * self._y_i)
        cpen_mix_x = np.sum(self._cpen * self._x_i)
        VolumeMy_y = Z_v * self._R * 1000 * self._t / self._p - cpen_mix_y
        VolumeMy_x = Z_l * self._R * 1000 * self._t / self._p - cpen_mix_x
        volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
        self._Volume = volume

        density_y = np.sum((self._mass * self._y_i) / VolumeMy_y) if VolumeMy_y != 0 else 0
        density_x = np.sum((self._mass * self._x_i) / VolumeMy_x) if VolumeMy_x != 0 else 0
        density = np.sum((self._mass * self._z) / self._Volume) if self._Volume != 0 else 0
        self._density = density

        return enthalpy, enthalpy_w, enthalpy_l, Cp, Cp_w, Cp_l, Cv, Cv_w, Cv_l, volume, VolumeMy_y, VolumeMy_x, density, density_y, density_x

    @abc.abstractmethod
    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        """
        Расчет энтальпии, Cp, Cv (есть во всех EOS),
        но если вы хотите, можете тоже сделать в базовом и переопределять только при необходимости.
        Предположим, что в данном примере общая логика одинаковая, 
        тогда вынести в базовый класс.
        """
        pass

    def calculate(self) -> equation_dto.EquationResultDTO:
        """
        Общий метод расчёта. Разбит на этапы.
        - Этап 1: Инициализировать параметры EOS
        - Этап 2: Предварительный расчёт фугактивности исходной смеси
        - Этап 3: Проверка стабильности
        - Этап 4: Flash-расчёт при необходимости
        - Этап 5: Расчет энтальпии, объёма, плотности
        - Этап 6: Формирование результата
        """

        logger.info("Начало флэш-расчёта EOSFlash.")

        # Этап 1: Инициализировать параметры EOS
        a_i, b_i, BIPs, psi_i, ac_i, c_params = self._init_eos_parameters()
        logger.debug("Параметры EOS инициализированы.")

        # Этап 2: Предварительный расчёт фугактивности для исходной смеси z
        K_i = self._initial_guess_K()
        fz_i, self._Z_init = self._calculate_fugacity(self._p, self._t, self._z, *c_params, isMax=True)
        logger.debug("Предварительная фугактивность рассчитана.")

        # Этап 3: Проверка стабильности
        Stable, TestPTF, K_iv, K_il, Sv, Sl = self._perform_stability_check(fz_i, self._Z_init, K_i, b_i, c_params)
        logger.debug(f"Стабильность: Stable={Stable}, TestPTF={TestPTF}")

        # Этап 4: Выполнение Flash-расчёта при необходимости
        W, Z_v, Z_l, stable_after_flash = self._perform_flash(Stable, TestPTF, K_iv, K_il, a_i, b_i, c_params)
        logger.debug(f"Flash расчёт завершён. W={W}, Z_v={Z_v}, Z_l={Z_l}")

        # Этап 5: Расчет энтальпии, объёма, плотности и других свойств
        (enthalpy, enthalpy_w, enthalpy_l, Cp, Cp_w, Cp_l, Cv, Cv_w, Cv_l,
         volume, VolumeMy_y, VolumeMy_x, density, density_y, density_x) = self._calculate_thermo_properties(W, Z_v, Z_l, a_i, BIPs, psi_i, ac_i, b_i, self._cpen)

        # Этап 6: Формирование результата
        result = equation_dto.EquationResultDTO(
            w=W,
            z_v=Z_v,
            z_l=Z_l,
            x_i=self._x_i.tolist(),
            y_i=self._y_i.tolist(),
            is_stable=Stable,
            m=self._MaxIterationFlash,  # Например, или реально считать итерации
            enthalpy=enthalpy,
            enthalpy_w=enthalpy_w,
            enthalpy_l=enthalpy_l,
            cp=Cp,
            cp_w=Cp_w,
            cp_l=Cp_l,
            cv=Cv,
            cv_w=Cv_w,
            cv_l=Cv_l,
            volume=volume,
            volume_my_y=VolumeMy_y,
            volume_my_x=VolumeMy_x,
            density=density,
            density_y=density_y,
            density_x=density_x,
        )

        logger.info("Флэш-расчёт завершён успешно.")
        return result
