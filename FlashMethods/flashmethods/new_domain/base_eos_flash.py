import abc
import logging
import functools
import numpy as np
from flashmethods.dtos import equation_dto
from flashmethods.new_domain.base_flash import Flash


logger = logging.getLogger(__name__)
# Настройка логгера
logging.basicConfig(level=logging.DEBUG, format="%(message)s")


import functools
import inspect
import os


def truncate_array(arr):
    if isinstance(arr, np.ndarray):
        if arr.size > 10:
            return f"{arr[:5].tolist()}, ..., {arr[-5:].tolist()}"
    return arr


def log_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Получаем информацию о файле
        file_path = inspect.getfile(func)
        file_name = os.path.basename(file_path)

        # Начало
        logger.info(f"\n{'='*50}")
        logger.info(f"Файл: {file_name}")
        logger.info(f"Функция: {func.__name__}")
        logger.info(f"Начало выполнения")

        # Входные данные
        if len(args) > 1:  # Пропускаем self
            truncated_args = tuple(truncate_array(arg) for arg in args[1:])
            # logger.info(f"Входные позиционные аргументы: {truncated_args}")

        if kwargs:
            truncated_kwargs = {k: truncate_array(v) for k, v in kwargs.items()}
            # logger.info(f"Входные именованные аргументы: {truncated_kwargs}")

        try:
            # Выполнение функции
            result = func(*args, **kwargs)

            # Выходные данные
            truncated_result = truncate_array(result)
            # logger.info(f"Выходные данные: {truncated_result}")
        except Exception as e:
            logger.exception(f"Неожиданная ошибка: {str(e)}")
            raise
        finally:
            logger.info(f"{'='*50}")

        return result

    return wrapper


class EOSFlash(Flash):
    """
    Базовый класс для всех уравнений состояния типа ПР, СРК, ПР с Пенеле, СРК с Пенеле.
    Содержит общую логику метода calculate, разбитую на этапы.
    Наследники переопределяют только методы _init_eos_parameters и _calculate_fugacity.
    """

    @abc.abstractmethod
    def _init_eos_parameters(self):
        """
        Инициализирует параметры уравнения состояния.

        Возвращает:
            tuple:
                - a_i (ndarray): Параметр "a" для каждого компонента
                - b_i (ndarray): Параметр "b" для каждого компонента
                - c_i (ndarray): Параметр Peneloux (для Peneloux-версий) или None (для обычных)
                - BIPs (ndarray): Бинарные параметры взаимодействия
                - psi_i (ndarray): Коэффициенты для расчёта параметра альфа
                - ac_i (ndarray): Значения параметра "a" при критических условиях
                - params (tuple): Дополнительные параметры для расчёта фугитивности
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

    @log_function
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
        """
        logger.debug("Начало проверки стабильности.")

        def check_phase_stability(K_i, fz_i, phase):
            m = 0
            TS_flag = False
            while m < self._MaxIterationStable:
                Yi = self._z * K_i if phase == "gas" else self._z / K_i
                S = np.sum(Yi)
                y_i = Yi / S

                f_i, Z = self._calculate_fugacity(
                    self._p, self._t, y_i, *c_params, isMax=(phase == "gas")
                )
                Ri = fz_i / (S * f_i) if phase == "gas" else S * f_i / fz_i
                Ri_sum = np.sum((Ri - 1) ** 2)

                logger.debug(f"Итерация стаб. {phase} {m}: Ri_sum={Ri_sum}")
                if Ri_sum < 1e-12:
                    break

                K_i *= Ri
                TS = np.sum(np.log(K_i) ** 2)
                if TS < 1e-4:
                    TS_flag = True
                    break
                m += 1

            return K_i, S, TS_flag

        # Проверка стабильности газовой фазы
        K_iv, Sv, TS_v_flag = check_phase_stability(K_i, fz_i, phase="gas")

        # Проверка стабильности жидкой фазы
        K_il, Sl, TS_l_flag = check_phase_stability(
            self._initial_guess_K(), fz_i, phase="liquid"
        )

        # Определение общей стабильности
        Stable = (
            (TS_l_flag and TS_v_flag)
            or (Sv <= 1 and TS_l_flag)
            or (Sl <= 1 and TS_v_flag)
            or (Sv < 1 and Sl <= 1)
        )
        TestPTF = Stable

        logger.debug(
            f"TS_l_flag={TS_l_flag}, TS_v_flag={TS_v_flag}, Sv={Sv}, Sl={Sl}, Stable={Stable}"
        )

        return Stable, TestPTF, K_iv, K_il, Sv, Sl

    @log_function
    def _initial_guess_K(self):
        """
        Возвращает начальные оценки K-факторов.
        Общая логика для всех EOS, можно вынести в отдельный метод.
        """
        return (
            np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / self._t))
            * self._Pkr
            / self._p
        ) ** 1.0

    @log_function
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
        if not Stable:
            logger.debug(
                f"Система нестабильна (Stable = {Stable}) или TestPTF=True (TestPTF = {TestPTF})"
            )

            K_i = self._choose_best_K(K_iv, K_il)
            logger.debug(f"K_i = {K_i}")
            m = 0
            eps_f = 1
            # Итерации Flash
            iteration_count = 0
            while eps_f > 1e-5 and m < self._MaxIterationFlash:
                logger.debug(f"iteration_count = {iteration_count}")
                W = self._findroot(K_i)
                # Расчет молярных долей жидкой фазы
                self._x_i = self._z / (1 + W * (K_i - 1))
                logger.debug(
                    f"Рассчитанные молярные доли жидкой фазы (до нормализации): {self._x_i}"
                )

                # Нормализация молярных долей жидкой фазы к единице
                sum_x = np.sum(self._x_i)
                if sum_x == 0:
                    logger.error(
                        "Сумма молярных долей жидкой фазы равна нулю. Невозможно нормализовать."
                    )
                    raise ZeroDivisionError(
                        "Сумма молярных долей жидкой фазы равна нулю."
                    )
                self._x_i = self._x_i / sum_x
                logger.debug(
                    f"Молярные доли жидкой фазы после нормализации: {self._x_i}"
                )

                # Расчет молярных долей газовой фазы
                self._y_i = K_i * self._x_i
                # Нормализация молярных долей газовой фазы к единице
                sum_y = np.sum(self._y_i)
                if sum_y == 0:
                    logger.error(
                        "Сумма молярных долей газовой фазы равна нулю. Невозможно нормализовать."
                    )
                    raise ZeroDivisionError(
                        "Сумма молярных долей газовой фазы равна нулю."
                    )
                self._y_i = self._y_i / sum_y
                logger.debug(
                    f"Молярные доли газовой фазы после нормализации: {self._y_i}"
                )

                logger.debug(f"W = {W}")
                logger.debug(f"self._x_i = {sum(self._x_i)}")
                logger.debug(f"self._y_i = {sum(self._y_i)}")
                fw_i, Z_v = self._calculate_fugacity(
                    self._p, self._t, self._y_i, *c_params, isMax=True
                )
                fl_i, Z_l = self._calculate_fugacity(
                    self._p, self._t, self._x_i, *c_params, isMax=False
                )

                Rr = fl_i / fw_i
                K_i *= Rr
                eps_f = np.max(np.abs(Rr - 1))
                logger.debug(f"Итерация Flash {m}: eps_f={eps_f}, W={W}")
                m += 1
                iteration_count += 1

            return W, Z_v, Z_l, False
        else:
            # Если стабильна или TestPTF
            # Аналогично логике из оригинала, выбор W, фаз.
            W, Z_v, Z_l = self._choose_phase_if_stable()
            return W, Z_v, Z_l, True

    @log_function
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
        Kst_v = np.sum((K_iv - 1) ** 2)
        Kst_l = np.sum((K_il - 1) ** 2)
        return K_il if Kst_l > Kst_v else K_iv

    @log_function
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

    @log_function
    def _calculate_thermo_properties(
        self, W, Z_v, Z_l, a_i, BIPs, psi_i, ac_i, b_i, cpen
    ):
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
        enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(
            self._p,
            self._t,
            self._y_i,
            a_i,
            BIPs,
            psi_i,
            ac_i,
            b_i,
            cpen,
            Z_v,
            constants,
        )
        enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(
            self._p,
            self._t,
            self._x_i,
            a_i,
            BIPs,
            psi_i,
            ac_i,
            b_i,
            cpen,
            Z_l,
            constants,
        )

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

        density_y = (
            np.sum((self._mass * self._y_i) / VolumeMy_y) if VolumeMy_y != 0 else 0
        )
        density_x = (
            np.sum((self._mass * self._x_i) / VolumeMy_x) if VolumeMy_x != 0 else 0
        )
        density = (
            np.sum((self._mass * self._z) / self._Volume) if self._Volume != 0 else 0
        )
        self._density = density

        return (
            enthalpy,
            enthalpy_w,
            enthalpy_l,
            Cp,
            Cp_w,
            Cp_l,
            Cv,
            Cv_w,
            Cv_l,
            volume,
            VolumeMy_y,
            VolumeMy_x,
            density,
            density_y,
            density_x,
        )

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

        logger.info("ЭТАП 0 Начало флэш-расчёта EOSFlash.\n")

        # Этап 1: Инициализировать параметры EOS
        a_i, b_i, c_i, BIPs, psi_i, ac_i, c_params = self._init_eos_parameters()
        logger.debug("ЭТАП 1 Параметры EOS инициализированы.\n")

        # Этап 2: Предварительный расчёт фугактивности для исходной смеси z
        K_i = self._initial_guess_K()
        fz_i, self._Z_init = self._calculate_fugacity(
            self._p, self._t, self._z, *c_params, isMax=True
        )
        logger.debug("ЭТАП 2 Предварительная фугактивность рассчитана.\n")

        # Этап 3: Проверка стабильности
        Stable, TestPTF, K_iv, K_il, Sv, Sl = self._perform_stability_check(
            fz_i, self._Z_init, K_i, b_i, c_params
        )
        logger.debug(f"ЭТАП 3 Стабильность: Stable={Stable}, TestPTF={TestPTF}\n")

        # Этап 4: Выполнение Flash-расчёта при необходимости
        W, Z_v, Z_l, stable_after_flash = self._perform_flash(
            Stable, TestPTF, K_iv, K_il, a_i, b_i, c_params
        )
        logger.debug(f"ЭТАП 4 Flash расчёт завершён. W={W}, Z_v={Z_v}, Z_l={Z_l}\n")

        # Этап 5: Расчет энтальпии, объёма, плотности и других свойств
        logger.info("\n=== ЭТАП 5: Расчет термодинамических свойств ===")
        (
            enthalpy,
            enthalpy_w,
            enthalpy_l,
            Cp,
            Cp_w,
            Cp_l,
            Cv,
            Cv_w,
            Cv_l,
            volume,
            VolumeMy_y,
            VolumeMy_x,
            density,
            density_y,
            density_x,
        ) = self._calculate_thermo_properties(
            W, Z_v, Z_l, a_i, BIPs, psi_i, ac_i, b_i, self._cpen
        )

        # Этап 6: Формирование результата
        logger.info("\n=== ЭТАП 6: Формирование итогового результата ===")
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
