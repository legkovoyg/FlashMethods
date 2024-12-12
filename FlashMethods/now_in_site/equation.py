import numpy as np
from time import time
from abc import ABC, abstractmethod
import math

class Flash(ABC):

    def __init__(self, Fluid):
        self._N = len(Fluid.components)  # количество компонент N
        self._R = 0.00831472998267014  # универсальная газовая МПа*м3/(кмоль*К)
        self._MaxIterationStable = (
            30  # Максимальное количество итераций для проверки 1
        )
        self._MaxIterationFlash = (
            400  # Максимальное количество итераций для Flash-расчета
        )
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
        if Fluid.BIPs is None:
            self._c = np.zeros((self._N, self._N), dtype=np.float64)
        else:
            self._c = Fluid.BIPs.to_numpy(dtype=np.float64)
        z_sum = 0.0
        for i, component in enumerate(Fluid.components):
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
            print(f'Warning: _findroot did not converge after {max_iter} iterations')
            return (a + b) / 2

        return (a + b) / 2

    def _solve_cubic(self, a, b, c, d):
        inv_a = 1. / a
        b_a = inv_a * b
        b_a2 = b_a * b_a
        c_a = inv_a * c
        d_a = inv_a * d

        Q = (3 * c_a - b_a2) / 9
        R = (9 * b_a * c_a - 27 * d_a - 2 * b_a * b_a2) / 54
        Q3 = Q * Q * Q
        D = Q3 + R * R
        b_a_3 = (1. / 3.) * b_a

        if Q == 0:
            if R == 0:
                x0 = x1 = x2 = - b_a_3
                return np.array([x0, x1, x2])
            else:
                cube_root = (2 * R) ** (1. / 3.)
                x0 = cube_root - b_a_3
                return np.array([x0])

        if D <= 0:
            theta = np.arccos(R / np.sqrt(-Q3))
            sqrt_Q = np.sqrt(-Q)
            x0 = 2 * sqrt_Q * np.cos(theta / 3.0) - b_a_3
            x1 = 2 * sqrt_Q * np.cos((theta + 2 * np.pi) / 3.0) - b_a_3
            x2 = 2 * sqrt_Q * np.cos((theta + 4 * np.pi) / 3.0) - b_a_3
            return np.array([x0, x1, x2])

        AD = 0.
        BD = 0.
        R_abs = np.fabs(R)
        if R_abs > 2.2204460492503131e-16:
            AD = (R_abs + np.sqrt(D)) ** (1. / 3.)
            AD = AD if R >= 0 else -AD
            BD = -Q / AD

        x0 = AD + BD - b_a_3
        return np.array([x0])
    
    @abstractmethod
    def vle(self, P, T):
        pass

# class PR_Flash (Flash):
#     def __init__(self, Fluid):
#         super().__init__(Fluid)  # Вызов конструктора родительского класса
#
#     def _calculate(self, P, T, molar_frac, b_i, c_a_i, isMax):
#         R = self._R
#
#         aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
#         bw = np.dot(molar_frac, b_i)
#         Aw = aw * P / (R**2 * T**2)
#         Bw = bw * P / (R * T)
#
#         # SRK - Peneloux cubic equation coefficients
#         a = 1
#         b = -(1 - Bw)
#         c = Aw - 3 * Bw**2 - 2 * Bw
#         d = -(Aw * Bw - Bw**2 - Bw**3)
#
#         roots = self._solve_cubic(a, b, c, d)
#         Z_v = np.max(roots) if isMax else np.min(roots)
#
#         avvv = np.dot(molar_frac, c_a_i)
#
#         log_coeff = (
#             np.log(molar_frac * P)
#             - np.log(Z_v - Bw)
#             + (b_i / bw) * (Z_v - 1)
#             - (Aw / (2 * math.sqrt(2) * Bw))
#             * ((2 * avvv / aw) - (b_i / bw))
#             * np.log((Z_v + (1 + math.sqrt(2)) * Bw) / (Z_v + (1 - math.sqrt(2)) * Bw))
#         )
#         # SRK - Peneloux
#         fz_i = np.exp(log_coeff)
#
#         return fz_i, Z_v
#
#     def _calculate_reid_enthalpy_cp(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#     ):
#         Tref = 273.15
#         dT_1 = T - Tref
#         dT_2 = T**2 - Tref**2
#         dT_3 = T**3 - Tref**3
#         dT_4 = T**4 - Tref**4
#         Hid = (
#             self._Cp1_id * (T - Tref)
#             + 0.5 * self._Cp2_id * (T**2 - Tref**2)
#             + self._Cp3_id * (T**3 - Tref**3) / 3
#             + 0.25 * self._Cp4_id * (T**4 - Tref**4)
#         )
#         dYi_dT = (
#             -psi_i
#             * np.power(self._Tkr * T, -0.5)
#             * (1 + psi_i * (1 - np.sqrt(T / self._Tkr)))
#         )
#         d2Yi_dT2 = (
#             0.5 * psi_i * np.power(self._Tkr, -0.5) * (1 + psi_i) * np.power(T, -1.5)
#         )
#         dai_dT = ac_i * dYi_dT
#         d2ai_dT2 = ac_i * d2Yi_dT2
#         R = self._R
#         am = (
#             (1 - BIPs)
#             * np.sqrt(a_i[:, None] * a_i[None, :])
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#         )
#         am = sum(sum(am))  # type: ignore
#         dam_dT = (
#             (1 - BIPs)
#             * np.power(a_i[:, None] * a_i[None, :], -0.5)
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#             * (dai_dT[None, :] * a_i[:, None] + dai_dT[:, None] * a_i[None, :])
#         )
#         dam_dT = 0.5 * sum(sum(dam_dT))  # type: ignore
#         dAm_dT = dam_dT * P / ((R**2) * (T**2))
#         bm = np.dot(molar_frac, np.array(b_i))
#         cm = np.dot(cpen, molar_frac)
#         Am = am * P / (R**2 * T**2)
#         Bm = bm * P / (R * T)
#         Cm = cm * P / (R * T)
#         H_residual = (
#             1000
#             * R
#             * T
#             * (Z - 1 - ((Am - T * dAm_dT) / Bm) * np.log((Z + Cm + Bm) / (Z + Cm)))
#         )
#         enthalpy = H_residual + np.dot(molar_frac, Hid)
#         Cid = (
#             self._Cp1_id
#             + self._Cp2_id * T
#             + self._Cp3_id * T**2
#             + self._Cp4_id * T**3
#         )
#         d2am_dT2 = (
#             molar_frac[None, :]
#             * molar_frac[:, None]
#             * (1 - BIPs)
#             * (
#                 (np.sqrt(a_i[:, None] / a_i[None, :]) * d2ai_dT2[None, :])
#                 + (np.sqrt(a_i[None, :] / a_i[:, None]) * d2ai_dT2[:, None])
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[:, None] / a_i[None, :]))
#                     * dai_dT[:, None]
#                     * (dai_dT[None, :] * a_i[:, None] - a_i[None, :] * dai_dT[:, None])
#                     / a_i[:, None] ** 2
#                 )
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[None, :] / a_i[:, None]))
#                     * dai_dT[None, :]
#                     * (dai_dT[:, None] * a_i[None, :] - a_i[:, None] * dai_dT[None, :])
#                     / a_i[None, :] ** 2
#                 )
#             )
#         )
#         d2am_dT2 = (sum(sum(d2am_dT2))) * 0.5  # type: ignore
#         d2Am_dT2 = d2am_dT2 * P / ((R**2) * (T**2))
#         dCp1 = R * T**2 * d2Am_dT2 * np.log((Z + Cm + Bm) / (Z + Cm)) / Bm
#         dCp2 = (
#             (P / T) * (1 / (Z + Cm - Bm) - T * dAm_dT / ((Z + Cm) * (Z + Cm + Bm)))
#         ) ** 2
#         dCp3 = (P**2 / (R * T)) * (
#             1 / (Z + Cm - Bm) ** 2
#             + (Am / Bm) * (1 / (Z + Cm + Bm) ** 2 - 1 / (Z + Cm) ** 2)
#         )
#         dCp = 1000 * (dCp1 + T * dCp2 / dCp3 - R)
#         dCv = 1000 * (dCp1 - R)
#         Cp = np.dot(molar_frac, Cid) + dCp
#         Cv = np.dot(molar_frac, Cid) + dCv
#         return enthalpy, Cp, Cv
#
#     def _calculate_enthalpy(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
#     ):
#         # Константы рассчитанные
#         if constants == "1":
#             self._Cp1_id = np.array(
#                 [
#                     31.1000000000,
#                     19.8000000000,
#                     19.2500000000,
#                     5.4100000000,
#                     -4.2200000000,
#                     -1.3900000000,
#                     9.4900000000,
#                     -9.5200000000,
#                     -3.6300000000,
#                     2.070882024,
#                     2.680427097,
#                     3.061998984,
#                     3.448838617,
#                     3.963983093,
#                     4.661621733,
#                     5.369993367,
#                     6.246422599,
#                     7.027781774,
#                     7.977955941,
#                     9.402230224,
#                     11.19448541,
#                     17.11742253,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.0135600000,
#                     0.0734300000,
#                     0.0521200000,
#                     0.1781000000,
#                     0.3063000000,
#                     0.3847000000,
#                     0.3313000000,
#                     0.5066000000,
#                     0.4873000000,
#                     0.537931391,
#                     0.57139549,
#                     0.657236586,
#                     0.748346359,
#                     0.874392967,
#                     1.052392918,
#                     1.239796259,
#                     1.475872416,
#                     1.68103,
#                     1.9364032,
#                     2.321552548,
#                     2.795743977,
#                     4.375101534,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     0.0000267900,
#                     -0.0000560200,
#                     0.0000119700,
#                     -0.0000693700,
#                     -0.0001586000,
#                     -0.0001846000,
#                     -0.0001108000,
#                     -0.0002729000,
#                     -0.0002580000,
#                     -0.000224192,
#                     -0.000234127,
#                     -0.000269476,
#                     -0.000307147,
#                     -0.000359429,
#                     -0.000433511,
#                     -0.000511723,
#                     -0.000610382,
#                     -0.000695956,
#                     -0.000802662,
#                     -0.000963668,
#                     -0.001161573,
#                     -0.001821109,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -0.0000000117,
#                     0.0000000172,
#                     -0.0000000113,
#                     0.0000000087,
#                     0.0000000322,
#                     0.0000000290,
#                     -0.0000000028,
#                     0.0000000572,
#                     0.0000000531,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # Константы заказчики
#         elif constants == "2":
#             self._Cp1_id = np.array(
#                 [
#                     31.15,
#                     19.79,
#                     19.25,
#                     5.41,
#                     -4.22,
#                     -1.39,
#                     9.49,
#                     -9.52,
#                     -3.63,
#                     -4.41,
#                     -11.93,
#                     -7.72,
#                     -3.26,
#                     0.86,
#                     4.53,
#                     6.94,
#                     10.03,
#                     12.45,
#                     15.64,
#                     20.74,
#                     26.86,
#                     54.16,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.014,
#                     0.073,
#                     0.052,
#                     0.178,
#                     0.306,
#                     0.385,
#                     0.331,
#                     0.507,
#                     0.487,
#                     0.582,
#                     0.539,
#                     0.606,
#                     0.679,
#                     0.788,
#                     0.948,
#                     1.124,
#                     1.348,
#                     1.538,
#                     1.783,
#                     2.150,
#                     2.598,
#                     4.344,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     2.68e-5,
#                     -5.60e-5,
#                     1.20e-5,
#                     -6.94e-5,
#                     -1.59e-4,
#                     -1.85e-4,
#                     -1.11e-4,
#                     -2.73e-4,
#                     -2.58e-4,
#                     -3.12e-4,
#                     -1.97e-4,
#                     -2.31e-4,
#                     -2.67e-4,
#                     -3.16e-4,
#                     -3.83e-4,
#                     -4.51e-4,
#                     -5.37e-4,
#                     -6.11e-4,
#                     -7.04e-4,
#                     -8.44e-4,
#                     -1.02e-3,
#                     -1.68e-3,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -1.17e-8,
#                     1.72e-8,
#                     -1.13e-8,
#                     8.71e-9,
#                     3.21e-8,
#                     2.90e-8,
#                     -2.82e-9,
#                     5.72e-8,
#                     5.3e-8,
#                     6.49e-8,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # PVTSIM тест
#         elif constants == "3":
#             self._Сp1_id = np.array(
#                 [
#                     3.11e01,
#                     1.98e01,
#                     1.93e01,
#                     5.41e00,
#                     -4.22e00,
#                     -1.39e00,
#                     9.49e00,
#                     -9.52e00,
#                     -3.63e00,
#                     -4.41e00,
#                     -4.95e01,
#                     -6.30e01,
#                     -7.11e01,
#                     -6.05e01,
#                     -2.76e01,
#                     -9.24e00,
#                     4.33e00,
#                     9.35e00,
#                     1.09e01,
#                     1.39e01,
#                     1.69e01,
#                     2.90e01,
#                 ]
#             )
#             self._Сp2_id = np.array(
#                 [
#                     -1.356e-02,
#                     7.343e-02,
#                     5.212e-02,
#                     1.781e-01,
#                     3.063e-01,
#                     3.847e-01,
#                     3.313e-01,
#                     5.066e-01,
#                     4.873e-01,
#                     5.819e-01,
#                     7.058e-01,
#                     8.025e-01,
#                     9.071e-01,
#                     9.931e-01,
#                     1.041e00,
#                     1.149e00,
#                     1.298e00,
#                     1.481e00,
#                     1.723e00,
#                     2.026e00,
#                     2.469e00,
#                     4.452e00,
#                 ]
#             )
#             self._Сp3_id = np.array(
#                 [
#                     2.679e-05,
#                     -5.602e-05,
#                     1.197e-05,
#                     -6.937e-05,
#                     -1.586e-04,
#                     -1.846e-04,
#                     -1.108e-04,
#                     -2.729e-04,
#                     -2.580e-04,
#                     -3.119e-04,
#                     -1.741e-04,
#                     -1.887e-04,
#                     -2.136e-04,
#                     -2.648e-04,
#                     -3.513e-04,
#                     -4.341e-04,
#                     -5.159e-04,
#                     -5.938e-04,
#                     -6.911e-04,
#                     -8.102e-04,
#                     -9.853e-04,
#                     -1.765e-03,
#                 ]
#             )
#             self._Сp4_id = np.array(
#                 [
#                     -1.168e-08,
#                     1.715e-08,
#                     -1.132e-08,
#                     8.712e-09,
#                     3.215e-08,
#                     2.895e-08,
#                     -2.822e-09,
#                     5.723e-08,
#                     5.305e-08,
#                     6.494e-08,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#
#     def vle(self, P, T):
#
#         self._x_i = np.full(self._N, -999, dtype=np.float64)
#         self._y_i = np.full(self._N, -999, dtype=np.float64)
#
#         # Пункт 2 Брусиловский - Уравнение ПР
#         ac_i = 0.457235 * self._R**2 * self._Tkr**2.0 / self._Pkr
#         psi_i = 0.37464 + 1.54226 * self._w - 0.26992 * self._w**2
#         alpha_i = (1 + psi_i * (1 - (T / self._Tkr) ** 0.5)) ** 2
#         a_i = ac_i * alpha_i
#         b_i = 0.077796 * self._R * self._Tkr / self._Pkr
#         # c_i = self._cpen
#         Z_v = 0
#         Z_l = 0
#         W = 0
#
#         K_i = (
#             np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
#         ) ** (1.0)
#         c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
#         fz_i, self._Z_init = self._calculate(P, T, self._z, b_i, c_a_i, 1)
#         m = 0
#         indx = 0
#         eps_f = 1
#         Ri_v = 1
#         TS_v_flag = 0
#         TS_l_flag = 0
#
#         # Поиск газовой фазы
#         while m < self._MaxIterationStable:
#             Yi_v = self._z * K_i
#             Sv = np.sum(Yi_v)
#             self._y_i = np.divide(Yi_v, Sv)
#
#             fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
#             Ri = fz_i / (Sv * fw_i)
#             Ri_v = np.sum((Ri - 1) ** 2)
#
#             if Ri_v < 10 ** (-12):
#                 m = self._MaxIterationStable
#
#             K_i = K_i * Ri
#             TS_v = np.sum(np.log(K_i) ** 2)
#
#             if TS_v < 10 ** (-4):
#                 TS_v_flag = 1
#                 m = self._MaxIterationStable
#             m = m + 1
#
#         K_iv = K_i
#         K_i = (
#             np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
#         ) ** 1.0
#         fz_i, Z_v = self._calculate(P, T, self._z, b_i, c_a_i, 1)
#         ml = 0
#         Ri_l = 1
#
#         # Поиск жидкой фазы
#         while ml < self._MaxIterationStable:
#             Yi_l = self._z / K_i
#             Sl = np.sum(Yi_l)
#             self._x_i = Yi_l / Sl
#
#             fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
#             Ri = Sl * fl_i / fz_i
#             Ri_l = np.sum((Ri - 1) ** 2)
#
#             if Ri_l < 10 ** (-12):
#                 m = self._MaxIterationStable
#
#             K_i = K_i * Ri
#             TS = np.sum(np.log(K_i) ** 2)
#
#             if TS < 10 ** (-4):
#                 TS_l_flag = 1
#                 m = self._MaxIterationStable
#             ml = ml + 1
#
#         K_il = K_i
#         if (
#             (TS_l_flag == 1 and TS_v_flag == 1)
#             or (Sv <= 1 and TS_l_flag == 1)
#             or (Sl <= 1 and TS_v_flag == 1)
#             or (Sv < 1 and Sl <= 1)
#         ):
#             Stable = 1
#             TestPTF = 1
#         else:
#             Stable = 0
#             TestPTF = 0
#         # Flash - расчет
#         if Stable == 0 or TestPTF == 1:
#             Kst_v = np.sum((K_iv - 1) ** 2)
#             Kst_l = np.sum((K_il - 1) ** 2)
#             if Kst_l > Kst_v:
#                 K_i = K_il
#             else:
#                 K_i = K_iv
#
#             m = 0
#             eps_f = 1
#             c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
#             Rr_old = np.zeros(self._N)
#             Crit3 = 0
#             Crit1 = 0
#             Crit2 = 0
#             W_old = 0
#             # while eps_f > 0.0000001 and m < self._MaxIterationFlash:
#             while eps_f > 0.00001:
#                 W = self._findroot(K_i)
#                 self._x_i = self._z / (1 + W * (K_i - 1))
#                 self._y_i = K_i * self._x_i
#
#                 fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
#                 fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
#                 df_lv = np.zeros(self._N)
#
#                 mask = fl_i != 0
#                 Rr = fl_i / fw_i
#                 df_lv = Rr - 1
#                 # if m <= self._N:
#                 #     K_i[mask] *= Rr[mask]
#
#                 # if m > 1:
#                 #     Crit3 = np.sum((Rr - 1) ** 2)
#                 #     Crit1 = Crit3 / np.sum((Rr_old - 1) ** 2)
#                 #     Crit2 = np.abs(W - W_old)
#
#                 # if m > self._N and Crit1 > 0.8 and Crit2 < 0.1 and Crit3 < 0.001:
#                 #     K_i[mask] *= (Rr[mask]) ** 6.0
#
#                 # elif m > self._N:
#                 #     K_i[mask] *= Rr[mask]
#                 # K_i[mask] *= Rr[mask]
#                 K_i *= fl_i / fw_i
#                 W_old = W
#                 Rr_old = Rr
#                 eps_f = np.max(np.abs(df_lv))
#                 m = m + 1
#
#                 # if (m > 5) and (np.abs(W) > 2):
#                 #     m = self._MaxIterationFlash
#                 #     Stable = 1
#                 # else:
#                 #     Stable = 0
#
#                 # if (m > 80) and (np.abs(W - 0.5) > 0.501):
#                 #     m = self._MaxIterationFlash
#                 #     Stable = 1
#                 # else:
#                 #     Stable = 0
#
#         if Stable == 1 or np.abs(W - 0.5) > 0.5:
#             Volume = 1000.0 * (self._Z_init * self._R * T / P)
#             if self._T_pkr < 260:
#                 W = 1
#                 self._x_i = np.zeros(self._N)
#                 self._y_i = self._z
#                 Z_v = self._Z_init
#                 Z_l = 0.0
#             else:
#                 test = Volume * T * T
#                 if Volume * T * T > self._Control_phase:
#                     W = 1
#                     self._x_i = np.zeros(self._N)
#                     self._y_i = self._z
#                     Z_v = self._Z_init
#                     Z_l = 0.0
#                 else:
#                     W = 0
#                     self._x_i = self._z
#                     self._y_i = np.zeros(self._N)
#                     Z_v = 0
#                     Z_l = self._Z_init
#         # Константы
#         # 1 - Константы PVTSim
#         # 2 - Константы Заказчики
#         constants = "3"
#         enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(
#             P,
#             T,
#             self._y_i,
#             a_i,
#             self._c,
#             psi_i,
#             ac_i,
#             b_i,
#             self._cpen,
#             Z_v,
#             constants,
#         )  # type: ignore
#         enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(
#             P,
#             T,
#             self._x_i,
#             a_i,
#             self._c,
#             psi_i,
#             ac_i,
#             b_i,
#             self._cpen,
#             Z_l,
#             constants,
#         )  # type: ignore
#         if W == 1:
#             enthalpy = enthalpy_w
#             Cp = Cp_w
#             Cv = Cv_w
#         elif W == 0:
#             enthalpy = enthalpy_l
#             Cp = Cp_l
#             Cv = Cv_l
#         else:
#             enthalpy = enthalpy_w * (W) + enthalpy_l * (1 - W)
#             Cp = Cp_w * (W) + Cp_l * (1 - W)
#             Cv = Cv_w * (W) + Cv_l * (1 - W)
#
#         # Объем
#         cpen_mix_y = sum(self._cpen * self._y_i)
#         cpen_mix_x = sum(self._cpen * self._x_i)
#         VolumeMy_y = Z_v * self._R * 1000 * T / P - cpen_mix_y
#         VolumeMy_x = Z_l * self._R * 1000 * T / P - cpen_mix_x
#         volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
#         self._Volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
#
#         # Плотность (оставил для отдельной фазы, если пригодится)
#         density_y = sum((self._mass * self._y_i) / VolumeMy_y)
#         density_x = sum((self._mass * self._x_i) / VolumeMy_x)
#         density = sum((self._mass * self._z / self._Volume))
#         self._density = sum((self._mass * self._z / self._Volume))
#
#         return (
#             W,
#             Z_v,
#             Z_l,
#             self._x_i,
#             self._y_i,
#             Stable,
#             m,
#             enthalpy,
#             enthalpy_w,
#             enthalpy_l,
#             Cp,
#             Cp_w,
#             Cp_l,
#             Cv,
#             Cv_w,
#             Cv_l,
#             volume,
#             VolumeMy_y,
#             VolumeMy_x,
#             density,
#             density_y,
#             density_x,
#         )

"""Дохлый скрипт, ответ none"""
# class PR_peneloux_Flash(Flash):
#     def __init__(self, Fluid):
#         super().__init__(Fluid)  # Вызов конструктора родительского класса

#     def _calculate(self, P, T, molar_frac, b_i, c_i, c_a_i, isMax):
#         R = self._R
#         cpen = self._cpen

#         aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
#         bw = np.dot(molar_frac, b_i)
#         Aw = aw * P / (R**2 * T**2)
#         Bw = bw * P / (R * T)
#         cw = np.dot(cpen, molar_frac)

#         Cw = cw * P / (R * T)
#         Biw = b_i * P / (R * T)
#         Ciw = c_i * P / (R * T)

#         # SRK - Peneloux cubic equation coefficients
#         a = 1
        
#         b = Bw + 3 * Cw - 1
        
#         c = Aw - 3 * Bw**2 - 2 * Bw + 2 * Bw * Cw + 3 * Cw**2 - 2 * Cw
#         d = (
#             Aw * Cw
#             - 2 * Bw * Cw
#             - Cw**2
#             + Cw**2 * Bw
#             + Cw**3
#             - 3 * Cw * Bw**2
#             - Aw * Bw
#             + Bw**2
#             + Bw**3
#         )

#         roots = self._solve_cubic(a, b, c, d)
#         Z_v = np.max(roots) if isMax else np.min(roots)

#         avvv = np.dot(molar_frac, c_a_i)

#         BC_matrix = np.outer(Biw, Ciw)
#         CB_matrix = np.outer(Ciw, Biw)

#         hardlog = sum((BC_matrix + CB_matrix) * molar_frac[None, :])

#         # Пенга-Робинсона без поправки
#         log_coeff = (
#             np.log(molar_frac * P)
#             - np.log(Z_v + Cw - Bw)
#             + (Biw - Ciw) / (Z_v + Cw - Bw)
#             - ((Z_v * Biw + 2 * Cw * Biw - hardlog) / (Bw))
#             * ((1 / (Z_v + Cw - Bw)) - 1)
#             - (Aw / (2 * np.sqrt(2) * Bw))
#             * (
#                 ((2 * avvv / aw) - (Biw / Bw))
#                 * np.log(
#                     (
#                         (Z_v + Cw + Bw * (1 + np.sqrt(2)))
#                         / (Z_v + Cw + Bw * (1 - np.sqrt(2)))
#                     )
#                 )
#             )
#         )
#         # SRK - Peneloux
#         fz_i = np.exp(log_coeff)

#         return fz_i, Z_v
    
#     def _calculate_reid_enthalpy_cp(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#     ):
#         Tref = 273.15
#         dT_1 = T - Tref
#         dT_2 = T**2 - Tref**2
#         dT_3 = T**3 - Tref**3
#         dT_4 = T**4 - Tref**4
#         Hid = (
#             self._Cp1_id * (T - Tref)
#             + 0.5 * self._Cp2_id * (T**2 - Tref**2)
#             + self._Cp3_id * (T**3 - Tref**3) / 3
#             + 0.25 * self._Cp4_id * (T**4 - Tref**4)
#         )
#         dYi_dT = (
#             -psi_i
#             * np.power(self._Tkr * T, -0.5)
#             * (1 + psi_i * (1 - np.sqrt(T / self._Tkr)))
#         )
#         d2Yi_dT2 = (
#             0.5 * psi_i * np.power(self._Tkr, -0.5) * (1 + psi_i) * np.power(T, -1.5)
#         )
#         dai_dT = ac_i * dYi_dT
#         d2ai_dT2 = ac_i * d2Yi_dT2
#         R = self._R
#         am = (
#             (1 - BIPs)
#             * np.sqrt(a_i[:, None] * a_i[None, :])
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#         )
#         am = sum(sum(am))  # type: ignore
#         dam_dT = (
#             (1 - BIPs)
#             * np.power(a_i[:, None] * a_i[None, :], -0.5)
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#             * (dai_dT[None, :] * a_i[:, None] + dai_dT[:, None] * a_i[None, :])
#         )
#         dam_dT = 0.5 * sum(sum(dam_dT))  # type: ignore
#         dAm_dT = dam_dT * P / ((R**2) * (T**2))
#         bm = np.dot(molar_frac, np.array(b_i))
#         cm = np.dot(cpen, molar_frac)
#         Am = am * P / (R**2 * T**2)
#         Bm = bm * P / (R * T)
#         Cm = cm * P / (R * T)
#         H_residual = (
#             1000
#             * R
#             * T
#             * (Z - 1 - ((Am - T * dAm_dT) / Bm) * np.log((Z + Cm + Bm) / (Z + Cm)))
#         )
#         enthalpy = H_residual + np.dot(molar_frac, Hid)
#         Cid = (
#             self._Cp1_id
#             + self._Cp2_id * T
#             + self._Cp3_id * T**2
#             + self._Cp4_id * T**3
#         )
#         d2am_dT2 = (
#             molar_frac[None, :]
#             * molar_frac[:, None]
#             * (1 - BIPs)
#             * (
#                 (np.sqrt(a_i[:, None] / a_i[None, :]) * d2ai_dT2[None, :])
#                 + (np.sqrt(a_i[None, :] / a_i[:, None]) * d2ai_dT2[:, None])
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[:, None] / a_i[None, :]))
#                     * dai_dT[:, None]
#                     * (dai_dT[None, :] * a_i[:, None] - a_i[None, :] * dai_dT[:, None])
#                     / a_i[:, None] ** 2
#                 )
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[None, :] / a_i[:, None]))
#                     * dai_dT[None, :]
#                     * (dai_dT[:, None] * a_i[None, :] - a_i[:, None] * dai_dT[None, :])
#                     / a_i[None, :] ** 2
#                 )
#             )
#         )
#         d2am_dT2 = (sum(sum(d2am_dT2))) * 0.5  # type: ignore
#         d2Am_dT2 = d2am_dT2 * P / ((R**2) * (T**2))
#         dCp1 = R * T**2 * d2Am_dT2 * np.log((Z + Cm + Bm) / (Z + Cm)) / Bm
#         dCp2 = (
#             (P / T) * (1 / (Z + Cm - Bm) - T * dAm_dT / ((Z + Cm) * (Z + Cm + Bm)))
#         ) ** 2
#         dCp3 = (P**2 / (R * T)) * (
#             1 / (Z + Cm - Bm) ** 2
#             + (Am / Bm) * (1 / (Z + Cm + Bm) ** 2 - 1 / (Z + Cm) ** 2)
#         )
#         dCp = 1000 * (dCp1 + T * dCp2 / dCp3 - R)
#         dCv = 1000 * (dCp1 - R)
#         Cp = np.dot(molar_frac, Cid) + dCp
#         Cv = np.dot(molar_frac, Cid) + dCv
#         return enthalpy, Cp, Cv
    
#     def _calculate_enthalpy(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
#     ):
#         # Константы рассчитанные
#         if constants == "1":
#             self._Cp1_id = np.array(
#                 [
#                     31.1000000000,
#                     19.8000000000,
#                     19.2500000000,
#                     5.4100000000,
#                     -4.2200000000,
#                     -1.3900000000,
#                     9.4900000000,
#                     -9.5200000000,
#                     -3.6300000000,
#                     2.070882024,
#                     2.680427097,
#                     3.061998984,
#                     3.448838617,
#                     3.963983093,
#                     4.661621733,
#                     5.369993367,
#                     6.246422599,
#                     7.027781774,
#                     7.977955941,
#                     9.402230224,
#                     11.19448541,
#                     17.11742253,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.0135600000,
#                     0.0734300000,
#                     0.0521200000,
#                     0.1781000000,
#                     0.3063000000,
#                     0.3847000000,
#                     0.3313000000,
#                     0.5066000000,
#                     0.4873000000,
#                     0.537931391,
#                     0.57139549,
#                     0.657236586,
#                     0.748346359,
#                     0.874392967,
#                     1.052392918,
#                     1.239796259,
#                     1.475872416,
#                     1.68103,
#                     1.9364032,
#                     2.321552548,
#                     2.795743977,
#                     4.375101534,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     0.0000267900,
#                     -0.0000560200,
#                     0.0000119700,
#                     -0.0000693700,
#                     -0.0001586000,
#                     -0.0001846000,
#                     -0.0001108000,
#                     -0.0002729000,
#                     -0.0002580000,
#                     -0.000224192,
#                     -0.000234127,
#                     -0.000269476,
#                     -0.000307147,
#                     -0.000359429,
#                     -0.000433511,
#                     -0.000511723,
#                     -0.000610382,
#                     -0.000695956,
#                     -0.000802662,
#                     -0.000963668,
#                     -0.001161573,
#                     -0.001821109,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -0.0000000117,
#                     0.0000000172,
#                     -0.0000000113,
#                     0.0000000087,
#                     0.0000000322,
#                     0.0000000290,
#                     -0.0000000028,
#                     0.0000000572,
#                     0.0000000531,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # Константы заказчики
#         elif constants == "2":
#             self._Cp1_id = np.array(
#                 [
#                     31.15,
#                     19.79,
#                     19.25,
#                     5.41,
#                     -4.22,
#                     -1.39,
#                     9.49,
#                     -9.52,
#                     -3.63,
#                     -4.41,
#                     -11.93,
#                     -7.72,
#                     -3.26,
#                     0.86,
#                     4.53,
#                     6.94,
#                     10.03,
#                     12.45,
#                     15.64,
#                     20.74,
#                     26.86,
#                     54.16,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.014,
#                     0.073,
#                     0.052,
#                     0.178,
#                     0.306,
#                     0.385,
#                     0.331,
#                     0.507,
#                     0.487,
#                     0.582,
#                     0.539,
#                     0.606,
#                     0.679,
#                     0.788,
#                     0.948,
#                     1.124,
#                     1.348,
#                     1.538,
#                     1.783,
#                     2.150,
#                     2.598,
#                     4.344,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     2.68e-5,
#                     -5.60e-5,
#                     1.20e-5,
#                     -6.94e-5,
#                     -1.59e-4,
#                     -1.85e-4,
#                     -1.11e-4,
#                     -2.73e-4,
#                     -2.58e-4,
#                     -3.12e-4,
#                     -1.97e-4,
#                     -2.31e-4,
#                     -2.67e-4,
#                     -3.16e-4,
#                     -3.83e-4,
#                     -4.51e-4,
#                     -5.37e-4,
#                     -6.11e-4,
#                     -7.04e-4,
#                     -8.44e-4,
#                     -1.02e-3,
#                     -1.68e-3,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -1.17e-8,
#                     1.72e-8,
#                     -1.13e-8,
#                     8.71e-9,
#                     3.21e-8,
#                     2.90e-8,
#                     -2.82e-9,
#                     5.72e-8,
#                     5.3e-8,
#                     6.49e-8,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # PVTSIM тест
#         elif constants == "3":
#             self._Сp1_id = np.array(
#                 [
#                     3.11e01,
#                     1.98e01,
#                     1.93e01,
#                     5.41e00,
#                     -4.22e00,
#                     -1.39e00,
#                     9.49e00,
#                     -9.52e00,
#                     -3.63e00,
#                     -4.41e00,
#                     -4.95e01,
#                     -6.30e01,
#                     -7.11e01,
#                     -6.05e01,
#                     -2.76e01,
#                     -9.24e00,
#                     4.33e00,
#                     9.35e00,
#                     1.09e01,
#                     1.39e01,
#                     1.69e01,
#                     2.90e01,
#                 ]
#             )
#             self._Сp2_id = np.array(
#                 [
#                     -1.356e-02,
#                     7.343e-02,
#                     5.212e-02,
#                     1.781e-01,
#                     3.063e-01,
#                     3.847e-01,
#                     3.313e-01,
#                     5.066e-01,
#                     4.873e-01,
#                     5.819e-01,
#                     7.058e-01,
#                     8.025e-01,
#                     9.071e-01,
#                     9.931e-01,
#                     1.041e00,
#                     1.149e00,
#                     1.298e00,
#                     1.481e00,
#                     1.723e00,
#                     2.026e00,
#                     2.469e00,
#                     4.452e00,
#                 ]
#             )
#             self._Сp3_id = np.array(
#                 [
#                     2.679e-05,
#                     -5.602e-05,
#                     1.197e-05,
#                     -6.937e-05,
#                     -1.586e-04,
#                     -1.846e-04,
#                     -1.108e-04,
#                     -2.729e-04,
#                     -2.580e-04,
#                     -3.119e-04,
#                     -1.741e-04,
#                     -1.887e-04,
#                     -2.136e-04,
#                     -2.648e-04,
#                     -3.513e-04,
#                     -4.341e-04,
#                     -5.159e-04,
#                     -5.938e-04,
#                     -6.911e-04,
#                     -8.102e-04,
#                     -9.853e-04,
#                     -1.765e-03,
#                 ]
#             )
#             self._Сp4_id = np.array(
#                 [
#                     -1.168e-08,
#                     1.715e-08,
#                     -1.132e-08,
#                     8.712e-09,
#                     3.215e-08,
#                     2.895e-08,
#                     -2.822e-09,
#                     5.723e-08,
#                     5.305e-08,
#                     6.494e-08,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv

#     def vle(self, P, T):
#         self._x_i = np.full(self._N, -999, dtype=np.float64)
#         self._y_i = np.full(self._N, -999, dtype=np.float64)

#         ac_i = 0.42747 * self._R**2 * self._Tkr**2.0 / self._Pkr
#         psi_i = 0.48 + 1.574 * self._w - 0.176 * self._w**2
#         alpha_i = (1 + psi_i * (1 - (T / self._Tkr) ** 0.5)) ** 2
#         a_i = ac_i * alpha_i
#         b_i = 0.08664 * self._R * self._Tkr / self._Pkr
#         c_i = self._cpen
#         Z_v = 0
#         Z_l = 0
#         W = 0

#         K_i = (
#             np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
#         ) ** 1.0
#         c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
#         fz_i, self._Z_init = self._calculate(P, T, self._z, b_i, c_i, c_a_i, True)
#         m = 0
#         TS_v_flag = 0
#         TS_l_flag = 0

#         # Поиск газовой фазы
#         while m < self._MaxIterationStable:
#             Yi_v = self._z * K_i
#             Sv = np.sum(Yi_v)
#             self._y_i = np.divide(Yi_v, Sv)

#             fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_i, c_a_i, True)
#             Ri = fz_i / (Sv * fw_i)
#             Ri_v = np.sum((Ri - 1) ** 2)

#             if Ri_v < 1e-12:
#                 break

#             K_i *= Ri
#             TS_v = np.sum(np.log(K_i) ** 2)

#             if TS_v < 1e-4:
#                 TS_v_flag = 1
#                 break
#             m += 1

#         K_iv = K_i.copy()
#         fz_i, Z_v = self._calculate(P, T, self._z, b_i, c_i, c_a_i, True)
#         m = 0
#         Ri_l = 1

#         # Поиск жидкой фазы
#         while m < self._MaxIterationStable:
#             Yi_l = self._z / K_i
#             Sl = np.sum(Yi_l)
#             self._x_i = Yi_l / Sl

#             fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_i, c_a_i, False)
#             Ri = Sl * fl_i / fz_i
#             Ri_l = np.sum((Ri - 1) ** 2)

#             if Ri_l < 1e-12:
#                 break

#             K_i *= Ri
#             TS = np.sum(np.log(K_i) ** 2)

#             if TS < 1e-4:
#                 TS_l_flag = 1
#                 break
#             m += 1

#         K_il = K_i.copy()
#         if (
#             (TS_l_flag == 1 and TS_v_flag == 1)
#             or (Sv <= 1 and TS_l_flag == 1)
#             or (Sl <= 1 and TS_v_flag == 1)
#             or (Sv < 1 and Sl <= 1)
#         ):
#             Stable = 1
#         else:
#             Stable = 0

#         if Stable == 0:
#             Kst_v = np.sum((K_iv - 1) ** 2)
#             Kst_l = np.sum((K_il - 1) ** 2)
#             K_i = K_il if Kst_l > Kst_v else K_iv

#             m = 0
#             eps_f = 1
#             Rr_old = np.zeros(self._N)
#             W_old = 0
#             while eps_f > 0.000001 and m < self._MaxIterationFlash:
#                 W = self._findroot(K_i)
#                 self._x_i = self._z / (1 + W * (K_i - 1))
#                 self._y_i = K_i * self._x_i

#                 fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_i, c_a_i, True)
#                 fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_i, c_a_i, False)

#                 mask = fl_i != 0
#                 Rr = fl_i / fw_i
#                 df_lv = np.zeros(self._N)
#                 df_lv[mask] = Rr[mask] - 1

#                 if m > 1:
#                     eps_f = np.max(np.abs(df_lv))

#                 if eps_f < 1e-6:
#                     break

#                 if m > 1:
#                     Crit1 = np.sum((Rr - 1) ** 2) / np.sum((Rr_old - 1) ** 2)
#                     Crit2 = np.abs(W - W_old)

#                     if Crit1 > 0.8 and Crit2 < 0.1:
#                         K_i[mask] *= Rr[mask] ** 6.0
#                     else:
#                         K_i[mask] *= Rr[mask]

#                 W_old = W
#                 Rr_old = Rr
#                 m += 1

#         constants = "3"
#         enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(
#             P, T, self._y_i, a_i, self._c, psi_i, ac_i, b_i, c_i, Z_v, constants
#         )
#         enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(
#             P, T, self._x_i, a_i, self._c, psi_i, ac_i, b_i, c_i, Z_l, constants
#         )

#         if W == 1:
#             enthalpy = enthalpy_w
#             Cp = Cp_w
#             Cv = Cv_w
#         elif W == 0:
#             enthalpy = enthalpy_l
#             Cp = Cp_l
#             Cv = Cv_l
#         else:
#             enthalpy = W * enthalpy_w + (1 - W) * enthalpy_l
#             Cp = W * Cp_w + (1 - W) * Cp_l
#             Cv = W * Cv_w + (1 - W) * Cv_l

#         cpen_mix_y = sum(self._cpen * self._y_i)
#         cpen_mix_x = sum(self._cpen * self._x_i)
#         VolumeMy_y = Z_v * self._R * 1000 * T / P - cpen_mix_y
#         VolumeMy_x = Z_l * self._R * 1000 * T / P - cpen_mix_x
#         volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
#         self._Volume = volume

#         density_y = sum((self._mass * self._y_i) / VolumeMy_y)
#         density_x = sum((self._mass * self._x_i) / VolumeMy_x)
#         density = sum((self._mass * self._z) / self._Volume)
#         self._density = density

#         return (
#             W,
#             Z_v,
#             Z_l,
#             self._x_i,
#             self._y_i,
#             Stable,
#             enthalpy,
#             enthalpy_w,
#             enthalpy_l,
#             Cp,
#             Cp_w,
#             Cp_l,
#             Cv,
#             Cv_w,
#             Cv_l,
#             volume,
#             VolumeMy_y,
#             VolumeMy_x,
#             density,
#             density_y,
#             density_x,
#         )
    

# class SRK_Flash(Flash):
#     def __init__(self, Fluid):
#         super().__init__(Fluid)
#
#     def _calculate(self, P, T, molar_frac, b_i, c_a_i, isMax):
#         R = self._R
#
#         aw = np.sum(molar_frac[:, None] * molar_frac[None, :] * c_a_i)
#         bw = np.dot(molar_frac, b_i)
#         Aw = aw * P / (R ** 2 * T ** 2)
#         Bw = bw * P / (R * T)
#
#         # SRK - Peneloux cubic equation coefficients
#         a = 1
#         b = - 1
#         c = (Aw - Bw - Bw**2)
#         d = -(Aw*Bw)
#
#         roots = self._solve_cubic(a, b, c, d)
#         Z_v = np.max(roots) if isMax else np.min(roots)
#
#         avvv = np.dot(molar_frac, c_a_i)
#
#         log_coeff = (
#                 np.log(molar_frac * P)
#                 - np.log(Z_v - Bw)
#                 + (b_i / bw) * (Z_v - 1)
#                 - (Aw / Bw)
#                 * ((2 * avvv / aw) - (b_i / bw))
#                 * np.log(1 + Bw/Z_v)
#         )
#         # SRK - Peneloux
#         fz_i = np.exp(log_coeff)
#
#         return fz_i, Z_v
#
#     def _calculate_reid_enthalpy_cp(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#     ):
#         Tref = 273.15
#         dT_1 = T - Tref
#         dT_2 = T**2 - Tref**2
#         dT_3 = T**3 - Tref**3
#         dT_4 = T**4 - Tref**4
#         Hid = (
#             self._Cp1_id * (T - Tref)
#             + 0.5 * self._Cp2_id * (T**2 - Tref**2)
#             + self._Cp3_id * (T**3 - Tref**3) / 3
#             + 0.25 * self._Cp4_id * (T**4 - Tref**4)
#         )
#         dYi_dT = (
#             -psi_i
#             * np.power(self._Tkr * T, -0.5)
#             * (1 + psi_i * (1 - np.sqrt(T / self._Tkr)))
#         )
#         d2Yi_dT2 = (
#             0.5 * psi_i * np.power(self._Tkr, -0.5) * (1 + psi_i) * np.power(T, -1.5)
#         )
#         dai_dT = ac_i * dYi_dT
#         d2ai_dT2 = ac_i * d2Yi_dT2
#         R = self._R
#         am = (
#             (1 - BIPs)
#             * np.sqrt(a_i[:, None] * a_i[None, :])
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#         )
#         am = sum(sum(am))  # type: ignore
#         dam_dT = (
#             (1 - BIPs)
#             * np.power(a_i[:, None] * a_i[None, :], -0.5)
#             * molar_frac[None, :]
#             * molar_frac[:, None]
#             * (dai_dT[None, :] * a_i[:, None] + dai_dT[:, None] * a_i[None, :])
#         )
#         dam_dT = 0.5 * sum(sum(dam_dT))  # type: ignore
#         dAm_dT = dam_dT * P / ((R**2) * (T**2))
#         bm = np.dot(molar_frac, np.array(b_i))
#         cm = np.dot(cpen, molar_frac)
#         Am = am * P / (R**2 * T**2)
#         Bm = bm * P / (R * T)
#         Cm = cm * P / (R * T)
#         H_residual = (
#             1000
#             * R
#             * T
#             * (Z - 1 - ((Am - T * dAm_dT) / Bm) * np.log((Z + Cm + Bm) / (Z + Cm)))
#         )
#         enthalpy = H_residual + np.dot(molar_frac, Hid)
#         Cid = (
#             self._Cp1_id
#             + self._Cp2_id * T
#             + self._Cp3_id * T**2
#             + self._Cp4_id * T**3
#         )
#         d2am_dT2 = (
#             molar_frac[None, :]
#             * molar_frac[:, None]
#             * (1 - BIPs)
#             * (
#                 (np.sqrt(a_i[:, None] / a_i[None, :]) * d2ai_dT2[None, :])
#                 + (np.sqrt(a_i[None, :] / a_i[:, None]) * d2ai_dT2[:, None])
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[:, None] / a_i[None, :]))
#                     * dai_dT[:, None]
#                     * (dai_dT[None, :] * a_i[:, None] - a_i[None, :] * dai_dT[:, None])
#                     / a_i[:, None] ** 2
#                 )
#                 + (
#                     0.5
#                     * (np.sqrt(a_i[None, :] / a_i[:, None]))
#                     * dai_dT[None, :]
#                     * (dai_dT[:, None] * a_i[None, :] - a_i[:, None] * dai_dT[None, :])
#                     / a_i[None, :] ** 2
#                 )
#             )
#         )
#         d2am_dT2 = (sum(sum(d2am_dT2))) * 0.5  # type: ignore
#         d2Am_dT2 = d2am_dT2 * P / ((R**2) * (T**2))
#         dCp1 = R * T**2 * d2Am_dT2 * np.log((Z + Cm + Bm) / (Z + Cm)) / Bm
#         dCp2 = (
#             (P / T) * (1 / (Z + Cm - Bm) - T * dAm_dT / ((Z + Cm) * (Z + Cm + Bm)))
#         ) ** 2
#         dCp3 = (P**2 / (R * T)) * (
#             1 / (Z + Cm - Bm) ** 2
#             + (Am / Bm) * (1 / (Z + Cm + Bm) ** 2 - 1 / (Z + Cm) ** 2)
#         )
#         dCp = 1000 * (dCp1 + T * dCp2 / dCp3 - R)
#         dCv = 1000 * (dCp1 - R)
#         Cp = np.dot(molar_frac, Cid) + dCp
#         Cv = np.dot(molar_frac, Cid) + dCv
#         return enthalpy, Cp, Cv
#
#     def _calculate_enthalpy(
#         self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
#     ):
#         # Константы рассчитанные
#         if constants == "1":
#             self._Cp1_id = np.array(
#                 [
#                     31.1000000000,
#                     19.8000000000,
#                     19.2500000000,
#                     5.4100000000,
#                     -4.2200000000,
#                     -1.3900000000,
#                     9.4900000000,
#                     -9.5200000000,
#                     -3.6300000000,
#                     2.070882024,
#                     2.680427097,
#                     3.061998984,
#                     3.448838617,
#                     3.963983093,
#                     4.661621733,
#                     5.369993367,
#                     6.246422599,
#                     7.027781774,
#                     7.977955941,
#                     9.402230224,
#                     11.19448541,
#                     17.11742253,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.0135600000,
#                     0.0734300000,
#                     0.0521200000,
#                     0.1781000000,
#                     0.3063000000,
#                     0.3847000000,
#                     0.3313000000,
#                     0.5066000000,
#                     0.4873000000,
#                     0.537931391,
#                     0.57139549,
#                     0.657236586,
#                     0.748346359,
#                     0.874392967,
#                     1.052392918,
#                     1.239796259,
#                     1.475872416,
#                     1.68103,
#                     1.9364032,
#                     2.321552548,
#                     2.795743977,
#                     4.375101534,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     0.0000267900,
#                     -0.0000560200,
#                     0.0000119700,
#                     -0.0000693700,
#                     -0.0001586000,
#                     -0.0001846000,
#                     -0.0001108000,
#                     -0.0002729000,
#                     -0.0002580000,
#                     -0.000224192,
#                     -0.000234127,
#                     -0.000269476,
#                     -0.000307147,
#                     -0.000359429,
#                     -0.000433511,
#                     -0.000511723,
#                     -0.000610382,
#                     -0.000695956,
#                     -0.000802662,
#                     -0.000963668,
#                     -0.001161573,
#                     -0.001821109,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -0.0000000117,
#                     0.0000000172,
#                     -0.0000000113,
#                     0.0000000087,
#                     0.0000000322,
#                     0.0000000290,
#                     -0.0000000028,
#                     0.0000000572,
#                     0.0000000531,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                     0,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # Константы заказчики
#         elif constants == "2":
#             self._Cp1_id = np.array(
#                 [
#                     31.15,
#                     19.79,
#                     19.25,
#                     5.41,
#                     -4.22,
#                     -1.39,
#                     9.49,
#                     -9.52,
#                     -3.63,
#                     -4.41,
#                     -11.93,
#                     -7.72,
#                     -3.26,
#                     0.86,
#                     4.53,
#                     6.94,
#                     10.03,
#                     12.45,
#                     15.64,
#                     20.74,
#                     26.86,
#                     54.16,
#                 ]
#             )
#             self._Cp2_id = np.array(
#                 [
#                     -0.014,
#                     0.073,
#                     0.052,
#                     0.178,
#                     0.306,
#                     0.385,
#                     0.331,
#                     0.507,
#                     0.487,
#                     0.582,
#                     0.539,
#                     0.606,
#                     0.679,
#                     0.788,
#                     0.948,
#                     1.124,
#                     1.348,
#                     1.538,
#                     1.783,
#                     2.150,
#                     2.598,
#                     4.344,
#                 ]
#             )
#             self._Cp3_id = np.array(
#                 [
#                     2.68e-5,
#                     -5.60e-5,
#                     1.20e-5,
#                     -6.94e-5,
#                     -1.59e-4,
#                     -1.85e-4,
#                     -1.11e-4,
#                     -2.73e-4,
#                     -2.58e-4,
#                     -3.12e-4,
#                     -1.97e-4,
#                     -2.31e-4,
#                     -2.67e-4,
#                     -3.16e-4,
#                     -3.83e-4,
#                     -4.51e-4,
#                     -5.37e-4,
#                     -6.11e-4,
#                     -7.04e-4,
#                     -8.44e-4,
#                     -1.02e-3,
#                     -1.68e-3,
#                 ]
#             )
#             self._Cp4_id = np.array(
#                 [
#                     -1.17e-8,
#                     1.72e-8,
#                     -1.13e-8,
#                     8.71e-9,
#                     3.21e-8,
#                     2.90e-8,
#                     -2.82e-9,
#                     5.72e-8,
#                     5.3e-8,
#                     6.49e-8,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                     0.000,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#         # PVTSIM тест
#         elif constants == "3":
#             self._Сp1_id = np.array(
#                 [
#                     3.11e01,
#                     1.98e01,
#                     1.93e01,
#                     5.41e00,
#                     -4.22e00,
#                     -1.39e00,
#                     9.49e00,
#                     -9.52e00,
#                     -3.63e00,
#                     -4.41e00,
#                     -4.95e01,
#                     -6.30e01,
#                     -7.11e01,
#                     -6.05e01,
#                     -2.76e01,
#                     -9.24e00,
#                     4.33e00,
#                     9.35e00,
#                     1.09e01,
#                     1.39e01,
#                     1.69e01,
#                     2.90e01,
#                 ]
#             )
#             self._Сp2_id = np.array(
#                 [
#                     -1.356e-02,
#                     7.343e-02,
#                     5.212e-02,
#                     1.781e-01,
#                     3.063e-01,
#                     3.847e-01,
#                     3.313e-01,
#                     5.066e-01,
#                     4.873e-01,
#                     5.819e-01,
#                     7.058e-01,
#                     8.025e-01,
#                     9.071e-01,
#                     9.931e-01,
#                     1.041e00,
#                     1.149e00,
#                     1.298e00,
#                     1.481e00,
#                     1.723e00,
#                     2.026e00,
#                     2.469e00,
#                     4.452e00,
#                 ]
#             )
#             self._Сp3_id = np.array(
#                 [
#                     2.679e-05,
#                     -5.602e-05,
#                     1.197e-05,
#                     -6.937e-05,
#                     -1.586e-04,
#                     -1.846e-04,
#                     -1.108e-04,
#                     -2.729e-04,
#                     -2.580e-04,
#                     -3.119e-04,
#                     -1.741e-04,
#                     -1.887e-04,
#                     -2.136e-04,
#                     -2.648e-04,
#                     -3.513e-04,
#                     -4.341e-04,
#                     -5.159e-04,
#                     -5.938e-04,
#                     -6.911e-04,
#                     -8.102e-04,
#                     -9.853e-04,
#                     -1.765e-03,
#                 ]
#             )
#             self._Сp4_id = np.array(
#                 [
#                     -1.168e-08,
#                     1.715e-08,
#                     -1.132e-08,
#                     8.712e-09,
#                     3.215e-08,
#                     2.895e-08,
#                     -2.822e-09,
#                     5.723e-08,
#                     5.305e-08,
#                     6.494e-08,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                     0.000e-01,
#                 ]
#             )
#             enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
#                 P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
#             )
#             return enthalpy, Cp, Cv
#
#     def vle(self, P, T):
#         self._x_i = np.full(self._N, -999, dtype=np.float64)
#         self._y_i = np.full(self._N, -999, dtype=np.float64)
#
#         ac_i = 0.42748 * self._R ** 2 * self._Tkr ** 2. / self._Pkr
#         psi_i = 0.48 + 1.574 * self._w - 0.176 * self._w ** 2
#         alpha_i = (1 + psi_i * (1 - (T / self._Tkr) ** .5)) ** 2
#         a_i = ac_i * alpha_i
#         b_i = 0.08664 * self._R * self._Tkr / self._Pkr
#
#         Z_v = 0
#         Z_l = 0
#         W = 0
#
#         K_i = (
#             np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
#         ) ** (1.0)
#         c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
#         fz_i, self._Z_init = self._calculate(P, T, self._z, b_i, c_a_i, 1)
#         m = 0
#         indx = 0
#         eps_f = 1
#         Ri_v = 1
#         TS_v_flag = 0
#         TS_l_flag = 0
#
#         # Vapor phase calculation
#         while m < self._MaxIterationStable:
#             Yi_v = self._z * K_i
#             Sv = np.sum(Yi_v)
#             self._y_i = np.divide(Yi_v, Sv)
#
#             fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
#             Ri = fz_i / (Sv * fw_i)
#             Ri_v = np.sum((Ri - 1) ** 2)
#
#             if Ri_v < 10 ** (-12):
#                 m = self._MaxIterationStable
#
#             K_i = K_i * Ri
#             TS_v = np.sum(np.log(K_i) ** 2)
#
#             if TS_v < 10 ** (-4):
#                 TS_v_flag = 1
#                 m = self._MaxIterationStable
#             m += 1
#
#         K_iv = K_i
#         K_i = (
#             np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
#         ) ** 1.0
#         fz_i, Z_v = self._calculate(P, T, self._z, b_i, c_a_i, 1)
#         ml = 0
#         Ri_l = 1
#
#         # Liquid phase calculation
#         while ml < self._MaxIterationStable:
#             Yi_l = self._z / K_i
#             Sl = np.sum(Yi_l)
#             self._x_i = Yi_l / Sl
#
#             fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
#             Ri = Sl * fl_i / fz_i
#             Ri_l = np.sum((Ri - 1) ** 2)
#
#             if Ri_l < 10 ** (-12):
#                 m = self._MaxIterationStable
#
#             K_i = K_i * Ri
#             TS = np.sum(np.log(K_i) ** 2)
#
#             if TS < 10 ** (-4):
#                 TS_l_flag = 1
#                 m = self._MaxIterationStable
#             ml += 1
#
#         K_il = K_i
#         if (
#             (TS_l_flag == 1 and TS_v_flag == 1)
#             or (Sv <= 1 and TS_l_flag == 1)
#             or (Sl <= 1 and TS_v_flag == 1)
#             or (Sv < 1 and Sl <= 1)
#         ):
#             Stable = 1
#             TestPTF = 1
#         else:
#             Stable = 0
#             TestPTF = 0
#
#         # Flash calculation
#         if Stable == 0 or TestPTF == 1:
#             Kst_v = np.sum((K_iv - 1) ** 2)
#             Kst_l = np.sum((K_il - 1) ** 2)
#             if Kst_l > Kst_v:
#                 K_i = K_il
#             else:
#                 K_i = K_iv
#
#             m = 0
#             eps_f = 1
#             c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
#             Rr_old = np.zeros(self._N)
#             W_old = 0
#             while eps_f > 0.00001:
#                 W = self._findroot(K_i)
#                 self._x_i = self._z / (1 + W * (K_i - 1))
#                 self._y_i = K_i * self._x_i
#
#                 fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
#                 fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
#                 Rr = fl_i / fw_i
#                 K_i *= fl_i / fw_i
#                 W_old = W
#                 eps_f = np.max(np.abs(Rr - 1))
#                 m += 1
#
#         if Stable == 1 or np.abs(W - 0.5) > 0.5:
#             Volume = 1000. * (self._Z_init * self._R * T / P)
#             if self._T_pkr < 260:
#                 W = 1
#                 self._x_i = np.zeros(self._N)
#                 self._y_i = self._z
#                 Z_v = self._Z_init
#                 Z_l = 0.
#             else:
#                 if Volume * T * T > self._Control_phase:
#                     W = 1
#                     self._x_i = np.zeros(self._N)
#                     self._y_i = self._z
#                     Z_v = self._Z_init
#                     Z_l = 0.
#                 else:
#                     W = 0
#                     self._x_i = self._z
#                     self._y_i = np.zeros(self._N)
#                     Z_v = 0
#                     Z_l = self._Z_init
#
#         constants = '3'
#         enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(P, T, self._y_i, a_i, self._c, psi_i, ac_i, b_i,
#                                                         self._cpen, Z_v, constants)
#         enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(P, T, self._x_i, a_i, self._c, psi_i, ac_i, b_i,
#                                                         self._cpen, Z_l, constants)
#         if W == 1:
#             enthalpy = enthalpy_w
#             Cp = Cp_w
#             Cv = Cv_w
#         elif W == 0:
#             enthalpy = enthalpy_l
#             Cp = Cp_l
#             Cv = Cv_l
#         else:
#             enthalpy = enthalpy_w * (W) + enthalpy_l * (1 - W)
#             Cp = Cp_w * (W) + Cp_l * (1 - W)
#             Cv = Cv_w * (W) + Cv_l * (1 - W)
#
#         cpen_mix_y = sum(self._cpen * self._y_i)
#         cpen_mix_x = sum(self._cpen * self._x_i)
#         VolumeMy_y = Z_v * self._R * 1000 * T / P - cpen_mix_y
#         VolumeMy_x = Z_l * self._R * 1000 * T / P - cpen_mix_x
#         volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
#         self._Volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
#
#         density_y = sum((self._mass * self._y_i) / VolumeMy_y)
#         density_x = sum((self._mass * self._x_i) / VolumeMy_x)
#         density = sum((self._mass * self._z / self._Volume))
#         self._density = sum((self._mass * self._z / self._Volume))
#
#         return W, Z_v, Z_l, self._x_i, self._y_i, Stable, m, enthalpy, enthalpy_w, enthalpy_l, Cp, Cp_w, Cp_l, Cv, Cv_w, Cv_l, volume, VolumeMy_y, VolumeMy_x, density, density_y, density_x



class SRK_peneloux_Flash(Flash):
    def __init__(self, Fluid):
        super().__init__(Fluid)
    
    def _calculate_reid_enthalpy_cp(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
    ):
        Tref = 273.15
        dT_1 = T - Tref
        dT_2 = T**2 - Tref**2
        dT_3 = T**3 - Tref**3
        dT_4 = T**4 - Tref**4
        Hid = (
            self._Cp1_id * (T - Tref)
            + 0.5 * self._Cp2_id * (T**2 - Tref**2)
            + self._Cp3_id * (T**3 - Tref**3) / 3
            + 0.25 * self._Cp4_id * (T**4 - Tref**4)
        )
        dYi_dT = (
            -psi_i
            * np.power(self._Tkr * T, -0.5)
            * (1 + psi_i * (1 - np.sqrt(T / self._Tkr)))
        )
        d2Yi_dT2 = (
            0.5 * psi_i * np.power(self._Tkr, -0.5) * (1 + psi_i) * np.power(T, -1.5)
        )
        dai_dT = ac_i * dYi_dT
        d2ai_dT2 = ac_i * d2Yi_dT2
        R = self._R
        am = (
            (1 - BIPs)
            * np.sqrt(a_i[:, None] * a_i[None, :])
            * molar_frac[None, :]
            * molar_frac[:, None]
        )
        am = sum(sum(am))  # type: ignore
        dam_dT = (
            (1 - BIPs)
            * np.power(a_i[:, None] * a_i[None, :], -0.5)
            * molar_frac[None, :]
            * molar_frac[:, None]
            * (dai_dT[None, :] * a_i[:, None] + dai_dT[:, None] * a_i[None, :])
        )
        dam_dT = 0.5 * sum(sum(dam_dT))  # type: ignore
        dAm_dT = dam_dT * P / ((R**2) * (T**2))
        bm = np.dot(molar_frac, np.array(b_i))
        cm = np.dot(cpen, molar_frac)
        Am = am * P / (R**2 * T**2)
        Bm = bm * P / (R * T)
        Cm = cm * P / (R * T)
        H_residual = (
            1000
            * R
            * T
            * (Z - 1 - ((Am - T * dAm_dT) / Bm) * np.log((Z + Cm + Bm) / (Z + Cm)))
        )
        enthalpy = H_residual + np.dot(molar_frac, Hid)
        Cid = (
            self._Cp1_id
            + self._Cp2_id * T
            + self._Cp3_id * T**2
            + self._Cp4_id * T**3
        )
        d2am_dT2 = (
            molar_frac[None, :]
            * molar_frac[:, None]
            * (1 - BIPs)
            * (
                (np.sqrt(a_i[:, None] / a_i[None, :]) * d2ai_dT2[None, :])
                + (np.sqrt(a_i[None, :] / a_i[:, None]) * d2ai_dT2[:, None])
                + (
                    0.5
                    * (np.sqrt(a_i[:, None] / a_i[None, :]))
                    * dai_dT[:, None]
                    * (dai_dT[None, :] * a_i[:, None] - a_i[None, :] * dai_dT[:, None])
                    / a_i[:, None] ** 2
                )
                + (
                    0.5
                    * (np.sqrt(a_i[None, :] / a_i[:, None]))
                    * dai_dT[None, :]
                    * (dai_dT[:, None] * a_i[None, :] - a_i[:, None] * dai_dT[None, :])
                    / a_i[None, :] ** 2
                )
            )
        )
        d2am_dT2 = (sum(sum(d2am_dT2))) * 0.5  # type: ignore
        d2Am_dT2 = d2am_dT2 * P / ((R**2) * (T**2))
        dCp1 = R * T**2 * d2Am_dT2 * np.log((Z + Cm + Bm) / (Z + Cm)) / Bm
        dCp2 = (
            (P / T) * (1 / (Z + Cm - Bm) - T * dAm_dT / ((Z + Cm) * (Z + Cm + Bm)))
        ) ** 2
        dCp3 = (P**2 / (R * T)) * (
            1 / (Z + Cm - Bm) ** 2
            + (Am / Bm) * (1 / (Z + Cm + Bm) ** 2 - 1 / (Z + Cm) ** 2)
        )
        dCp = 1000 * (dCp1 + T * dCp2 / dCp3 - R)
        dCv = 1000 * (dCp1 - R)
        Cp = np.dot(molar_frac, Cid) + dCp
        Cv = np.dot(molar_frac, Cid) + dCv
        return enthalpy, Cp, Cv

    def _calculate_enthalpy(
        self, P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z, constants
    ):
        # Константы рассчитанные
        if constants == "1":
            self._Cp1_id = np.array(
                [
                    31.1000000000,
                    19.8000000000,
                    19.2500000000,
                    5.4100000000,
                    -4.2200000000,
                    -1.3900000000,
                    9.4900000000,
                    -9.5200000000,
                    -3.6300000000,
                    2.070882024,
                    2.680427097,
                    3.061998984,
                    3.448838617,
                    3.963983093,
                    4.661621733,
                    5.369993367,
                    6.246422599,
                    7.027781774,
                    7.977955941,
                    9.402230224,
                    11.19448541,
                    17.11742253,
                ]
            )
            self._Cp2_id = np.array(
                [
                    -0.0135600000,
                    0.0734300000,
                    0.0521200000,
                    0.1781000000,
                    0.3063000000,
                    0.3847000000,
                    0.3313000000,
                    0.5066000000,
                    0.4873000000,
                    0.537931391,
                    0.57139549,
                    0.657236586,
                    0.748346359,
                    0.874392967,
                    1.052392918,
                    1.239796259,
                    1.475872416,
                    1.68103,
                    1.9364032,
                    2.321552548,
                    2.795743977,
                    4.375101534,
                ]
            )
            self._Cp3_id = np.array(
                [
                    0.0000267900,
                    -0.0000560200,
                    0.0000119700,
                    -0.0000693700,
                    -0.0001586000,
                    -0.0001846000,
                    -0.0001108000,
                    -0.0002729000,
                    -0.0002580000,
                    -0.000224192,
                    -0.000234127,
                    -0.000269476,
                    -0.000307147,
                    -0.000359429,
                    -0.000433511,
                    -0.000511723,
                    -0.000610382,
                    -0.000695956,
                    -0.000802662,
                    -0.000963668,
                    -0.001161573,
                    -0.001821109,
                ]
            )
            self._Cp4_id = np.array(
                [
                    -0.0000000117,
                    0.0000000172,
                    -0.0000000113,
                    0.0000000087,
                    0.0000000322,
                    0.0000000290,
                    -0.0000000028,
                    0.0000000572,
                    0.0000000531,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ]
            )
            enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
                P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
            )
            return enthalpy, Cp, Cv
        # Константы заказчики
        elif constants == "2":
            self._Cp1_id = np.array(
                [
                    31.15,
                    19.79,
                    19.25,
                    5.41,
                    -4.22,
                    -1.39,
                    9.49,
                    -9.52,
                    -3.63,
                    -4.41,
                    -11.93,
                    -7.72,
                    -3.26,
                    0.86,
                    4.53,
                    6.94,
                    10.03,
                    12.45,
                    15.64,
                    20.74,
                    26.86,
                    54.16,
                ]
            )
            self._Cp2_id = np.array(
                [
                    -0.014,
                    0.073,
                    0.052,
                    0.178,
                    0.306,
                    0.385,
                    0.331,
                    0.507,
                    0.487,
                    0.582,
                    0.539,
                    0.606,
                    0.679,
                    0.788,
                    0.948,
                    1.124,
                    1.348,
                    1.538,
                    1.783,
                    2.150,
                    2.598,
                    4.344,
                ]
            )
            self._Cp3_id = np.array(
                [
                    2.68e-5,
                    -5.60e-5,
                    1.20e-5,
                    -6.94e-5,
                    -1.59e-4,
                    -1.85e-4,
                    -1.11e-4,
                    -2.73e-4,
                    -2.58e-4,
                    -3.12e-4,
                    -1.97e-4,
                    -2.31e-4,
                    -2.67e-4,
                    -3.16e-4,
                    -3.83e-4,
                    -4.51e-4,
                    -5.37e-4,
                    -6.11e-4,
                    -7.04e-4,
                    -8.44e-4,
                    -1.02e-3,
                    -1.68e-3,
                ]
            )
            self._Cp4_id = np.array(
                [
                    -1.17e-8,
                    1.72e-8,
                    -1.13e-8,
                    8.71e-9,
                    3.21e-8,
                    2.90e-8,
                    -2.82e-9,
                    5.72e-8,
                    5.3e-8,
                    6.49e-8,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                    0.000,
                ]
            )
            enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
                P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
            )
            return enthalpy, Cp, Cv
        # PVTSIM тест
        elif constants == "3":
            self._Сp1_id = np.array(
                [
                    3.11e01,
                    1.98e01,
                    1.93e01,
                    5.41e00,
                    -4.22e00,
                    -1.39e00,
                    9.49e00,
                    -9.52e00,
                    -3.63e00,
                    -4.41e00,
                    -4.95e01,
                    -6.30e01,
                    -7.11e01,
                    -6.05e01,
                    -2.76e01,
                    -9.24e00,
                    4.33e00,
                    9.35e00,
                    1.09e01,
                    1.39e01,
                    1.69e01,
                    2.90e01,
                ]
            )
            self._Сp2_id = np.array(
                [
                    -1.356e-02,
                    7.343e-02,
                    5.212e-02,
                    1.781e-01,
                    3.063e-01,
                    3.847e-01,
                    3.313e-01,
                    5.066e-01,
                    4.873e-01,
                    5.819e-01,
                    7.058e-01,
                    8.025e-01,
                    9.071e-01,
                    9.931e-01,
                    1.041e00,
                    1.149e00,
                    1.298e00,
                    1.481e00,
                    1.723e00,
                    2.026e00,
                    2.469e00,
                    4.452e00,
                ]
            )
            self._Сp3_id = np.array(
                [
                    2.679e-05,
                    -5.602e-05,
                    1.197e-05,
                    -6.937e-05,
                    -1.586e-04,
                    -1.846e-04,
                    -1.108e-04,
                    -2.729e-04,
                    -2.580e-04,
                    -3.119e-04,
                    -1.741e-04,
                    -1.887e-04,
                    -2.136e-04,
                    -2.648e-04,
                    -3.513e-04,
                    -4.341e-04,
                    -5.159e-04,
                    -5.938e-04,
                    -6.911e-04,
                    -8.102e-04,
                    -9.853e-04,
                    -1.765e-03,
                ]
            )
            self._Сp4_id = np.array(
                [
                    -1.168e-08,
                    1.715e-08,
                    -1.132e-08,
                    8.712e-09,
                    3.215e-08,
                    2.895e-08,
                    -2.822e-09,
                    5.723e-08,
                    5.305e-08,
                    6.494e-08,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                    0.000e-01,
                ]
            )
            enthalpy, Cp, Cv = self._calculate_reid_enthalpy_cp(
                P, T, molar_frac, a_i, BIPs, psi_i, ac_i, b_i, cpen, Z
            )
            return enthalpy, Cp, Cv

    def vle(self, P, T):
        self._x_i = np.full(self._N, -999, dtype=np.float64)
        self._y_i = np.full(self._N, -999, dtype=np.float64)

        ac_i = 0.42748 * self._R ** 2 * self._Tkr ** 2. / self._Pkr
        psi_i = 0.48 + 1.574 * self._w - 0.176 * self._w ** 2
        alpha_i = (1 + psi_i * (1 - (T / self._Tkr) ** .5)) ** 2
        a_i = ac_i * alpha_i
        b_i = 0.08664 * self._R * self._Tkr / self._Pkr

        Z_v = 0
        Z_l = 0
        W = 0

        K_i = (
            np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
        ) ** (1.0)
        c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
        fz_i, self._Z_init = self._calculate(P, T, self._z, b_i, c_a_i, 1)
        m = 0
        indx = 0
        eps_f = 1
        Ri_v = 1
        TS_v_flag = 0
        TS_l_flag = 0

        # Vapor phase calculation
        while m < self._MaxIterationStable:
            Yi_v = self._z * K_i
            Sv = np.sum(Yi_v)
            self._y_i = np.divide(Yi_v, Sv)

            fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
            Ri = fz_i / (Sv * fw_i)
            Ri_v = np.sum((Ri - 1) ** 2)

            if Ri_v < 10 ** (-12):
                m = self._MaxIterationStable

            K_i = K_i * Ri
            TS_v = np.sum(np.log(K_i) ** 2)

            if TS_v < 10 ** (-4):
                TS_v_flag = 1
                m = self._MaxIterationStable
            m += 1

        K_iv = K_i
        K_i = (
            np.exp(5.373 * (1 + self._w) * (1 - self._Tkr / T)) * self._Pkr / P
        ) ** 1.0
        fz_i, Z_v = self._calculate(P, T, self._z, b_i, c_a_i, 1)
        ml = 0
        Ri_l = 1

        # Liquid phase calculation
        while ml < self._MaxIterationStable:
            Yi_l = self._z / K_i
            Sl = np.sum(Yi_l)
            self._x_i = Yi_l / Sl

            fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
            Ri = Sl * fl_i / fz_i
            Ri_l = np.sum((Ri - 1) ** 2)

            if Ri_l < 10 ** (-12):
                m = self._MaxIterationStable

            K_i = K_i * Ri
            TS = np.sum(np.log(K_i) ** 2)

            if TS < 10 ** (-4):
                TS_l_flag = 1
                m = self._MaxIterationStable
            ml += 1

        K_il = K_i
        if (
            (TS_l_flag == 1 and TS_v_flag == 1)
            or (Sv <= 1 and TS_l_flag == 1)
            or (Sl <= 1 and TS_v_flag == 1)
            or (Sv < 1 and Sl <= 1)
        ):
            Stable = 1
            TestPTF = 1
        else:
            Stable = 0
            TestPTF = 0

        # Flash calculation
        if Stable == 0 or TestPTF == 1:
            Kst_v = np.sum((K_iv - 1) ** 2)
            Kst_l = np.sum((K_il - 1) ** 2)
            if Kst_l > Kst_v:
                K_i = K_il
            else:
                K_i = K_iv

            m = 0
            eps_f = 1
            c_a_i = (1 - self._c) * np.sqrt(a_i[:, None] * a_i[None, :])
            Rr_old = np.zeros(self._N)
            W_old = 0
            while eps_f > 0.00001:
                W = self._findroot(K_i)
                self._x_i = self._z / (1 + W * (K_i - 1))
                self._y_i = K_i * self._x_i

                fw_i, Z_v = self._calculate(P, T, self._y_i, b_i, c_a_i, 1)
                fl_i, Z_l = self._calculate(P, T, self._x_i, b_i, c_a_i, 0)
                Rr = fl_i / fw_i
                K_i *= fl_i / fw_i
                W_old = W
                eps_f = np.max(np.abs(Rr - 1))
                m += 1

        if Stable == 1 or np.abs(W - 0.5) > 0.5:
            Volume = 1000. * (self._Z_init * self._R * T / P)
            if self._T_pkr < 260:
                W = 1
                self._x_i = np.zeros(self._N)
                self._y_i = self._z
                Z_v = self._Z_init
                Z_l = 0.
            else:
                if Volume * T * T > self._Control_phase:
                    W = 1
                    self._x_i = np.zeros(self._N)
                    self._y_i = self._z
                    Z_v = self._Z_init
                    Z_l = 0.
                else:
                    W = 0
                    self._x_i = self._z
                    self._y_i = np.zeros(self._N)
                    Z_v = 0
                    Z_l = self._Z_init

        constants = '3'
        enthalpy_w, Cp_w, Cv_w = self._calculate_enthalpy(P, T, self._y_i, a_i, self._c, psi_i, ac_i, b_i,
                                                        self._cpen, Z_v, constants)
        enthalpy_l, Cp_l, Cv_l = self._calculate_enthalpy(P, T, self._x_i, a_i, self._c, psi_i, ac_i, b_i,
                                                        self._cpen, Z_l, constants)
        if W == 1:
            enthalpy = enthalpy_w
            Cp = Cp_w
            Cv = Cv_w
        elif W == 0:
            enthalpy = enthalpy_l
            Cp = Cp_l
            Cv = Cv_l
        else:
            enthalpy = enthalpy_w * (W) + enthalpy_l * (1 - W)
            Cp = Cp_w * (W) + Cp_l * (1 - W)
            Cv = Cv_w * (W) + Cv_l * (1 - W)

        cpen_mix_y = sum(self._cpen * self._y_i)
        cpen_mix_x = sum(self._cpen * self._x_i)
        VolumeMy_y = Z_v * self._R * 1000 * T / P - cpen_mix_y
        VolumeMy_x = Z_l * self._R * 1000 * T / P - cpen_mix_x
        volume = VolumeMy_y * W + VolumeMy_x * (1 - W)
        self._Volume = VolumeMy_y * W + VolumeMy_x * (1 - W)

        density_y = sum((self._mass * self._y_i) / VolumeMy_y)
        density_x = sum((self._mass * self._x_i) / VolumeMy_x)
        density = sum((self._mass * self._z / self._Volume))
        self._density = sum((self._mass * self._z / self._Volume))

        return W, Z_v, Z_l, self._x_i, self._y_i, Stable, m, enthalpy, enthalpy_w, enthalpy_l, Cp, Cp_w, Cp_l, Cv, Cv_w, Cv_l, volume, VolumeMy_y, VolumeMy_x, density, density_y, density_x


'''
Примечание

1
Пока что методы _calculate_reid_enthalpy_cp (где считается энтальпия) и _calculate_enthalpy (где задаются константы) не известно как реализовать
Мы сделали расчетную модель для задающихся Cp (поэтому там 3 варианта constants).
На данный момент известно, что значения Cp для C1-C7 могут задаваться заранее из БД (например)
Однако не понятно что насчет тяжелых компонент. Надо будет их считать или нет.

2
Также, у _calculate_enthalpy в разных уравнениях состояния (SRK_Flash, PR_Flash) скорее всего будут свои формулы расчета энтальпии
На данный момент реализовано правильно для SRK_peneloux_Flash

3
Уравнения состояния различаются между собой, наиболее наглядный пример - _calculate в разных уравнениях

4
Проверка стабильности зашита уже в расчетки - см. # Liquid phase calculation и # Поиск жидкой фазы

5
Поговорил с нашим матанщиком, он сказал что разделять проверку стабильности и расчетки - неправильно
Поэтому оставил их тут вшитыми
'''