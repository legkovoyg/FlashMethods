import numpy as np
component_name = np.array(
    # Старый
    # ["N2", "CO2", "C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6", "C7", "C8", "C9", "C10-C11", "C12-C13", "C14-C15",
    # "C16-C18", "C19-C20", "C21-C24", "C25-C29", "C30-C35", "C36-C80"]
    # PVTSim тест Объем
    [
        "N2", "CO2", "C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6", "C7", "C8", "C9", "C10-C11", "C12-C13",
        "C14-C15", "C16-C17", "C18-C20", "C21-C23", "C24-C28", "C29-C34", "C35-C80"]
)

# z - мольные доли исходного состава, %
z = np.array(
    # [0.228, 0.605, 83.578, 7.400, 3.345, 0.755, 0.962, 0.338, 0.316, 0.356, 0.483, 0.593, 0.275, 0.297, 0.163, 0.112,
    #  0.090, 0.035, 0.039, 0.021, 8.00E-3, 1.13E-3]
    # [0.109, 0.574, 52.685, 6.788, 3.872, 1.008, 1.406, 0.575, 0.560, 0.730, 1.569, 2.193, 1.153, 1.484, 1.080, 1.036, 1.338, 0.834, 1.550, 2.200, 2.96, 14.3]
    # [0.486, 0.458, 82.225, 7.269, 4.210, 0.574, 1.599, 0.342, 0.446, 0.310, 0.561, 0.512, 0.248, 0.294, 0.168, 0.110, 0.088, 0.035, 0.037, 0.018, 7.46E-3, 5.87E-4]
    # [0.258, 0.457, 55.057, 6.964, 5.087, 0.808, 2.442, 0.614, 0.829, 0.675, 1.839, 1.933, 1.074, 1.536, 1.192, 1.108, 1.446, 0.927, 1.681, 2.213, 3.08, 8.78]
    # [0.999, 0.001]
    # [0.386, 0.871, 73.444, 7.677, 4.522, 0.944, 1.656, 0.554, 0.665, 0.754, 1.372, 1.596, 0.818, 1.066, 0.814, 0.677,
    #  0.641, 0.345, 0.464, 0.366, 0.259, 0.109]
    # [0.11, 0.672,	83.51,	7.13,	2.906,	0.705,	0.866,	0.33,	0.299,	0.382,	0.52,	0.682,	0.368,	0.461,	0.263,	0.205,	0.211,	0.089,	0.122,	0.088,	0.047,	0.034] #14
    # [0.302,	0.833,	75.697,	7.624,	4.187,	0.903,	1.477,	0.509,	0.577,	0.672,	1.155,	1.381,	0.718,	0.932,	0.676,	0.558,	0.541,	0.28,	0.381,	0.298,	0.201,	0.098] #27
    # [0.328,	0.846,	74.922,	7.646,	4.3,	0.918,	1.537,	0.525,	0.606,	0.7,	1.227,	1.454,	0.753,	0.979,	0.723,	0.599,	0.577,	0.302,	0.41,	0.322,	0.222,	0.104] #28
    # [0.355,	0.858,	74.172,	7.663,	4.412,	0.931,	1.597,	0.54,	0.636,	0.727,	1.3,	1.526,	0.786,	1.024,	0.769,	0.639,	0.61,	0.324,	0.438,	0.345,	0.241,	0.107] #29
    # [0.386, 0.871, 73.444, 7.677, 4.522, 0.944, 1.656, 0.554, 0.665, 0.754, 1.372, 1.596, 0.818, 1.066, 0.814, 0.677,
    #  0.641, 0.345, 0.464, 0.366, 0.259, 0.109]
    
    # PVTSIM тест Объем
    [0.386, 0.871, 73.444, 7.677, 4.522, 0.944, 1.656, 0.554, 0.665, 0.754, 1.372, 1.596, 0.818, 1.066, 0.814, 0.677,
     0.427, 0.559, 0.348, 0.409, 0.289, 0.152]
)

# mass - молекулярная масса
mass = np.array(
    # [16.043, 30.07]
    # [28.014, 44.01, 16.043, 30.07, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 91.925, 105.718, 120.343, 140.56,
    # 169.086, 199.098, 236.892, 269.752, 310.637, 372.292, 448.232, 701.123])
    # PVTSim Тест Объем
    [28.014, 44.010, 16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 96.000, 107.000, 121.000, 140.500,
     168.000, 198.000, 229.500, 262.116, 304.667, 356.427, 433.363, 716.655])


# N - количество компонентов в смеси
N = z.shape[0]
# Pkr - критическое давление
# Pkr = np.array()
Pkr = np.array(
    # [3.394388, 7.376459, 4.600155, 4.883864, 4.245517, 3.647701, 3.799688, 3.384255, 3.374123, 2.968823, 3.574564,
    #  3.122868, 2.771944, 2.50968, 2.19141, 1.972617, 1.792949, 1.694122, 1.573101, 1.488414, 1.424685, 1.351062])
    
    # Хз откуда состав
    # [3.394388, 7.376459, 4.600155, 4.883864, 4.245517, 3.647701, 3.799688, 3.384255, 3.374123,
    # 2.968823, 3.574564, 3.122868, 2.771944, 2.50968, 2.19141, 1.972617, 1.792949, 1.694122,
    # 1.573101, 1.488414, 1.424685, 1.351062])
    
    # PVTSim тест объем
    [33.94, 73.76, 46.00, 48.84, 42.46, 36.48, 38.00, 33.84, 33.74, 29.69, 28.40, 28.25, 27.46, 26.11, 24.50, 23.04,
     21.72, 20.67, 19.68, 18.78, 17.85, 16.24]) / 10
# Tkr - критическая температура
Tkr = np.array(
    # [126.2, 304.19995, 190.6, 305.39996, 369.79996, 408.1001, 425.2001, 460.3999, 469.6, 507.4001,
    #             546.6401, 569.9407, 591.9464, 600.8699, 633.8047, 664.5171, 699.4725, 727.1604, 767.4636,
    #             813.2253, 865.1088, 1039.5869]
    
    # PVTSim тест объем
    [-146.950, 31.050, -82.550, 32.250, 96.650, 134.950, 152.050, 187.250, 196.450, 234.250, 306.000, 342.982, 371.337,
     398.768, 436.074, 468.256, 486.628, 504.100, 528.712, 552.337, 580.831, 661.141]) + 273.15

# Vkr - критический объем
Vkr = np.array(
    # [89.8, 94, 99.00002, 148, 203, 263, 255, 306, 304.0001, 370.0001, 425.4445, 471.2443, 523.4991, 600.5144, 712.9235,
    #  838.2478, 1003.597, 1147.239, 1337.256, 1629.407, 1997.867, 3509.115]
    [89.80, 94.00, 99.00, 148.00, 203.00, 263.00, 255.00, 306.00, 304.00, 370.00, 475.55, 487.84, 529.18, 601.62,
     700.63, 818.96, 959.11, 1107.50, 1301.81, 1551.87, 1930.58, 3751.42]
)

# w - ацентрический фактор
w = np.array(
    # [0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.436319, 0.472542,
    #           0.510256, 0.616047, 0.692018, 0.768374, 0.859284, 0.931732, 0.972036, 1.076324, 1.174966,
    #           1.109354]
    # PVTSim тест объем
    [0.0400, 0.2250, 0.0080, 0.0980, 0.1520, 0.1760, 0.1930, 0.2270, 0.2510, 0.2960, 0.3952, 0.4208, 0.4374, 0.4521,
     0.4764, 0.4981, 0.5075, 0.5175, 0.5344, 0.5501, 0.5687, 0.6169
     ]
    )

# cpen - поправка Пенеле
cpen = np.array(
    # [0.92, 3.03, 0.63, 2.63, 5.06, 7.29, 7.86, 10.93, 12.18, 17.98, 6.72, 13.03, 19.41,
    #              15.12, 19.66, 20.81, 17.47, 10.18, 9.11, -14.21, -49.83, -186.03]
    #     PVTSIM тест объем
    [0.92, 3.03, 0.63, 2.63, 5.06, 7.29, 7.86, 10.93, 12.18, 17.98, 37.69, 36.71, 33.13, 27.54, 22.26, 15.15, 1.02,
     -14.82, -34.82, -63.49, -111.42, -310.16]
    )

# T_boil - температура кипения компонент
T_boil = np.array(
    # [-195.750, -78.500, -161.550, -88.550, -42.050, -11.750, -0.450, 27.850, 36.050, 68.750, 91.950,
    # 116.750, 142.250, 175.324, 218.058, 256.457, 298.203, 330.989, 367.594, 416.965, 464.691, 596.610]
    
    # PVTSIM тест
    [-195.750, -78.500, -161.550, -88.550, -42.050, -11.750, -0.450, 27.850, 36.050, 68.750, 91.950, 116.750, 142.250,
     177.045, 218.194, 256.542, 291.628, 324.667, 363.205, 406.306, 457.057, 609.141])

# density_liq_phase - плотность жидкой фазы компонент (г/см3)
density_liq_phase = np.array(
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6640, 0.763, 0.7738, 0.7833, 0.7953, 0.8101, 0.8223, 0.8351, 0.8464, 0.8579, 0.8728,
    #  0.8877, 0.9256]
    
    # PVTSim тест
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6640, 0.7380, 0.7650, 0.7810, 0.7941, 0.8177, 0.8391, 0.8469, 0.8559, 0.8728, 0.8885,
     0.9076, 0.9606])
