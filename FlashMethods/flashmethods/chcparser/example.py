import numpy as np
import matplotlib.pyplot as plt
import PVTSimParser as prs
from time import time
from stable_and_flash import PTFlash
from openpyxl import Workbook


# Инициализация парсера
file_path = "C:/Users/User/PycharmProjects/FlashMethods/FlashMethods/flashmethods/testdomain/mixtures/FLASHtest10.CHC"
extractor = prs.MixtureDataExtractor(file_path)
data_srk = extractor.extract_all_srk()
extractor.print_data_srk()

# Формирование Fluid
Fluid = []
N = len(data_srk["component_name"])
for i in range(N):
    Fluid.append([
        data_srk["component_name"][i],
        data_srk["z"][i],
        data_srk["mass"][i],
        data_srk["Pkr"][i],
        data_srk["Tkr"][i],
        data_srk["Vkr"][i],
        data_srk["w"][i],
        data_srk["cpen_srk"][i],
        data_srk["T_boil"][i],
        data_srk["density_liq_phase"][i]
    ])
    
Ptest = 15
Ttest = 473.15
v = PTFlash(Fluid, BIPs=data_srk["BIPs_SRK"])
(
    W,
    Z_v,
    Z_l,
    x_i,
    y_i,
    Stable,
    m,
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
    volume_y,
    volume_x,
    density,
    density_y,
    density_x,
) = v.vle(Ptest, Ttest)

# Список с названиями переменных
needed_results_names = [
    "Ptest",
    "Ttest",
    "W",
    "Z_v",
    "Z_l",
    "Stable",
    "enthalpy",
    "enthalpy_w",
    "enthalpy_l",
    "Cp",
    "Cp_w",
    "Cp_l",
    "Cv",
    "Cv_w",
    "Cv_l",
    "volume",
    "volume_y",
    "volume_x",
    "density",
    "density_y",
    "density_x",
]

# Автоматический вывод названия переменной и её значения
for name in needed_results_names:
    print(f"{name} = {eval(name)}")


