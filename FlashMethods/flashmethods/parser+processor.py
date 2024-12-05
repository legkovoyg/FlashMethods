import pandas as pd
import json
from openpyxl import load_workbook
from openpyxl.utils.exceptions import InvalidFileException
from abc import ABC, abstractmethod
import numpy as np
from typing import List, Optional, Dict
from dtos.fluid_dto import FluidComponentDTO, FluidDTO

from new_domain.pr_flash import PRFlash
from domain.srk_flash import SRKFlash
from domain.srk_peneloux_flash import SRK_peneloux_Flash


# === Абстрактный парсер ===
class Parser(ABC):
    @abstractmethod
    def load_data(self, file_path: str):
        pass


# === Парсер данных из Excel ===
class ExcelDataParser(Parser):
    def load_data(self, file_path: str) -> Dict:
        try:
            wb = load_workbook(file_path)
            ws = wb.active
            data = {}
            # Итерация по колонкам
            for col_idx, col in enumerate(ws.iter_cols(values_only=True), start=1):
                header = ws.cell(row=1, column=col_idx).value
                if header is not None:
                    # Фильтрация данных, исключая None и заголовок
                    col_data = [cell for cell in col[1:] if cell is not None]
                    data[header] = col_data
            return data
        except InvalidFileException:
            print("ExcelDataParser: ошибка при открытии файла")
            return None


# === Парсер матрицы BIPs из Excel ===
class ExcelBIPSParser(Parser):
    def load_data(self, file_path: str) -> pd.DataFrame:
        try:
            BIPS = pd.read_excel(file_path, header=None)
            return BIPS
        except InvalidFileException:
            print("ExcelBIPSParser: ошибка при открытии файла")
            return None


# === Маппер для данных ===
class FluidMapper:
    @staticmethod
    def map_to_fluid_dto(data: Dict, p: float, t:float, bips: Optional[pd.DataFrame] = None) -> FluidDTO:
        components = [
            FluidComponentDTO(
                component_name=data["component_name"][i],
                z=float(data["z"][i]),
                mass=float(data["mass"][i]),
                Pkr=float(data["Pkr"][i]),
                Tkr=float(data["Tkr"][i]),
                Vkr=float(data["Vkr"][i]),
                w=float(data["w"][i]),
                cpen=float(data["cpen"][i]),
                T_boil=float(data["T_boil"][i]),
                density_liq_phase=float(data["density_liq_phase"][i]),
            )
            for i in range(len(data["component_name"]))
        ]
        return FluidDTO(components=components, BIPs=bips, p = p, t = t)


# === Универсальный обработчик данных ===
class FluidDataProcessor:
    def __init__(self, parser: Parser):
        self.parser = parser

    def process(self, file_path: str, p:float, t:float, bips_path: Optional[str] = None) -> FluidDTO:
        # Загрузка данных
        data = self.parser.load_data(file_path)
        if not data:
            raise ValueError("Не удалось загрузить данные.")

        # Загрузка BIPs, если указан путь
        bips = None
        if bips_path:
            bips_parser = ExcelBIPSParser()
            bips = bips_parser.load_data(bips_path)

        # Маппинг данных в FluidDTO
        return FluidMapper.map_to_fluid_dto(data, p, t, bips)


# === Основной блок ===
if __name__ == "__main__":
    # обработчик
    excel_parser = ExcelDataParser()
    processor = FluidDataProcessor(parser=excel_parser)
    # Присваиваем в переменную объект dto
    fluid = processor.process(
        file_path="data/input_data.xlsx",
        p = 30.1,
        t = 250.1,  
        bips_path="data/BIPSS.xlsx"
    )
    # Кидаем в нужный класс объект dto (пользователь сам выбирает уравнение состояния)
    calculated_PR = PRFlash(fluid)
    calculated_SRK = SRKFlash(fluid)
    calculated_SRK_Peneloux = SRK_peneloux_Flash(fluid)
    # Такой вот результат расчета
    pr = calculated_PR.calculate()
    srk = calculated_SRK.calculate()
    srk_peneloux = calculated_SRK_Peneloux.calculate()
    print(pr)

 