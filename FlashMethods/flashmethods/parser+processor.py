# === Импорты ===
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Optional, Dict
import re
import numpy as np
import logging
import os
from pathlib import Path

from flashmethods.new_domain.srk_peneloux_flash import SRKPenelouxFlash
from flashmethods.new_domain.pr_flash import PRFlash
from flashmethods.new_domain.srk_flash import SRKFlash


# === Настройка Логирования ===
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

import sys
from pathlib import Path

# Добавляем корневую директорию проекта в PYTHONPATH
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))


# === Вспомогательная функция для определения путей ===
def get_project_paths():
    """
    Определяет основные пути проекта.

    Returns:
        tuple: (project_root, mixtures_path, data_path)
    """
    current_file = Path(__file__).resolve()
    project_root = current_file.parent
    mixtures_path = project_root / "chcparser" / "mixtures"
    data_path = project_root / "data"
    mixtures_path.mkdir(parents=True, exist_ok=True)
    data_path.mkdir(parents=True, exist_ok=True)
    return project_root, mixtures_path, data_path


# === DTO Классы ===
class FluidComponentDTO:
    def __init__(
        self,
        component_name: str,
        z: float,
        mass: float,
        Pkr: float,
        Tkr: float,
        Vkr: float,
        w: float,
        cpen: float,
        T_boil: float,
        density_liq_phase: float,
    ):
        self.component_name = component_name
        self.z = z
        self.mass = mass
        self.Pkr = Pkr
        self.Tkr = Tkr
        self.Vkr = Vkr
        self.w = w
        self.cpen = cpen
        self.T_boil = T_boil
        self.density_liq_phase = density_liq_phase

    def __repr__(self):
        return (
            f"FluidComponentDTO(component_name={self.component_name}, z={self.z}, mass={self.mass}, "
            f"Pkr={self.Pkr}, Tkr={self.Tkr}, Vkr={self.Vkr}, w={self.w}, "
            f"cpen={self.cpen}, T_boil={self.T_boil}, density_liq_phase={self.density_liq_phase})"
        )


class FluidDTO:
    def __init__(
        self,
        components: List[FluidComponentDTO],
        BIPs: Optional[pd.DataFrame],
        p: float,
        t: float,
    ):
        self.components = components
        self.BIPs = BIPs
        self.p = p
        self.t = t


# === Абстрактный парсер ===
class Parser(ABC):
    @abstractmethod
    def load_data(self, file_path: str) -> Optional[Dict]:
        pass


# === Вспомогательные парсеры для CHC ===


class ComponentCountParser:
    def __init__(self, content):
        self.content = content

    def parse(self) -> Optional[int]:
        # Извлекаем количество компонентов из строки <Short Comp Name>
        match = re.search(r"<Short Comp Name>\s*(\d+)\s+\d+", self.content)
        if match:
            component_count = int(match.group(1))
            logger.debug(f"Извлечено количество компонентов: {component_count}")
            return component_count
        else:
            logger.warning(
                "Не удалось найти количество компонентов в теге <Short Comp Name>."
            )
            return None


class BIPsSRKParser:
    def __init__(self, content, component_count):
        self.content = content
        self.component_count = component_count

    def parse(self):
        # Ищем блок данных для SRK kij
        match = re.search(r"<kij SRK>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            expected_length = (self.component_count * (self.component_count - 1)) // 2
            if len(raw_values) < 2 + expected_length:
                logger.warning("Недостаточно данных для BIPs_SRK после тега <kij SRK>.")
                return np.zeros(
                    (self.component_count, self.component_count), dtype=np.float64
                )
            bips_values = np.array(
                [float(value) for value in raw_values[2 : 2 + expected_length]],
                dtype=np.float64,
            )
            if len(bips_values) != expected_length:
                logger.warning(
                    f"BIPs_SRK имеют длину {len(bips_values)}, ожидается {expected_length}."
                )
                return np.zeros(
                    (self.component_count, self.component_count), dtype=np.float64
                )
            c_matrix = np.zeros(
                (self.component_count, self.component_count), dtype=np.float64
            )
            k = 0
            for i in range(1, self.component_count):
                for j in range(0, i):
                    c_matrix[i, j] = bips_values[k]
                    c_matrix[j, i] = bips_values[k]
                    k += 1
            return c_matrix
        else:
            logger.warning("Тег <kij SRK> не найден в файле.")
            return np.zeros(
                (self.component_count, self.component_count), dtype=np.float64
            )


class BIPsPRParser:
    def __init__(self, content, component_count):
        self.content = content
        self.component_count = component_count

    def parse(self):
        # Ищем блок данных для PR kij
        match = re.search(r"<kij PR>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            expected_length = (self.component_count * (self.component_count - 1)) // 2
            if len(raw_values) < 2 + expected_length:
                logger.warning("Недостаточно данных для BIPs_PR после тега <kij PR>.")
                return np.zeros(
                    (self.component_count, self.component_count), dtype=np.float64
                )
            bips_values = np.array(
                [float(value) for value in raw_values[2 : 2 + expected_length]],
                dtype=np.float64,
            )
            if len(bips_values) != expected_length:
                logger.warning(
                    f"BIPs_PR имеют длину {len(bips_values)}, ожидается {expected_length}."
                )
                return np.zeros(
                    (self.component_count, self.component_count), dtype=np.float64
                )
            c_matrix = np.zeros(
                (self.component_count, self.component_count), dtype=np.float64
            )
            k = 0
            for i in range(1, self.component_count):
                for j in range(0, i):
                    c_matrix[i, j] = bips_values[k]
                    c_matrix[j, i] = bips_values[k]
                    k += 1
            return c_matrix
        else:
            logger.warning("Тег <kij PR> не найден в файле.")
            return np.zeros(
                (self.component_count, self.component_count), dtype=np.float64
            )


# === Конкретные Экстракторы CHC ===
class CHCExtractor(ABC):
    def __init__(self, content: str):
        self.content = content

    @abstractmethod
    def extract(self):
        pass


class ShortCompNameExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[str]]:
        match = re.search(r"<Short Comp Name>\s*([\s\S]*?)\n(?=<|$)", self.content)
        if match:
            raw_list = re.findall(r'"([^"]+)"', match.group(1))
            component_names = [item.strip() for item in raw_list]
            logger.info(
                f"ShortCompNameExtractor: успешно извлечены {len(component_names)} component_name"
            )
            return {"component_name": component_names}
        else:
            logger.warning("Тег <Short Comp Name> не найден в файле.")
            return {"component_name": []}


class CompAmountExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Comp Amount>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            z = [float(value) for value in raw_values][2:]
            logger.info(f"CompAmountExtractor: успешно извлечены {len(z)}z")
            return {"z": z}
        else:
            logger.warning("Тег <Comp Amount> не найден в файле.")
            return {"z": []}


class MolecularWeightsExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Mwn>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            mass = [float(value) for value in raw_values][2:]
            logger.info(
                f"MolecularWeightsExtractor: успешно извлечены {len(mass)} mass"
            )
            return {"mass": mass}
        else:
            logger.warning("Тег <Mwn> не найден в файле.")
            return {"mass": []}


class CriticalPressureExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Pc>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            pressures_atm = [float(value) for value in raw_values][2:]
            pressures_mpa = [value * 0.101325 for value in pressures_atm]
            logger.info(
                f"CriticalPressureExtractor: успешно извлечены {len(pressures_mpa)} Pkr"
            )
            return {"Pkr": pressures_mpa}
        else:
            logger.warning("Тег <Pc> не найден в файле.")
            return {"Pkr": []}


class CriticalTemperatureExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Tc>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            temperatures = [float(value) for value in raw_values][2:]
            logger.info(
                f"CriticalTemperatureExtractor: успешно извлечены {len(temperatures)} Tkr"
            )
            return {"Tkr": temperatures}
        else:
            logger.warning("Тег <Tc> не найден в файле.")
            return {"Tkr": []}


class CriticalVolumeExtractor(CHCExtractor):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)

    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Vc>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            vc_over_r = [float(value) for value in raw_values][2:]
            vc = [value * self.R * 1000 for value in vc_over_r]
            logger.info(f"CriticalVolumeExtractor: успешно извлечены {len(vc)} Vkr")
            return {"Vkr": vc}
        else:
            logger.warning("Тег <Vc> не найден в файле.")
            return {"Vkr": []}


class AcentricFactorExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(
            r"<Acentric Factor>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content
        )
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            w = [float(value) for value in raw_values][2:]
            logger.info(f"AcentricFactorExtractor: успешно извлечены {len(w)} w")
            return {"w": w}
        else:
            logger.warning("Тег <Acentric Factor> не найден в файле.")
            return {"w": []}


class BoilingTemperatureExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Tb>[\s\S]*?([\d.,\s]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.]+", match.group(1))
            T_boil = [float(value) for value in raw_values][2:]
            logger.info(
                f"BoilingTemperatureExtractor: успешно извлечены {len(T_boil)} T_boil"
            )
            return {"T_boil": T_boil}
        else:
            logger.warning("Тег <Tb> не найден в файле.")
            return {"T_boil": []}


class LiquidDensityExtractor(CHCExtractor):
    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Density>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            density = [float(value) if float(value) > 0 else 0 for value in raw_values][
                2:
            ]
            logger.info(
                f"LiquidDensityExtractor: успешно извлечены {len(density)} density_liq_phase"
            )
            return {"density_liq_phase": density}
        else:
            logger.warning("Тег <Density> не найден в файле.")
            return {"density_liq_phase": []}


class CpenSRKExtractor(CHCExtractor):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)

    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Cpen SRK>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            cpen_srk_over_r = [float(value) for value in raw_values][2:]
            cpen = [value * self.R * 1000 for value in cpen_srk_over_r]
            logger.info(f"CpenSRKExtractor: успешно извлечены {len(cpen)} cpen")
            return {"cpen": cpen}
        else:
            logger.warning("Тег <Cpen SRK> не найден в файле.")
            return {"cpen": []}


class CpenPRExtractor(CHCExtractor):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)

    def extract(self) -> Dict[str, List[float]]:
        match = re.search(r"<Cpen PR>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            cpen_pr_over_r = [float(value) for value in raw_values][2:]
            cpen_pr = [value * self.R * 1000 for value in cpen_pr_over_r]
            logger.info(f"CpenPRExtractor: успешно извлечены {len(cpen_pr)} cpen")
            return {"cpen": cpen_pr}
        else:
            logger.warning("Тег <Cpen PR> не найден в файле.")
            return {"cpen": []}


# === Парсеры PVTSim для SRK и PR ===
class PVTSimParserSRK(Parser):
    def load_data(self, file_path: str) -> Optional[Dict]:
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                content = file.read()
        except FileNotFoundError:
            logger.error(f"{self.__class__.__name__}: файл {file_path} не найден.")
            return None
        except Exception as e:
            logger.error(
                f"{self.__class__.__name__}: ошибка при чтении файла {file_path}: {e}"
            )
            return None

        # Определяем число компонентов
        component_count_parser = ComponentCountParser(content)
        component_count = component_count_parser.parse()
        if component_count is None:
            return None

        extractors = [
            ShortCompNameExtractor(content),
            CompAmountExtractor(content),
            MolecularWeightsExtractor(content),
            CriticalPressureExtractor(content),
            CriticalTemperatureExtractor(content),
            CriticalVolumeExtractor(content),
            AcentricFactorExtractor(content),
            BoilingTemperatureExtractor(content),
            LiquidDensityExtractor(content),
            CpenSRKExtractor(content),
        ]

        data = {}
        for extractor in extractors:
            extracted_data = extractor.extract()
            if isinstance(extracted_data, dict):
                data.update(extracted_data)

        # Извлекаем BIPs для SRK
        bips_srk_parser = BIPsSRKParser(content, component_count)
        c_matrix_srk = bips_srk_parser.parse()
        # Преобразуем в DataFrame
        if c_matrix_srk.size > 0:
            data["BIPs"] = pd.DataFrame(c_matrix_srk)
        else:
            data["BIPs"] = None

        return data


class PVTSimParserPR(Parser):
    def load_data(self, file_path: str) -> Optional[Dict]:
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                content = file.read()
        except FileNotFoundError:
            logger.error(f"{self.__class__.__name__}: файл {file_path} не найден.")
            return None
        except Exception as e:
            logger.error(
                f"{self.__class__.__name__}: ошибка при чтении файла {file_path}: {e}"
            )
            return None

        # Определяем число компонентов
        component_count_parser = ComponentCountParser(content)
        component_count = component_count_parser.parse()
        if component_count is None:
            return None

        extractors = [
            ShortCompNameExtractor(content),
            CompAmountExtractor(content),
            MolecularWeightsExtractor(content),
            CriticalPressureExtractor(content),
            CriticalTemperatureExtractor(content),
            CriticalVolumeExtractor(content),
            AcentricFactorExtractor(content),
            BoilingTemperatureExtractor(content),
            LiquidDensityExtractor(content),
            CpenPRExtractor(content),
        ]

        data = {}
        for extractor in extractors:
            extracted_data = extractor.extract()
            if isinstance(extracted_data, dict):
                data.update(extracted_data)

        # Извлекаем BIPs для PR
        bips_pr_parser = BIPsPRParser(content, component_count)
        c_matrix_pr = bips_pr_parser.parse()
        # Преобразуем в DataFrame
        if c_matrix_pr.size > 0:
            data["BIPs"] = pd.DataFrame(c_matrix_pr)
        else:
            data["BIPs"] = None

        return data


# === Маппер для данных ===
class FluidMapper:
    @staticmethod
    def map_to_fluid_dto(
        data: Dict, p: float, t: float, bips: Optional[pd.DataFrame] = None
    ) -> FluidDTO:
        component_names = data.get("component_name", [])
        z = data.get("z", [])
        mass = data.get("mass", [])
        Pkr = data.get("Pkr", [])
        Tkr = data.get("Tkr", [])
        Vkr = data.get("Vkr", [])
        w = data.get("w", [])
        cpen = data.get("cpen", [])
        T_boil = data.get("T_boil", [])
        density_liq_phase = data.get("density_liq_phase", [])

        num_components = len(component_names)
        # Расширение массива 'w' ведущими нулями, если его длина меньше количества компонентов
        if len(w) < num_components:
            padding_length = num_components - len(w)
            padding = [0.0] * padding_length
            w = padding + w
            logger.debug(
                f"Массив 'w' был расширен {padding_length} ведущими нулями: {w}"
            )
        # Проверка согласованности данных
        for key, value in [
            ("z", z),
            ("mass", mass),
            ("Pkr", Pkr),
            ("Tkr", Tkr),
            ("Vkr", Vkr),
            ("cpen", cpen),
            ("T_boil", T_boil),
            ("density_liq_phase", density_liq_phase),
        ]:
            if len(value) != num_components:
                raise ValueError(
                    f"Длина данных для {key} ({len(value)}) не совпадает с количеством компонентов ({num_components})."
                )

        components = [
            FluidComponentDTO(
                component_name=component_names[i],
                z=float(z[i]),
                mass=float(mass[i]),
                Pkr=float(Pkr[i]),
                Tkr=float(Tkr[i]),
                Vkr=float(Vkr[i]),
                w=float(w[i]),
                cpen=float(cpen[i]),
                T_boil=float(T_boil[i]),
                density_liq_phase=float(density_liq_phase[i]),
            )
            for i in range(num_components)
        ]
        logger.info("FluidMapper: успешно преобразованы данные в FluidDTO")
        return FluidDTO(components=components, BIPs=bips, p=p, t=t)


# === Универсальный обработчик данных ===
class FluidDataProcessor:
    def __init__(self, parser: Parser):
        self.parser = parser

    def process(self, file_path: str, p: float, t: float) -> FluidDTO:
        # Загрузка данных
        data = self.parser.load_data(file_path)
        if not data:
            raise ValueError("Не удалось загрузить данные.")

        bips = data.get("BIPs")

        # Маппинг данных в FluidDTO
        return FluidMapper.map_to_fluid_dto(data, p, t, bips)


# === Основной блок ===
if __name__ == "__main__":
    # Получаем пути
    project_root, mixtures_path, data_path = get_project_paths()

    # Формируем полный путь к CHC-файлу
    chc_file_path = mixtures_path / "FLASHtest4.CHC"

    # Проверяем существование файла CHC
    if not chc_file_path.exists():
        logger.error(f"Файл {chc_file_path} не существует!")
        logger.info(f"Текущая директория: {os.getcwd()}")
        logger.info(f"Ожидаемый путь к файлу: {chc_file_path}")
        raise FileNotFoundError(f"Файл {chc_file_path} не найден")

    # Пример использования PVTSimParserSRK
    pvtsim_srk_parser = PVTSimParserSRK()
    processor_srk = FluidDataProcessor(parser=pvtsim_srk_parser)

    # try:
    # Процессинг данных из CHC
    fluid_srk = processor_srk.process(file_path=str(chc_file_path), p=15, t=323.15)
    print(fluid_srk.components)
    # Передаем DTO в расчетный метод
    calculated_PR = PRFlash(fluid_srk)
    pr = calculated_PR.calculate()

    calculated_SRK = SRKFlash(fluid_srk)
    srk = calculated_SRK.calculate()

    calculated_SRK_Peneloux = SRKPenelouxFlash(fluid_srk)
    srk_peneloux_srk = calculated_SRK_Peneloux.calculate()

    print(f"W (SRK_Peneloux): {srk_peneloux_srk}\n")
    print(f"W (PR): {pr}\n")
    print(f"W (SRK): {srk}\n")

    # except ValueError as ve:
    #     logger.error(f"Ошибка при обработке CHC данных (SRK): {ve}")
    # except Exception as e:
    #     logger.error(f"Неожиданная ошибка: {e}")
