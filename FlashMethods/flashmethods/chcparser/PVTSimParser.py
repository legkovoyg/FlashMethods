import re
import numpy as np
from abc import ABC, abstractmethod

class BaseParser(ABC):
    def __init__(self, content):
        self.content = content

    @abstractmethod
    def parse(self):
        pass

class ShortCompNameParser(BaseParser):
    def parse(self):
        match = re.search(r"<Short Comp Name>\s*\d+\s+\d+\s*([\s\S]*?)(?=<|$)", self.content)
        if match:
            raw_list = re.findall(r'"([^"]+)"', match.group(1))
            return np.array([item.strip() for item in raw_list])
        else:
            print("Тег <Short Comp Name> не найден в файле.")
            return np.array([])

class ComponentCountParser(BaseParser):
    def parse(self):
        # Извлекаем количество компонентов из тега <Short Comp Name>
        match = re.search(r"<Short Comp Name>\s*(\d+)\s+\d+", self.content)
        if match:
            component_count = int(match.group(1))
            print(f"Извлечено количество компонентов: {component_count}")
            return component_count
        else:
            print("Не удалось найти количество компонентов в теге <Short Comp Name>.")
            return None

class CompAmountParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Comp Amount>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Comp Amount>.")
                return np.array([])
            z = np.array([float(value) for value in raw_values])[2:]
            if len(z) != self.component_count:
                print(f"Мольные доли компонентов имеют длину {len(z)}, ожидается {self.component_count}.")
            return z
        else:
            print("Тег <Comp Amount> не найден в файле.")
            return np.array([])

class MolecularWeightsParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Mwn>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Mwn>.")
                return np.array([])
            mass = np.array([float(value) for value in raw_values])[2:]
            if len(mass) != self.component_count:
                print(f"Молекулярные массы имеют длину {len(mass)}, ожидается {self.component_count}.")
            return mass
        else:
            print("Тег <Mwn> не найден в файле.")
            return np.array([])

class CriticalPressureParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Pc>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Pc>.")
                return np.array([])
            pressures_atm = np.array([float(value) for value in raw_values])[2:]
            if len(pressures_atm) != self.component_count:
                print(f"Критические давления имеют длину {len(pressures_atm)}, ожидается {self.component_count}.")
            pressures_mpa = pressures_atm * 0.101325
            return pressures_mpa
        else:
            print("Тег <Pc> не найден в файле.")
            return np.array([])

class CriticalTemperatureParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Tc>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Tc>.")
                return np.array([])
            temperatures = np.array([float(value) for value in raw_values])[2:]
            if len(temperatures) != self.component_count:
                print(f"Критические температуры имеют длину {len(temperatures)}, ожидается {self.component_count}.")
            return temperatures
        else:
            print("Тег <Tc> не найден в файле.")
            return np.array([])

class CriticalVolumeParser(BaseParser):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)
    
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Vc>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Vc>.")
                return np.array([])
            vc_over_r = np.array([float(value) for value in raw_values])[2:]
            if len(vc_over_r) != self.component_count:
                print(f"Критические объемы имеют длину {len(vc_over_r)}, ожидается {self.component_count}.")
            vc = vc_over_r * self.R * 1000  # Перевод в нужные единицы
            return vc
        else:
            print("Тег <Vc> не найден в файле.")
            return np.array([])

class AcentricFactorParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Acentric Factor>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Acentric Factor>.")
                return np.array([])
            w = np.array([float(value) for value in raw_values])[2:]
            if len(w) != self.component_count:
                print(f"Ацентрические факторы имеют длину {len(w)}, ожидается {self.component_count}.")
            return w
        else:
            print("Тег <Acentric Factor> не найден в файле.")
            return np.array([])

class BoilingTemperatureParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Tb>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Tb>.")
                return np.array([])
            T_boil = np.array([float(value) for value in raw_values])[2:]
            if len(T_boil) != self.component_count:
                print(f"Температуры кипения имеют длину {len(T_boil)}, ожидается {self.component_count}.")
            return T_boil
        else:
            print("Тег <Tb> не найден в файле.")
            return np.array([])

class LiquidDensityParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Density>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Density>.")
                return np.array([])
            density = np.array([float(value) if float(value) > 0 else 0 for value in raw_values])[2:]
            if len(density) != self.component_count:
                print(f"Плотности жидкой фазы имеют длину {len(density)}, ожидается {self.component_count}.")
            return density
        else:
            print("Тег <Density> не найден в файле.")
            return np.array([])

class CpenSRKParser(BaseParser):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)
    
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Cpen SRK>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Cpen SRK>.")
                return np.array([])
            cpen_srk_over_r = np.array([float(value) for value in raw_values])[2:]
            if len(cpen_srk_over_r) != self.component_count:
                print(f"Поправки Пенеле для SRK имеют длину {len(cpen_srk_over_r)}, ожидается {self.component_count}.")
            cpen_srk = cpen_srk_over_r * self.R * 1000  # Перевод в нужные единицы
            return cpen_srk
        else:
            print("Тег <Cpen SRK> не найден в файле.")
            return np.array([])

class CpenPRParser(BaseParser):
    R = 0.082057  # Газовая постоянная в L·atm/(K·mol)
    
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        match = re.search(r"<Cpen PR>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            raw_values = re.findall(r"[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2:
                print("Недостаточно данных для <Cpen PR>.")
                return np.array([])
            cpen_pr_over_r = np.array([float(value) for value in raw_values])[2:]
            if len(cpen_pr_over_r) != self.component_count:
                print(f"Поправки Пенеле для PR имеют длину {len(cpen_pr_over_r)}, ожидается {self.component_count}.")
            cpen_pr = cpen_pr_over_r * self.R * 1000  # Перевод в нужные единицы
            return cpen_pr
        else:
            print("Тег <Cpen PR> не найден в файле.")
            return np.array([])

# Новые парсеры для BIPs
class BIPsSRKParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        # Ищем блок данных внутри тега <kij SRK>
        match = re.search(r"<kij SRK>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            # Извлекаем все числовые значения, включая отрицательные и экспоненциальную запись
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2 + (self.component_count * (self.component_count -1)) // 2:
                print("Недостаточно данных для BIPs_SRK после тега <kij SRK>.")
                return np.zeros((self.component_count, self.component_count), dtype=np.float64)
            # Пропускаем первые два значения
            bips_values = np.array([float(value) for value in raw_values[2:2 + (self.component_count * (self.component_count -1)) // 2]], dtype=np.float64)
            expected_length = (self.component_count * (self.component_count -1)) // 2
            if len(bips_values) != expected_length:
                print(f"BIPs_SRK имеют длину {len(bips_values)}, ожидается {expected_length}.")
                return np.zeros((self.component_count, self.component_count), dtype=np.float64)
            # Создаём симметричную матрицу
            c_matrix = np.zeros((self.component_count, self.component_count), dtype=np.float64)
            k = 0
            for i in range(1, self.component_count):
                for j in range(0, i):
                    c_matrix[i, j] = bips_values[k]
                    c_matrix[j, i] = bips_values[k]
                    k += 1
            return c_matrix
        else:
            print("Тег <kij SRK> не найден в файле.")
            return np.zeros((self.component_count, self.component_count), dtype=np.float64)

class BIPsPRParser(BaseParser):
    def __init__(self, content, component_count):
        super().__init__(content)
        self.component_count = component_count

    def parse(self):
        # Ищем блок данных внутри тега <kij PR>
        match = re.search(r"<kij PR>[\s\S]*?([\d.eE,+\s-]+)\n(?=<|$)", self.content)
        if match:
            # Извлекаем все числовые значения, включая отрицательные и экспоненциальную запись
            raw_values = re.findall(r"-?[\d.eE+-]+", match.group(1))
            if len(raw_values) < 2 + (self.component_count * (self.component_count -1)) // 2:
                print("Недостаточно данных для BIPs_PR после тега <kij PR>.")
                return np.zeros((self.component_count, self.component_count), dtype=np.float64)
            # Пропускаем первые два значения
            bips_values = np.array([float(value) for value in raw_values[2:2 + (self.component_count * (self.component_count -1)) // 2]], dtype=np.float64)
            expected_length = (self.component_count * (self.component_count -1)) // 2
            if len(bips_values) != expected_length:
                print(f"BIPs_PR имеют длину {len(bips_values)}, ожидается {expected_length}.")
                return np.zeros((self.component_count, self.component_count), dtype=np.float64)
            # Создаём симметричную матрицу
            c_matrix = np.zeros((self.component_count, self.component_count), dtype=np.float64)
            k = 0
            for i in range(1, self.component_count):
                for j in range(0, i):
                    c_matrix[i, j] = bips_values[k]
                    c_matrix[j, i] = bips_values[k]
                    k += 1
            return c_matrix
        else:
            print("Тег <kij PR> не найден в файле.")
            return np.zeros((self.component_count, self.component_count), dtype=np.float64)

class MixtureDataExtractor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.content = self._read_file()
        # Извлекаем количество компонентов
        component_count_parser = ComponentCountParser(self.content)
        self.component_count = component_count_parser.parse()
        if self.component_count is None:
            raise ValueError("Не удалось определить количество компонентов.")
        # Размерность матрицы BIPs (для N компонентов, количество BIPs = N*(N-1)/2)
        self.N = self.component_count
        # Общие парсеры для всех данных с учётом количества компонентов
        self.common_parsers = {
            'component_name': ShortCompNameParser(self.content),
            'z': CompAmountParser(self.content, self.component_count),
            'mass': MolecularWeightsParser(self.content, self.component_count),
            'Pkr': CriticalPressureParser(self.content, self.component_count),
            'Tkr': CriticalTemperatureParser(self.content, self.component_count),
            'Vkr': CriticalVolumeParser(self.content, self.component_count),
            'w': AcentricFactorParser(self.content, self.component_count),
            'T_boil': BoilingTemperatureParser(self.content, self.component_count),
            'density_liq_phase': LiquidDensityParser(self.content, self.component_count),
            'cpen_srk': CpenSRKParser(self.content, self.component_count),
            'cpen_pr': CpenPRParser(self.content, self.component_count),
        }
        self.data = {}
        self._c_srk = None  # Матрица BIPs для SRK
        self._c_pr = None   # Матрица BIPs для PR

    def _read_file(self):
        try:
            with open(self.file_path, "r", encoding="utf-8") as file:
                return file.read()
        except FileNotFoundError:
            print(f"Файл {self.file_path} не найден.")
            return ""

    def extract_all_srk(self):
        self.data = {}
        # Извлекаем общие данные (кроме cpen_pr и BIPs_PR)
        for key, parser in self.common_parsers.items():
            if key != 'cpen_pr':  # cpen_pr относится к PR, пропускаем здесь
                self.data[key] = parser.parse()
        # Извлекаем данные для SRK
        self.data['cpen_srk'] = self.common_parsers['cpen_srk'].parse()
        # Извлекаем BIPs для SRK
        bips_srk_parser = BIPsSRKParser(self.content, self.component_count)
        self._c_srk = bips_srk_parser.parse()
        self.data['BIPs_SRK'] = self._c_srk
        return self.data

    def extract_all_pr(self):
        self.data = {}
        # Извлекаем общие данные (кроме cpen_srk и BIPs_SRK)
        for key, parser in self.common_parsers.items():
            if key != 'cpen_srk':  # cpen_srk относится к SRK, пропускаем здесь
                self.data[key] = parser.parse()
        # Извлекаем данные для PR
        self.data['cpen_pr'] = self.common_parsers['cpen_pr'].parse()
        # Извлекаем BIPs для PR
        bips_pr_parser = BIPsPRParser(self.content, self.component_count)
        self._c_pr = bips_pr_parser.parse()
        self.data['BIPs_PR'] = self._c_pr
        return self.data

    def get_c_matrix_srk(self):
        return self._c_srk if self._c_srk is not None else np.zeros((self.N, self.N), dtype=np.float64)

    def get_c_matrix_pr(self):
        return self._c_pr if self._c_pr is not None else np.zeros((self.N, self.N), dtype=np.float64)

    def print_data_srk(self):
        self.print_common_data()
        print("\nПоправка Пенеле (cpen_srk):")
        print(self.data.get('cpen_srk', np.array([])))
        print("\nМатрица BIPs (SRK):")
        print(self.data.get('BIPs_SRK', np.zeros((self.N, self.N), dtype=np.float64)))

    def print_data_pr(self):
        self.print_common_data()
        print("\nПоправка Пенеле для PR (cpen_pr):")
        print(self.data.get('cpen_pr', np.array([])))
        print("\nМатрица BIPs (PR):")
        print(self.data.get('BIPs_PR', np.zeros((self.N, self.N), dtype=np.float64)))

    def print_common_data(self):
        print("Извлечённые компоненты:")
        print(self.data.get('component_name', np.array([])))
        
        print("\nМольные доли исходного состава (z):")
        print(self.data.get('z', np.array([])))
        
        print("\nМолекулярная масса (mass):")
        print(self.data.get('mass', np.array([])))
        
        print("\nКритические давления (Pkr) в МПа:")
        print(self.data.get('Pkr', np.array([])))
        
        print("\nКритические температуры (Tkr):")
        print(self.data.get('Tkr', np.array([])))
        
        print("\nКритические объемы (Vkr):")
        print(self.data.get('Vkr', np.array([])))
        
        print("\nАцентрические факторы (w):")
        print(self.data.get('w', np.array([])))
        
        print("\nТемпературы кипения (T_boil):")
        print(self.data.get('T_boil', np.array([])))
        
        print("\nПлотности жидкой фазы (density_liq_phase):")
        print(self.data.get('density_liq_phase', np.array([])))

if __name__ == "__main__":
    file_path = r"mixtures\FLASHtest3.CHC"
    extractor = MixtureDataExtractor(file_path)
    
    # Извлечение данных для SRK
    print("Извлечение данных для SRK:")
    data_srk = extractor.extract_all_srk()
    extractor.print_data_srk()

    # Извлечение данных для PR
    # print("\nИзвлечение данных для PR:")
    # data_pr = extractor.extract_all_pr()
    # extractor.print_data_pr()
    
    # Пример доступа к матрицам BIPs
    # c_matrix_srk = extractor.get_c_matrix_srk()
    # print("\nМатрица _c для SRK:")
    # print(len(c_matrix_srk))
    # for each_elem in c_matrix_srk:
    #     print(len(each_elem))
    
    # c_matrix_pr = extractor.get_c_matrix_pr()
    # print("\nМатрица _c для PR:")
    # print(len(c_matrix_pr))
    # for each_elem in c_matrix_pr:
    #     print(len(each_elem))
