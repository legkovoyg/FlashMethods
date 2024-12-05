from pydantic import BaseModel, PrivateAttr
from typing import List, Optional, Dict
import numpy as np
import pandas as pd


# === DTO для компонента жидкости ===
class FluidComponentDTO(BaseModel):
    """
    Каждый компонент описывается по вот этой штуке
    По идее не все обязательные - например cpen не нужен если используется НЕ peneloux версия уравнений
    T_boil и density_liq_phase могут не использоваться (вроде) если не считается энтальпия, Cp и Cv
    """

    component_name: str
    z: float
    mass: float
    Pkr: float
    Tkr: float
    Vkr: float
    w: float
    cpen: float
    T_boil: float
    density_liq_phase: float


# === DTO для жидкости ===
class FluidDTO(BaseModel):
    """
    Каюсь использовал typing но в Pydantic чето не пошло тут, потом могу переделать
    BIPS могут быть опциональны, из-за них меняется форма фигуры (см. graphs)
    """

    components: List[FluidComponentDTO]
    BIPs: Optional[pd.DataFrame] = None
    _N: int = PrivateAttr(default=0)  # Приватное поле для количества компонентов
    _c: np.ndarray = PrivateAttr(
        default_factory=lambda: np.zeros((0, 0), dtype=np.float64)
    )  # Приватная матрица BIPs
    t: float
    p: float

    class Config:
        arbitrary_types_allowed = True  # Позволяет использовать DataFrame и ndarray

    def __init__(self, **data):
        super().__init__(**data)
        # Установим количество компонентов
        self._N = len(self.components)
        # Преобразуем DataFrame в numpy-матрицу или инициализируем нулевую
        if self.BIPs is not None and not self.BIPs.empty:
            self._c = self.BIPs.to_numpy()
        else:
            self._c = np.zeros((self._N, self._N), dtype=np.float64)

