
from pydantic import Field
from pydantic import BaseModel as Base
from pydantic import ConfigDict


class BaseModel(Base):
    """
    Базовая модель для моделей pydantic
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, populate_by_name=True)

class EquationResultDTO(BaseModel):
    """
    DTO, содержащее результаты расчета по уравнению состояния
    """

    w: float = Field(description="Коэффициент разделения фаз (доля газовой фазы)")
    z_v: float = Field(description="Фактор сжимаемости газовой фазы")
    z_l: float = Field(description="Фактор сжимаемости жидкой фазы")
    x_i: list[float] = Field(description="Мольные доли компонентов в жидкой фазе")
    y_i: list[float] = Field(description="Мольные доли компонентов в газовой фазе")
    is_stable: bool = Field(description="Флаг стабильности")
    m: int = Field(description="Количество итераций, выполненных при расчетах")
    enthalpy: float = Field(description="Энтальпия смеси")
    enthalpy_w: float = Field(description="Энтальпия газовой фазы")
    enthalpy_l: float = Field(description="Энтальпия жидкой фазы")
    cp: float = Field(description="Теплоемкость смеси при постоянном давлении")
    cp_w: float = Field(description="Теплоемкость газовой фазы при постоянном давлении")
    cp_l: float = Field(description="Теплоемкость жидкой фазы при постоянном давлении")
    cv: float = Field(description="Теплоемкость смеси при постоянном объеме")
    cv_w: float = Field(description="Теплоемкость газовой фазы при постоянном объеме")
    cv_l: float = Field(description="Теплоемкость жидкой фазы при постоянном объеме")
    volume: float = Field(description="Удельный объем смеси")
    volume_my_y: float = Field(description="Удельный объем газовой фазы")
    volume_my_x: float = Field(description="Удельный объем жидкой фазы")
    density: float = Field(description="Плотность смеси")
    density_y: float = Field(description="Плотность газовой фазы")
    density_x: float = Field(description="Плотность жидкой фазы")