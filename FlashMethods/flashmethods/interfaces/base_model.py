from pydantic import BaseModel as Base
from pydantic import ConfigDict


class BaseModel(Base):
    """
    Базовая модель для моделей pydantic
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, populate_by_name=True)
