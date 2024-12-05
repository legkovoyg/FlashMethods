import abc


class Equation(abc.ABC):
    """
    Интерфейс класса уравнения
    """

    @abc.abstractmethod
    def calculate(self, *args, **kwargs) -> any:
        """
        Провести расчет
        """

        raise NotImplementedError
