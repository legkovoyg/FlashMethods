import abc


class BaseFactory(abc.ABC):
    """
    Базовый класс фабрики
    """

    @abc.abstractmethod
    def create(self, *args, **kwargs) -> any:
        """
        Создать объект
        """

        raise NotImplementedError
