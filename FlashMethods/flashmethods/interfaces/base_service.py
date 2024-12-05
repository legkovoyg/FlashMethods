import abc


class AbstractService(abc.ABC):
    """
    Абстрактный класс сервиса
    """

    @abc.abstractmethod
    def execute(self, *args, **kwargs) -> any:
        """
        Выполнить логику сервиса
        """

        raise NotImplementedError
