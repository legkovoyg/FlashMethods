import abc


class BaseAlchemySession(abc.ABC):
    """
    Базовый класс сессии Алхимии
    """

    @abc.abstractmethod
    def SessionFactory(self):  # noqa
        raise NotImplementedError

    @abc.abstractmethod
    def get_session(self):
        raise NotImplementedError

    @abc.abstractmethod
    def _build_engine(self):
        raise NotImplementedError
