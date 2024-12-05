import abc
from contextlib import AbstractContextManager
from typing import Callable

from sqlalchemy.orm import Session


class AbstractRepository(abc.ABC):
    """
    Абстрактный класс репозитория
    """

    pass


class CreateMixin:
    """
    Миксин репозитория, содержащего метод create
    """

    @abc.abstractmethod
    def create(self, *args, **kwargs) -> any:
        """
        Создать записи
        """

        raise NotImplementedError


class RetrieveMixin:
    """
    Миксин репозитория, содержащего метод retrieve
    """

    @abc.abstractmethod
    def retrieve(self, *args, **kwargs) -> any:
        """
        Получить записи
        """

        raise NotImplementedError


class UpdateMixin:
    """
    Миксин репозитория, содержащего метод update
    """

    @abc.abstractmethod
    def update(self, *args, **kwargs) -> any:
        """
        Обновить записи
        """

        raise NotImplementedError


class DeleteMixin:
    """
    Миксин репозитория, содержащего метод delete
    """

    @abc.abstractmethod
    def delete(self, *args, **kwargs) -> any:
        """
        Удалить записи
        """

        raise NotImplementedError


class AbstractAlchemyRepository(AbstractRepository):
    """
    Абстрактный класс репозитория Алхимии
    """

    def __init__(self, session_factory: Callable[..., AbstractContextManager[Session]]) -> None:
        """
        Инициализировать переменные
        :param session_factory: фабрика сессий
        """

        self.session_factory = session_factory
