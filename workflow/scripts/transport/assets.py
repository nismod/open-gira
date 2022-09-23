"""
Logic and data structures for infrastructure assets involved in direct damage.

N.B. This module is imported and used by the Snakefile, so is parsed as part of
initialising the workflow. This means it is executed in a different context to
the traditional 'scripts' (with a different sys.path). If it is ever run as a
script, it will need to be robust to executing with different sys.path entries.
"""

from enum import Enum


class StrEnum(str, Enum):
    """
    Enum with strings as values.
    """
    pass


class Assets(StrEnum):
    """
    Typology of infrastructure assets for direct damage estimation (and more)
    """

    @classmethod
    def values(cls) -> list[str]:
        """
        Return a list of values registered to this enum.
        """
        return list(map(lambda c: c.value, cls))

    @classmethod
    def implemented_assets(cls) -> set[str]:
        """
        Return a set of all asset values implemented in subclasses of Assets
        """
        return {asset for subclass in cls.__subclasses__() for asset in subclass.values()}

    @classmethod
    def valid_selection(cls, requested_assets: set) -> None:
        """
        Check we have functionality to process the requested assets.

        Arguments:
            requested_assets (set): Asset types requested in configuration.

        Raises:
            ValueError if selection is not valid
        """

        if not requested_assets.issubset(cls.implemented_assets()):
            missing_assets = requested_assets - cls.implemented_assets()
            raise ValueError(f"{missing_assets=} are not implemented, please amend request")


class RailAssets(Assets):
    """
    Typology of rail assets
    """
    RAILWAY = "rail_railway",
    BRIDGE = "rail_bridge",
    STATION = "rail_station"


class RoadAssets(Assets):
    """
    Typology of road assets
    """
    UNPAVED = "road_unpaved",
    PAVED = "road_paved",
    BRIDGE = "road_bridge"
