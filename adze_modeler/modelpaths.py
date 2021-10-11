from pathlib import Path
import sys

__all__ = ['ModelDir']

class ModelDir:
    BASE = Path(__file__)
    MEDIA = BASE / "media"
    DATA = BASE / "data"
    RESOURCES = BASE / "resources"
    SNAPSHOTS = BASE / "snapshots"
    DEFAULTS = BASE / "defaults"

    @classmethod
    def set_base(cls, base_):
        base_ = Path(base_)
        
        cls.BASE = base_.parent
        cls.MEDIA = cls.BASE / "media"
        cls.DATA = cls.BASE / "data"
        cls.RESOURCES = cls.BASE / "resources"
        cls.SNAPSHOTS = cls.BASE / "snapshots"
        cls.DEFAULTS = cls.BASE / "defaults"

    @classmethod
    def get_dirs(cls):
        yield cls.BASE
        yield cls.MEDIA
        yield cls.DATA
        yield cls.RESOURCES
        yield cls.SNAPSHOTS
        yield cls.DEFAULTS
        


