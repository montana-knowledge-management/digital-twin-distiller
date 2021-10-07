from pathlib import Path

__all__ = ['ModelDir']

class ModelDir:
    BASE = Path(__file__)
    MEDIA = BASE / "media"
    DATA = BASE / "data"
    RESOURCES = BASE / "resources"
    SNAPSHOTS = BASE / "snapshots"

    @classmethod
    def set_base(cls, base_):
        base_ = Path(base_)
        if base_.is_dir():
            cls.BASE = base_
        else:
            cls.BASE = base_.parent

        cls.MEDIA = cls.BASE / "media"
        cls.DATA = cls.BASE / "data"
        cls.RESOURCES = cls.BASE / "resources"
        cls.SNAPSHOTS = cls.BASE / "snapshots"




