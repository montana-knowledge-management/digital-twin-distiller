from pathlib import Path

__all__ = ["ModelDir"]


class ModelDir:
    BASE = Path(__file__).parent
    MEDIA = BASE / "media"
    DATA = BASE / "data"
    RESOURCES = BASE / "resources"
    SNAPSHOTS = BASE / "snapshots"
    DEFAULTS = BASE / "defaults"
    DOCS = BASE / "docs"

    @classmethod
    def set_base(cls, base_):
        base_ = Path(base_).resolve()

        if base_.suffix:
            base_ = base_.parent

        cls.BASE = base_
        cls.MEDIA = cls.BASE / "media"
        cls.DATA = cls.BASE / "data"
        cls.RESOURCES = cls.BASE / "resources"
        cls.SNAPSHOTS = cls.BASE / "snapshots"
        cls.DEFAULTS = cls.BASE / "defaults"
        cls.DOCS = cls.BASE / "docs"

    @classmethod
    def get_dirs(cls):
        yield cls.BASE
        yield cls.MEDIA
        yield cls.DATA
        yield cls.RESOURCES
        yield cls.SNAPSHOTS
        yield cls.DEFAULTS
        yield cls.DOCS
