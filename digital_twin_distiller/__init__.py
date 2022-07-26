from digital_twin_distiller.utils import *
from digital_twin_distiller.objects import CircleArc, Line, Node, Rectangle
from digital_twin_distiller.modelpaths import ModelDir

from digital_twin_distiller.boundaries import (
    AntiPeriodicAirGap,
    AntiPeriodicBoundaryCondition,
    DirichletBoundaryCondition,
    NeumannBoundaryCondition,
    PeriodicAirGap,
    PeriodicBoundaryCondition,
)
from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata, Agros2DMetadata, NgSolveMetadata, NgElectrostaticMetadata
from digital_twin_distiller.modelpiece import ModelPiece
from digital_twin_distiller.platforms import *
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller.model import BaseModel

__all__ = ["new"]
