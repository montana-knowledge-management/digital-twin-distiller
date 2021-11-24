from copy import copy

from digital_twin_distiller.platforms.ngsolve import NgSolve


class NgMagnetostatics(NgSolve):
    """
    Magnetostatic field solver with NgSolve
    """

    def __copy__(self):
        return NgMagnetostatics(copy(self.metadata))
