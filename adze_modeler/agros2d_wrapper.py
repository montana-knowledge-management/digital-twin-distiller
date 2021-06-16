"""
to execture a script use:
    agros2d_solver -s script.py

execute a script with gui:
    agros2d -s script.py
"""


class Agros2DWrapper:
    def __init__(self):
        self.problem_type = None
        self.coordinate_type = "planar"

    def set_problem_type(self, ptype):
        """
        This function is used to set up the type of the problem:
        Supported problems:
            - Electrostatics
            - Magnetics

        :param ptype: Type of the problem. Possible values: 'electrostatic', 'magnetic'
        """

        if ptype.lower() in {"electrostatic", "magnetic"}:
            self.problem_type = ptype.lower()
        else:
            raise NotImplementedError(f"This type of problem ('{ptype}') is not implemented. Typo maybe?")

    def set_coordinate_type(self, coordinate_type):
        """ """
        if coordinate_type.lower() in {"planar", "axisymmetric"}:
            self.coordinate_type = coordinate_type.lower()
        else:
            raise NotImplementedError(
                f"This type of coordinate type ('{coordinate_type}') is not implemented. " f"Typo maybe?"
            )


class Agros2DExecutor:
    def __init__(self):
        pass


if __name__ == "__main__":
    ag = Agros2DWrapper()
    ag.set_problem_type("alma")
