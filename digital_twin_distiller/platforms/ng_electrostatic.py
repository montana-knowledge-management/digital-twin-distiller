from copy import copy

from digital_twin_distiller.platforms.ngsolve import NgSolve


class NgElectrostatics(NgSolve):
    """
    Electrostatic field solver with NgSolve
    """

    def __copy__(self):
        return NgElectrostatics(copy(self.metadata))

    def ng_export_metadata(self):
        self.comment("METADATA", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_material_definitions(self):
        self.comment("MATERIAL DEFINITIONS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_block_labels(self):
        self.comment("BLOCK LABELS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_boundary_definitions(self):
        self.comment("BOUNDARY DEFINITIONS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_solving_steps(self):
        self.comment("SOLVER", 1)
        self.write(f'fes = H1(mesh, order=3, dirichlet="gnd")')

        self.write("u = fes.TrialFunction()")
        self.write("v = fes.TestFunction()")
        self.write("f = LinearForm(fes)")
        self.write("f += 32 * (y*(1-y)+x*(1-x)) * v * dx")

        self.write("a = BilinearForm(fes, symmetric=True)")
        self.write("a += grad(u)*grad(v)*dx")

        self.write("a.Assemble()")
        self.write("f.Assemble()")

        self.write("gfu = GridFunction(fes)")
        self.write('gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec')

        self.write("Draw (gfu)")

        self.newline(2)

    def ng_export_postprocessing(self):
        self.comment("POSTPROCESSING AND EXPORTING", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_closing_steps(self):
        self.comment("CLOSING STEPS", 1)
        self.comment("empty", 1)

        self.newline(2)
