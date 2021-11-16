from digital_twin_distiller import Material
from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.platforms.platform import Platform

class NgElectrostatics(Platform):
    def __copy__(self):
        pass

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"# {str_}\n")
        self.newline(nb_newline)

    def export_preamble(self):
        ...

    def export_metadata(self):
        ...

    def export_material_definition(self, mat: Material):
        ...

    def export_block_label(self, x, y, mat: Material):
        ...

    def export_boundary_definition(self, b: BoundaryCondition):
        ...

    def export_geometry_element(self, e, boundary=None):
        ...

    def export_solving_steps(self):
        ...

    def export_results(self, action, entity, variable):
        ...

    def export_closing_steps(self):
        ...

    def execute(self, cleanup=False, timeout=10):
        self.open()
        self.ng_export()
        self.close()

    ######################################################################
    
    def ng_export(self):
        self.ng_export_preamble()
        self.ng_export_metadata()
        self.ng_export_material_definitions()
        self.ng_export_boundary_definitions()
        self.ng_export_geometry()
        self.ng_export_block_labels()
        self.ng_export_solving_steps()
        self.ng_export_postprocessing()
        self.ng_export_closing_steps()
        

    def ng_export_preamble(self):
        self.comment("PREAMBLE", 1)
        self.comment("empty", 1)
        
        self.newline(2)

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

    def ng_export_geometry(self):
        self.comment("GEOMETRY", 1)
        self.comment("empty", 1)
        
        self.newline(2)

    def ng_export_solving_steps(self):
        self.comment("SOLVER")
        self.comment("empty", 1)
        
        self.newline(2)

    def ng_export_postprocessing(self):
        self.comment("POSTPROCESSING AND EXPORTING", 1)
        self.comment("empty", 1)
        
        self.newline(2)


    def ng_export_closing_steps(self):
        self.comment('CLOSING STEPS', 1)
        self.comment("empty", 1)
        
        self.newline(2)

    ###################################################################

    def render_geo(self):
        plt.figure()
        for nodes, attrs in self.H.edges.items():
            n1, n2 = nodes
            cp = n1.mean(n2)
            ul = 0.5 * n1.unit_to(n2).rotate(pi/2)
            ur = 0.5 * n1.unit_to(n2).rotate(-pi/2)

            plt.text(*(cp+ul), attrs['leftdomain'], horizontalalignment='right', verticalalignment='bottom')
            plt.text(*(cp+ur), attrs['rightdomain'], horizontalalignment='left', verticalalignment='top')
            
            plt.plot([n1.x, n2.x], [n1.y, n2.y], 'k-')
            plt.scatter(*n1, c='red', zorder=12)
            plt.scatter(*n2, c='red', zorder=12)

        # render material labels
        for name, mat_i in self.mat.items():
            for xi, yi in mat_i.assigned:
                plt.scatter(xi, yi, c='green', s=80, zorder=100)
                plt.text(xi+0.1, yi+0.1, name, verticalalignment='bottom', horizontalalignment='left')


        plt.show()

