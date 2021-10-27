from copy import copy

from adze_modeler.metadata import Metadata
from adze_modeler.platforms.platform import Platform


class FenicsElectrostatics(Platform):
    def __init__(self, m: Metadata):
        super().__init__(m)

    def __copy__(self):
        return FenicsElectrostatics(copy(self.metadata))

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"# {str_}")
        self.newline(nb_newline)

    def export_preamble(self):
        self.write("import agros2d as a2d")
