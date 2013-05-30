class OrcaParser:
    """
    A parser for ORCA input and output files.

    :param file_name: the ORCA input file name
    :type file_name: str

    >>> molecule = OrcaParser('caffeine.inp')
    >>> print molecule.charge
    >>> molecule.scf['guess'] = 'pmodel'
    >>> molecule.inputblock += 'blyp def2-svp def2-svp/c'
    """
    def __init__(self, file_name):
        pass

    def __str__(self):

    def reset(self):
        pass

    def load(self, file_name):
        pass

    def _parse_molecule(self):
        pass

    def _parse_inputblock(self):
        pass

    def _parse_fields(self):
        pass

    def _parse_charges(self):
        pass

    def _extract(self, keyword):
        pass

    def format_inputblock(self):
        pass

    def format_molecule(self):
        pass

    def format_multipole_fields(self):
        pass

    def to_string(self):
        pass

    def save(self, output_name):
        handle = open(output_name, 'w')
        handle.write(self.to_string())
        handle.close()

    def interatomic_distance(self, atoms_index):
        """
        returns distance between two atoms in file

        :param atoms_index: contains index of positions of atoms of interest, 0 is first atom
        ;type atoms_index: list
        """

        pass

class OrcaInputParser(OrcaParser):
    """
    A parser for ORCA input files.
    """
    def __init__(self, file_name):
        OrcaParser.__init__(self, file_name)

class OrcaOutputParser(OrcaParser):
    """
    A parser for ORCA output files.
    """
    def __init__(self, file_name):
        OrcaParser.__init(self, file_name)

    def get_energy(self, string_to_search = None):
        pass

    def _set_output_keywords(self):
        pass

    def check_method_type(self):
        pass

    
