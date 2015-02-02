from scripts.periodic_table import Name as s2n
from scripts.periodic_table import AtomicNum
from scripts.gamess_basis_rename import invert_dict


def convert_basis_int(basis):
    """
    Return the given basis (in GAMESS-US format from the EMSL database)
    in the ORCA format suitable for adding to an input file.
    """
    n2s = invert_dict(s2n)
    output = ''
    # break up the file
    tmp = [line.split() for line in basis.splitlines()]
    basis_name = tmp[0][1]
    output += '# ' + basis_name + '\n'
    # discard unwanted lines, like comments
    tmp = [line for line in tmp if line != []]
    tmp = [line for line in tmp if line[0] != '!']
    tmp = [line for line in tmp if line[0][0] != '$']

    for index, line in enumerate(tmp):
        # There are 3 types of lines we must handle:
        # 1. Those that contain an element name (HYDROGEN, CARBON, COPPER, etc.) [length == 1]
        # 2. Those that contain shell info (S 3, L 1, etc.) [length == 2]
        # 3. Those that contain primitives (3 columns, first is an integer) [length >= 3]
        if len(line) == 1:
            if index > 0:
                output += ' end\n'
            # Need to convert 'CARBON' -> 6, etc.
            output += '{} {}\n'.format('newgto', AtomicNum[n2s[line[0].lower()]])
        if len(line) == 2:
            output += ' {} {}\n'.format(*line)
        if len(line) == 3:
            output += '  {} {} {}\n'.format(*line).replace('D', 'E')
        if len(line) == 4:
            output += '  {} {} {} {}\n'.format(*line).replace('D', 'E')
    output += ' end\n'
    return output

def convert_basis_ext(basis):
    """
    Return the given basis (in GAMESS-US format from the EMSL database)
    in the ORCA format suitable for reading from an external file.
    """
    output = ''
    # break up the file
    tmp = [line.split() for line in basis.splitlines()]
    basis_name = tmp[0][1]
    output += '# ' + basis_name + '\n'
    # discard unwanted lines, like comments
    tmp = [line for line in tmp if line != []]
    tmp = [line for line in tmp if line[0] != '!']
    tmp = [line for line in tmp if line[0][0] != '$']

    for line in tmp:
        # There are 3 types of lines we must handle:
        # 1. Those that contain an element name (HYDROGEN, CARBON, COPPER, etc.) [length == 1]
        # 2. Those that contain shell info (S 3, L 1, etc.) [length == 2]
        # 3. Those that contain primitives (3 columns, first is an integer) [length >= 3]
        if len(line) == 1:
            output += '{}\n'.format(*line)
        if len(line) == 2:
            output += ' {} {}\n'.format(*line)
        if len(line) == 3:
            output += '  {} {} {}\n'.format(*line).replace('D', 'E')
        if len(line) == 4:
            output += '  {} {} {} {}\n'.format(*line).replace('D', 'E')
    output += 'STOP\n'
    return output
