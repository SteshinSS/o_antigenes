from collections import namedtuple

class Monosaccharide:
    def __init__(self, name, source_atom, dest_atom, next_monosacc, text_pos=-1, edge_type = ""):
        # e.g. "aLRha(1->6)"
        self.name = name  # aLRha
        self.source_atom = source_atom  # 1
        self.dest_atom = dest_atom  # 6
        self.next_monosacc = next_monosacc  # id of an adjacent monosacc in the cycle

        # position in line
        # need only for find adjacent monosacc in the cycle for noncycle monosaccharides
        self.text_pos = text_pos

        # edge_type:
        # "" if (1->6)
        # P if (1-P-6)
        # N if (1-N-6)
        self.edge_type = edge_type


class Antigene:
    def __init__(self, lines, name):
        self.name = name
        self.cycle = []
        self.additional = []

        if lines is None:
            return

        cycle_index = self.__find_cycle__(lines)
        self.__parse_cycle__(lines[cycle_index])

        if len(lines) == 3:
            if cycle_index == 0:
                self.__parse_additional__(line=lines[2], edges=lines[1])
            if cycle_index == 2:
                self.__parse_additional__(line=lines[0], edges=lines[1])

        if len(lines) == 5:
            self.__parse_additional__(line=lines[0], edges=lines[1])
            self.__parse_additional__(line=lines[4], edges=lines[3])

    # return index of line with cycle
    # cycle begins with "->"
    def __find_cycle__(self, lines):
        answer = -1
        for i in range(len(lines)):
            for char in lines[i]:
                if char == ' ':
                    continue

                if char == '-':
                    if answer == -1:
                        answer = i
                    else:
                        raise NameError('There are two lines with cycle ' + self.name)
                break
        if answer == -1:
            raise NameError('There is antigene without cycle ' + self.name)
        return answer

    def __parse_cycle__(self, line):
        parts = line.split(")")

        # find first digit occurence in cycle
        # it looks like "->4)" and need for last monosaccharide
        last_destination = parts[0][-1]

        current_pos = 0
        for i, part in enumerate(parts):
            if i != 0:
                new_monosacc = part.split("(")  # ['aLRha', '1->5']
                name = new_monosacc[0]
                source_atom = new_monosacc[1][0]
                dest_atom = -1
                next_monosacc = i
                if len(new_monosacc[1]) == 4:
                    dest_atom = new_monosacc[1][3]
                else:
                    if i + 1 == len(parts):
                        # it's last monosacc in cycle
                        dest_atom = last_destination
                        next_monosacc = 0
                    else:
                        raise NameError('Unusual monosacc ' + part + ' in ' + self.name)

                if not (source_atom.isdigit() and dest_atom.isdigit()):
                    raise NameError('Unusual monosacc ' + part + ' in ' + self.name)
                self.cycle.append(Monosaccharide(name, source_atom, dest_atom, next_monosacc, current_pos))

            current_pos += len(part) + 1  # 1 because of cutted out ")"

    def __parse_additional__(self, line, edges):
        parts = line.split("->")

        # position of last seen "|" in edge line
        last_delimiter = -1
        for i, part in enumerate(parts):
            # cut out edge from previous monosacc "3)"
            if i != 0:
                part = part[2:]

            # erase leading spaces
            first_non_space = part.rfind(" ")
            if first_non_space != -1:
                part = part[first_non_space + 1:]

            if part == "":
                break

            new_monosacc = part.split("(")  # ['aLRha', '1']
            name = new_monosacc[0]
            source_atom = new_monosacc[1]
            dest_atom = parts[i + 1][0]
            if not (source_atom.isdigit() and dest_atom.isdigit()):
                raise NameError('Unusual monosacc ' + part + ' in ' + self.name)

            edge_position = edges.find("|", last_delimiter + 1)
            last_delimiter = edge_position
            monosacc_position = -1
            for i in range(len(self.cycle)):
                if (edge_position >= self.cycle[i].text_pos and
                        ( (i + 1 == len(self.cycle)) or (edge_position < self.cycle[i + 1].text_pos))):
                    monosacc_position = i
                    break

            self.additional.append(Monosaccharide(name, source_atom, dest_atom, monosacc_position))

    def show(self):
        print(self.name)
        monosaccharides = self.cycle
        monosaccharides.extend(self.additional)

        for monosacc in monosaccharides:
            if monosacc.edge_type == "":
                print('{}({}->{}) to the {}'.format(monosacc.name,
                                                    monosacc.source_atom,
                                                    monosacc.dest_atom,
                                                    monosacc.next_monosacc))
            else:
                print('{}({}-{}-{}) to the {}'.format(monosacc.name,
                                                      monosacc.source_atom,
                                                      monosacc.edge_type,
                                                      monosacc.dest_atom,
                                                      monosacc.next_monosacc))

def ReadRawAntigenes(path):
    data = open(path)
    RawAntigene = namedtuple("RawAntigen", "name lines")
    RawAntigenes = []

    last_name = ""
    antigene_lines = []
    for i, line in enumerate(data):
        # delete '\n' symbol
        line = line[:-1]
        # if we found new antigene name
        if line[0] == 'O':
            # if we already read something
            if antigene_lines:
                RawAntigenes.append(RawAntigene(last_name, antigene_lines))
            last_name = line
            antigene_lines = []
        else:
            antigene_lines.append(line)
    return RawAntigenes

def GenerateAntigenes():
    raw_antigenes = ReadRawAntigenes('database.txt')
    antigenes = []
    for raw_antigene in raw_antigenes:

        if raw_antigene.name == "O28_ab":
            O28_ab = Antigene(None, name="O28_ab")
            O28_ab.cycle.append(Monosaccharide("bDGlc","1","3","1"))
            O28_ab.cycle.append(Monosaccharide("DGro","1","4","2",edge_type="P"))
            O28_ab.cycle.append(Monosaccharide("bDGlcNAc","1","3","3"))
            O28_ab.cycle.append(Monosaccharide("aDGlcNAc","1","3","0"))
            antigenes.append(O28_ab)
            continue

        if raw_antigene.name == "O28_ac":
            O28_ac = Antigene(None, name="O28_ac")
            O28_ac.cycle.append(Monosaccharide("DGro", "1", "4", "1", edge_type="P"))
            O28_ac.cycle.append(Monosaccharide("bDGlcNAc", "1", "3", "2"))
            O28_ac.cycle.append(Monosaccharide("aDGlcNAc", "1", "3", "0"))
            antigenes.append(O28_ac)
            continue

        if raw_antigene.name == "O29":
            O29 = Antigene(None, name="O29")
            O29.cycle.append(Monosaccharide("DGro", "1", "6", "1", edge_type="P"))
            O29.cycle.append(Monosaccharide("bDGlc", "1","4","2"))
            O29.cycle.append(Monosaccharide("aLFucNAc", "1","3","3"))
            O29.cycle.append(Monosaccharide("bDGlcNAc", "1","3", "0"))
            O29.additional.append(Monosaccharide("aDGlc", "1","6","5"))
            O29.additional.append(Monosaccharide("aDGal","1","3","2"))
            antigenes.append(O29)
            continue

        if raw_antigene.name == "O37":
            O37 = Antigene(None, name="O37")
            O37.__parse_cycle__(raw_antigene.lines[0])
            O37.additional.append(Monosaccharide("DGro", "1","3","0",edge_type="P"))
            antigenes.append(O37)
            continue

        if raw_antigene.name == "O42":
            O42 = Antigene(None, name="O42")
            O42.cycle.append(Monosaccharide("DGro", "1","4","1",edge_type="P"))
            O42.cycle.append(Monosaccharide("bDGlcNAc", "1","3","2"))
            O42.cycle.append(Monosaccharide("bDGalf2Ac","1","3","3"))
            O42.cycle.append(Monosaccharide("aDGlcNAc","1","2","0"))
            O42.additional.append(Monosaccharide("aDGlc","1","3","1"))
            antigenes.append(O42)
            continue

        if raw_antigene.name == "O82":
            O82 = Antigene(None, name="O82")
            O82.__parse_cycle__(raw_antigene.lines[2])
            O82.additional.append(Monosaccharide("DGroA","2","6","0",edge_type="P"))
            antigenes.append(O82)
            continue

        if raw_antigene.name == "O100":
            O100 = Antigene(None, name="O100")
            O100.__parse_cycle__(raw_antigene.lines[0])
            O100.additional.append(Monosaccharide("DGro","1","6","0",edge_type="P"))
            antigenes.append(O100)
            continue

        if raw_antigene.name == "O112_ac":
            O112_ac = Antigene(None, name="O112_ac")
            O112_ac.__parse_cycle__(raw_antigene.lines[2])
            O112_ac.additional.append(Monosaccharide("bDGlcNAc4,6(S)Pyr","1","3","0"))
            antigenes.append(O112_ac)
            continue

        if raw_antigene.name == "O118":
            O118 = Antigene(None, name="O118")
            O118.cycle.append(Monosaccharide("DRibitol","5","6","1",edge_type="P"))
            O118.cycle.append(Monosaccharide("aDGal","1","3","2"))
            O118.cycle.append(Monosaccharide("aLFucNAm","1","3","3"))
            O118.cycle.append(Monosaccharide("bDGlcNAc","1","3","0"))
            antigenes.append(O118)
            continue

        if raw_antigene.name == "O130":
            O130 = Antigene(None, name="O130")
            O130.__parse_cycle__(raw_antigene.lines[2])
            O130.additional.append(Monosaccharide("Gro","2","4","4",edge_type="P"))
            O130.additional.append(Monosaccharide("bDGalNAc","1","3","0"))
            antigenes.append(O130)
            continue

        if raw_antigene.name == "O143":
            O143 = Antigene(None, name="O143")
            O143.__parse_cycle__(raw_antigene.lines[2])
            O143.additional.append((Monosaccharide("Gro","2","6","0",edge_type="N")))
            antigenes.append(O143)
            continue

        if raw_antigene.name == "O149":
            O149 = Antigene(None, name="O149")
            O149.cycle.append(Monosaccharide("bDGlcNAc4,6(S)Pyr","1","3","1"))
            O149.cycle.append(Monosaccharide("bLRha","1","4","2"))
            O149.cycle.append(Monosaccharide("bDGlcNAc","1","3","0"))
            antigenes.append(O149)
            continue

        if raw_antigene.name == "O151":
            O151 = Antigene(None, name="O151")
            O151.cycle.append(Monosaccharide("DRibitol","5","6","1",edge_type="P"))
            O151.cycle.append(Monosaccharide("aDGal","1","3","2"))
            O151.cycle.append(Monosaccharide("aLFucNAm","1","3","3"))
            O151.cycle.append(Monosaccharide("bDGlcNAc","1","2","0"))
            O151.additional.append(Monosaccharide("bDGlcNAc","1","4","3"))
            antigenes.append(O151)
            continue

        if raw_antigene.name == "O152":
            a = Antigene(None, name="O152")
            a.cycle.append(Monosaccharide("aDGlcNAc","1","6","1",edge_type="P"))
            a.cycle.append(Monosaccharide("aDGlc","1","2","2"))
            a.cycle.append(Monosaccharide("bDGlc","1","3","3"))
            a.cycle.append(Monosaccharide("bDGlcNAc","1","3","0"))
            a.additional.append(Monosaccharide("bLRha","1","4","0"))
            antigenes.append(a)
            continue

        if raw_antigene.name == "O156":
            b = Antigene(None, name="O156")
            b.__parse_cycle__(raw_antigene.lines[0])
            b.additional.append(Monosaccharide("aDGal4,6(R)Pyr","1","3","0"))
            antigenes.append(b)
            continue

        if raw_antigene.name == "O160":
            c = Antigene(None, name="O160")
            c.cycle.append(Monosaccharide("bDGlcNAc","1","3","1"))
            c.cycle.append(Monosaccharide("aDGal","1","6","2",edge_type="P"))
            c.cycle.append(Monosaccharide("bDGal","1","3","3"))
            c.cycle.append(Monosaccharide("bDGalNAc","1","4","0"))
            c.additional.append(Monosaccharide("bDGlc","1","6","0"))
            antigenes.append(c)
            continue

        if raw_antigene.name == "O172":
            d = Antigene(None, name="O172")
            d.cycle.append(Monosaccharide("aLFucNAc", "1","4","1"))
            d.cycle.append(Monosaccharide("aDGlc6Ac","1","4","2",edge_type="P"))
            d.cycle.append(Monosaccharide("aDGlc","1","3","3"))
            d.cycle.append(Monosaccharide("aLFucNAc","1","3","4"))
            d.cycle.append(Monosaccharide("aDGlcNAc","1","3","0"))
            antigenes.append(d)
            continue

        if raw_antigene.name == "O173":
            e = Antigene(None, name="O173")
            e.cycle.append(Monosaccharide("aDGlc", "1","6","1",edge_type="P"))
            e.cycle.append(Monosaccharide("aDGlc","1","2","2"))
            e.cycle.append(Monosaccharide("bDGlc","1","3","3"))
            e.cycle.append(Monosaccharide("bDGlcNAc","1","3","0"))
            e.additional.append((Monosaccharide("aLFuc","1","4","0")))
            antigenes.append(e)
            continue

        if raw_antigene.name == "O181":
            f = Antigene(None, name="O181")
            f.cycle.append(Monosaccharide("aDGalNAc6Ac","1","6","1"))
            f.cycle.append(Monosaccharide("aDGlc","1","4","2",edge_type="P"))
            f.cycle.append(Monosaccharide("aLQuiNAc","1","3","3"))
            f.cycle.append(Monosaccharide("bDGlcNAc","1","4","0"))
            f.additional.append(Monosaccharide("aLQuiNAc","1","3","0"))
            antigenes.append(f)
            continue

        antigenes.append(Antigene(raw_antigene.lines, raw_antigene.name))
    return antigenes

if __name__=='__main__':
    antigenes = GenerateAntigenes()
    for gene in antigenes:
        gene.show()
        print()


