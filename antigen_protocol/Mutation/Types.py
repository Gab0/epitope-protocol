 
class MutationSummary():
    """

    Holds information about a single aminoacid position that shows multiple
    residues across different alignments.

    """
    def __init__(self, baseAA, position):
        self.baseAA = baseAA
        self.position = position
        self.variations = {}

        self.SpecialSwap = []

    def addvar(self, base, organism_id):
        if base not in self.variations.keys():
            self.variations[base] = []

        self.variations[base].append(organism_id)

    def __str__(self):
        header = "%i -> %s" % (self.position, self.baseAA)
        message = [header]
        for b in sorted(self.variations.keys()):
            message.append("\t%s ->" % b)
            for organism in sorted(self.variations[b]):
                message.append("\t\t%s" % organism)

        if self.SpecialSwap:
            message.append("Warning: Special aminoacid change!\n" +
                           " & ".join(self.SpecialSwap))

        return "\n".join(message)

    def get_nbvar(self):
        return len(self.variations.keys())

    def get_nbseq(self):
        return sum([
            len(self.variations[k])
            for k in self.variations
        ])


class OutputMutationFile():
    def __init__(self, mutations, filename_prefix, output_filepath):
        self.mutations = mutations
        self.filename_prefix = filename_prefix
        self.output_filepath = output_filepath

    def write_recipe(self, filepath):
        with open(filepath, 'w') as f:
            content = "\n".join([str(m) for m in self.mutations])
            f.write(content)


class Mutation():
    """
    Information on a single mutation.
    """
    pos: int
    fromAA: str
    toAA: str

    def __init__(self, p, f, t):
        self.pos = p
        self.fromAA = f
        self.toAA = t

    def __str__(self):
        return self.show()

    def show_spread(self) -> str:
        return f"{self.pos}:{self.fromAA} -> {self.toAA}"

    def show(self) -> str:
        return "".join([self.fromAA, str(self.pos), self.toAA])

    @classmethod
    def read(cls, text):
        text = text.strip()
        f = text[0]
        t = text[-1]
        p = int(text[1:-1])
        return cls(p, f, t)
