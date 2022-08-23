class Point(list):
    def __new__(self, value, name=None, values=None):
        s = super(Point, self).__new__(self, value)
        return s

    @staticmethod
    def create(values, peptide):
        point = Point(values)
        point.peptide = peptide
        return point
