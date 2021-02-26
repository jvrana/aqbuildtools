class CellValue(str):
    """String class that maintains 2-dimensional index.

    Primarily used for relaying parsing errors.
    """

    def add_rc(self, r, c):
        self.row = r
        self.col = c

    def rc(self):
        return self.row, self.col

    def to_cell_values(values):
        new_values = []
        for r, row in enumerate(values):
            new_values.append([])
            for c, val in enumerate(row):
                cell_val = CellValue(val)
                cell_val.add_rc(r, c)
                new_values[-1].append(cell_val)
        return new_values

    def _apply(self, fxn_name: str, args=None, kwargs=None):
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}
        x = CellValue(getattr(super(), fxn_name)(*args, **kwargs))
        x.add_rc(*self.rc())
        return x

    def lower(self):
        return self._apply("lower")

    def upper(self):
        return self._apply("upper")

    def strip(self):
        return self._apply("strip")
