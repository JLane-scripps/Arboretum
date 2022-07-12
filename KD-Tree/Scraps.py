def add(self, psm: PSM) -> None:
    if len(self.tree) == 0:
        self.tree.append(psm)
        count + 1
        return

    if psm < self.tree[0]:
        self.tree.append(0)
        for i in reversed(range(count)):
            self.tree[i] = self.tree[i - 1]
        self.tree[0] = psm
        count + 1
        return

    # If new psm is the same value as current self.mid, compares rt values to place new psm above or below self.mid and adjust
    # the rest of the tree.
    if psm == self.tree[self.mid]:
        if psm.rt < self.tree[self.mid].rt:
            self.tree = self.tree[:self.mid - 1] + [psm] + self.tree[self.mid - 1:]
        if psm.rt > self.tree[self.mid].rt:
            self.tree = self.tree[:self.mid] + [psm] + self.tree[self.mid:]
        count + 1
        return
        # end

    if psm.mz > self.tree[-1]:
        self.tree.append(psm)
        count + 1
        self.mid = count / 2
        return

    # if new psm fits between first and self.mid, find closer point, iterate by 2, move everything up and insert.
    # quarters are an estimation of how many psms are likely to be between first and self.mid, not exactly how many.
    if self.tree[0].mz <= psm.mz <= self.tree[self.mid].mz:
        quarter1 = psm.mz - self.tree[0].mz  # this is not an index location, but a difference in mz value
        quarter2 = self.tree[self.mid].mz - psm.mz  # not an index location, just likely starting point
        if quarter2 > quarter1:
            for i in reversed(range(self.mid, 0, 2)):
                if self.tree[i] < psm:
                    if self.tree[i + 1] < psm:
                        self.tree = self.tree[:i] + [psm] + self.tree[i:]
                    else:
                        self.tree = self.tree[:i - 1] + [psm] + self.tree[i - 1:]
        else:  # if new psm likely to be closer to beginning
            for i in range(0, self.mid, 2):
                if self.tree[i] < psm:
                    if self.tree[i - 1] < psm:
                        self.tree = self.tree[:i] + [psm] + self.tree[i:]
                # Is there an appropriate throw condition for here? To ensure serach doesn't go out of range?
        count + 1
        self.mid = count / 2
        return
        # end of adding before self.mid

    if self.tree[self.mid].mz <= psm.mz <= self.tree[-1].mz:
        quarter3 = psm.mz - self.tree[self.mid].mz  # this is not an index location, but a difference in mz value
        quarter4 = self.tree[-1].mz - psm.mz  # not an index location, just likely starting point
        if quarter4 > quarter3:  # if new psm likely to be closer to the end of list
            for i in reversed(range(-1, self.mid, 2)):
                if self.tree[i] < psm:
                    if self.tree[i + 1] < psm:
                        self.tree = self.tree[:i] + [psm] + self.tree[i:]
                    else:
                        self.tree = self.tree[:i - 1] + [psm] + self.tree[i - 1:]
        else:  # if new psm likely to be closer to self.mid
            for i in range(self.mid, -1, 2):
                if self.tree[i] < psm:
                    if self.tree[i - 1] < psm:
                        self.tree = self.tree[:i] + [psm] + self.tree[i:]
                # Is there an appropriate throw condition for here? To ensure serach doesn't go out of range?
        count + 1
        self.mid = count / 2
        return
        # end of adding after self.mid
    # end of add