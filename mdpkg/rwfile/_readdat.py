

def read_dat(file, dtype=float):
    data = {}
    labels = []
    with open(file, 'r') as fd:
        line = fd.readline().split()[1:]
        for label in line:
            data[label] = []
            labels.append(label)
        line = fd.readline().split()
        while line:
            for i, label in enumerate(labels):
                data[label].append( dtype(line[i]) )
            line = fd.readline().split()
    return data
