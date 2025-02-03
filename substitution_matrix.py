class SubstitutionMatrix:
    def __init__(self, filename):
        with open(filename) as reader:
            line = reader.readline()
            characters = line.strip().split()
            self.values = {}
            for char in characters:
                self.values[char] = {}

            for line in reader:
                splitted_line = line.strip().split()
                char1 = splitted_line[0]
                for i in range(1, len(splitted_line)):
                    self.values[char1][characters[i - 1]] = float(splitted_line[i])


    def __getitem__(self, key):
        key1, key2 = key
        return self.values[key1][key2]


if __name__ == "__main__":
    matrix = SubstitutionMatrix("matrices/BLOSUM62.txt")
    print(matrix["M", "L"])
