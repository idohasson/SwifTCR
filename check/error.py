class Errors:
    def __init__(self):
        self.data = {}

    def ADD_LIST(self, errors_list, cont_name):
        if cont_name not in self.data:
            self.data[cont_name] = []

        for entry in errors_list:
            if entry != [[]] and entry != []:
                self.data[cont_name].append(entry)

    def MERGE(self, new_errors):

        for cont_name in list(new_errors.data.keys()):

            if cont_name not in self.data:
                self.data[cont_name] = []

            for entry in new_errors.data[cont_name]:
                self.data[cont_name].append(entry)

    def SORT(self):

        for cont_name in list(self.data.keys()):

            if len(self.data[cont_name]) > 1:
                temp_sorted = sorted(self.data[cont_name], key=lambda inter: inter[0], reverse=False)

                self.data[cont_name] = temp_sorted

    def PRINT(self):
        for cont_name in list(self.data.keys()):
            print(cont_name)

            for entry in self.data[cont_name]:
                print(entry)
            print()

    def CLASS_TO_LIST(self):
        error_list = []
        for cont_name in list(self.data.keys()):

            for entry in self.data[cont_name]:
                error_list.append(entry)

        return error_list


# test every method
if __name__ == "__main__":
    # test ADD_LIST
    errors = Errors()
    errors.ADD_LIST([['a', 'b', 'c'], ['d', 'e', 'f']], 'cont1')
    errors.ADD_LIST([['g', 'h', 'i'], ['j', 'k', 'l']], 'cont2')
    errors.ADD_LIST([['m', 'n', 'o'], ['p', 'q', 'r']], 'cont3')
    errors.PRINT()

    # test MERGE
    errors2 = Errors()
    errors2.ADD_LIST([['a', 'b', 'c'], ['d', 'e', 'f']], 'cont1')
    errors2.ADD_LIST([['g', 'h', 'i'], ['j', 'k', 'l']], 'cont2')
    errors2.ADD_LIST([['m', 'n', 'o'], ['p', 'q', 'r']], 'cont3')
    errors.MERGE(errors2)
    errors.PRINT()

    # test SORT
    errors.SORT()
    errors.PRINT()

    # test CLASS_TO_LIST
    print(errors.CLASS_TO_LIST())

    # test PRINT
    errors.PRINT()


