class FilePrint(object):

    def __init__(self):
        pass

    def print_array(self, array, filename):
        fout = open(filename, 'w')

        for row in range(len(array)):
            fout.write('%f\n' % (array[row]))

        fout.close()

    def print_multi_array(self, array, filename):
        fout = open(filename, 'w')

        for row in range(len(array)):
            for col in range(len(array[row])):
                fout.write('%f\t' % array[row][col])
            fout.write('\n')

        fout.close()
