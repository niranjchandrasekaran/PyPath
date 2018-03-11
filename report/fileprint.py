class FilePrint(object):

    def __init__(self):
        """
        Methods for printing arrays to file
        """
        pass

    def print_array(self, array, filename):
        """
        Prints a single dimensional array to file

        :param array: 1D array
        :param filename: output file name
        """
        fout = open(filename, 'w')

        for row in range(len(array)):
            fout.write('%f\n' % (array[row]))

        fout.close()

    def print_multi_array(self, array, filename):
        """
        Prints multidimesional array to file

        :param array: Multidimensional array
        :param filename: output file name
        """
        fout = open(filename, 'w')

        for row in range(len(array)):
            for col in range(len(array[row])):
                fout.write('%f\t' % array[row][col])
            fout.write('\n')

        fout.close()
