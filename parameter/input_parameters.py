class Parameters(object):

    def __init__(self, args):
        """
        Assigns the input parameters for the simulations to variables
        """
        self.n_conf = args.nconf

        if args.calpha:
            self.c_alpha = True
        else:
            self.c_alpha = False

        if args.eval:
            self.eval = True
        else:
            self.eval = False

        if args.torsion:
            self.torsion = True
        else:
            self.torsion = False
