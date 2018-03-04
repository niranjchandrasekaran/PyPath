class Parameters(object):

    def __init__(self, args):
        self.n_conf = args.nconf

        if args.calpha:
            self.c_alpha = True
        else:
            self.c_alpha = False

        try:
            self.mode = args.mode
        except AttributeError:
            self.mode = None

        try:
            with open(args.cons, 'r') as fopen:
                self.cons = [line.rstrip() for line in fopen]
        except TypeError:
            self.cons = None

        try:
            self.exaggeration = args.exag
        except AttributeError:
            self.exaggeration = None