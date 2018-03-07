class Parameters(object):

    def __init__(self, args):
        self.n_conf = args.nconf

        if args.calpha:
            self.c_alpha = True
        else:
            self.c_alpha = False

        if args.eval:
            self.eval = True
        else:
            self.eval = False
