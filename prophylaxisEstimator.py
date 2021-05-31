#!/usr/bin/env python3

import sys


from scripts.argsparser import parse_arguments


def checkfile(filename):
    if os.path.isfile(filename):
        ans = input('{} already exist. Continue to '
                    'overwrite? [Y/N] '.format(filename))
        if ans in 'Nn':
            sys.stderr.write('Program breaks\n')
            exit(1)

if __name__ == '__main__':
    args = parse_arguments()
    from scripts.dataprocess import *
    method = args.command
    print("------------")
    if method == 'extrande':
        from scripts.extrande import run_extrande
        checkfile(args.outputfile)
        data = run_extrande(args.mdose, args.cdose, args.adh, args.csimul,
                            args.tps, args.vload, args.inputfile)
        df = extrande_out(data, args.outputfile)
        if args.ifpdf:
            filename = args.outputfile[:-4] + '.pdf'
            extrande_plot(df, filename, args.mdose)
    else:
        if args.te <= args.ts:
            sys.stderr.write('Usage: end time should after start time\n')
            exit(1)

        if method == 'ntm':
            from scripts.ntm import run_ntm
            checkfile(args.outputfile)
            p1, p2, p3 = run_ntm(args.mdose, args.cdose, args.adh,
                                 args.inputfile, args.ts, args.te,
                                 args.timesteps, args.rate, args.ifll)
            df = deterministics_output(p1, p2, p3, args.outputfile,
                                       args.tspanres, args.timesteps,
                                       args.ts, args.te)

        elif method == 'ctsm':
            from scripts.ctsm import run_ctsm
            checkfile(args.outputfile)
            p1, p2, p3 = run_ctsm(args.mdose, args.cdose, args.adh,
                                  args.inputfile, args.ts, args.te,
                                  args.timesteps, args.ifll)
            df = deterministics_output(p1, p2, p3, args.outputfile,
                                       args.tspanres, args.timesteps,
                                       args.ts, args.te)

        elif method == 'pgs':
            from scripts.pgs import run_pgs
            checkfile(args.outputfile)
            p1, p2, p3 = run_pgs(args.mdose, args.cdose, args.adh,
                                 args.inputfile, args.ts, args.te,
                                 args.tspanres, args.ifll)
            df = deterministics_output(p1, p2, p3, args.outputfile,
                                       args.tspanres, [args.tspanres*60],
                                       args.ts, args.te)

        else:
            sys.stderr.write('No such command \n')
            exit(1)

        if args.ifpdf:
            filename = args.outputfile[:-4] + '.pdf'
            deterministics_plot(df, filename, args.mdose)

    print('Done')
