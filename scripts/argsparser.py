#!/usr/bin/env python3
# Lanxin Zhang

import argparse
import sys


def number(ts):
    try:
        ts = float(ts)
    except ValueError:
        try:
            ts = eval(ts)
        except NameError:
            msg = "{} is not a correct number".format(ts)
            raise argparse.ArgumentTypeError(msg)
    return ts



def parse_arguments():
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--mdose', type=float, default=50,
                        help='mass of DTG dose [mg] (default: 50)')
    parent_parser.add_argument('--adh', type=float, default=1,
                        help='adherence of drug taken (default: 1)')
    parent_parser.add_argument('--cdose', type=int, default=3,
                        help='count of doses in a regimen (default: 3)')
    parent_parser.add_argument('--inputfile', type=str, default='PK_default.csv',
                        help='specify the inputfile that contains the \
                        pharmacokinetic parameters. For deterministic methods \
                        please input only one set of pk parameters \
                        (default: PK_default.csv)')
    parent_parser.add_argument('--outputfile', type=str,
                        default='{}.csv'.format(sys.argv[1]),
                        help=' specify the name of output csv file \
                        (default: {}.csv)'.format(sys.argv[1]))
    parent_parser.add_argument('--ifpdf', type=bool, default=False,
                        help='if output the graphic of results in pdf \
                        (default: False)')

    parent_parser_deterministic = argparse.ArgumentParser(add_help=False)
    parent_parser_deterministic.add_argument('--ts', type=float, default=0,
                        help='start time point to run the computation, \
                        relative to timing of first dose [hr] (default: 0)')
    parent_parser_deterministic.add_argument('--te', type=float, default=240,
                        help='end time point of the computation, \
                        relative to timing of first dose [hr] (default: 240)')
    parent_parser_deterministic.add_argument('--tspanres', type=float, default=1,
                        help='timespan in the output data [hr] \
                        should not be shorter than timesteps (default: 1)')
    parent_parser_deterministic.add_argument('--ifll', type=bool, default=False,
                        help='if long-lived and latently infected cells are \
                        considered (default: False)')

    parent_parser_deterministic12 = argparse.ArgumentParser(add_help=False)
    parent_parser_deterministic12.add_argument('--timesteps', type=number,
                        default=[1, 1, 1],
                        nargs=3, help='time steps for V, T1, T2,  [min] \
                        Fraction like 1/6 is allowed. For V and T1 the time \
                        steps should be same and cannot be shorter than T2 \
                        (default: 1 1 1)')

    parser = argparse.ArgumentParser(
        description='Compute prophylactic efficacy for DTG',
        usage='{} <command> [<options>]'.format(sys.argv[0]),
        epilog='Run \'{} COMMAND --help\' for more information\
                on a command'.format(sys.argv[0]))
    subparsers = parser.add_subparsers(
                    title='command', dest='command',
                    description='Choose one method for computing')
    parser_extr = subparsers.add_parser(
                    'extrande',
                    help='run method extrande',
                    description='hybrid stochastic simulation method',
                    usage='{} extrande [<options>]'.format(sys.argv[0]),
                    formatter_class=argparse.MetavarTypeHelpFormatter,
                    parents=[parent_parser])
    parser_ntm = subparsers.add_parser(
                    'ntm',
                    help='run next transition method',
                    description='next transition method',
                    usage='{} ntm [<options>]'.format(sys.argv[0]),
                    formatter_class=argparse.MetavarTypeHelpFormatter,
                    parents=[parent_parser, parent_parser_deterministic, parent_parser_deterministic12])

    parser_ctsm = subparsers.add_parser(
                    'ctsm',
                    help='run constant time step method',
                    description='constant time step method',
                    usage='{} ctsm [<options>]'.format(sys.argv[0]),
                    formatter_class=argparse.MetavarTypeHelpFormatter,
                    parents=[parent_parser, parent_parser_deterministic, parent_parser_deterministic12])
    parser_pgs = subparsers.add_parser(
                    'pgs',
                    help='run probability generating system',
                    description='probability generating system',
                    usage='{} pgs [<options>]'.format(sys.argv[0]),
                    formatter_class=argparse.MetavarTypeHelpFormatter,
                    parents=[parent_parser, parent_parser_deterministic])

    parser_extr.add_argument('--csimul', type=int, default=5000,
                        help='count of simulations (default: 5000)')
    parser_extr.add_argument('--tps', type=list, default=[1, 6, 18, 23],
                        help='time points to begin simulation, i.e. \
                        time of viral exposure relative to timing \
                        of first dose [hr] (default: [1,6,18,23])')
    parser_extr.add_argument('--vload', type=str, default='1',
                        help='int or str {"homo", "hetero"}, initial viral load,\
                        can be fixed number or generated randomly according to\
                        the transition mode (default: 1)')
    parser_ntm.add_argument('--rate', type=float, default=0.999,
                        help='percentage of the observed next transition \
                        (default: 0.999)')

    args = parser.parse_args()
    return args
