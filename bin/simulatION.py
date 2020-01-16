#!/usr/bin/env python3
import argparse
import simulation.version


def simulator(parser, args):
    if args.command == 'simulate':
        from simulation import simulate
        simulate.run(parser, args)
    else:
        print("Please use simulate command!")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='simulatION')
    parser.add_argument("-v", "--version", help="Installed poretools version", action="version", version="%(prog)s " + str(simulation.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    parser_simulate = subparsers.add_parser(name='simulate',description="SimulatION simulation tool")
    parser_simulate.add_argument("-r", "--ref", required=True, help="Path to reference sequence file")
    parser_simulate.add_argument("-m", "--model", required=True, help="Path to pore model file")
    parser_simulate.add_argument("-c", "--config", required=True, help="Path to config file. Needs to be generated before using: simulation/WriteConfig.py")
    parser_simulate.add_argument("-e", "--error_rate", type=float, default=0.0, help="Error rate in decimal format (0.1 = 10%)")
    parser_simulate.add_argument("-n", "--reads", type=int, default=100, help="Number of reads that should be simulated")
    parser_simulate.add_argument("-o", "--output", type=str, required=True, help="Path to where output files should be written to.")
    parser_simulate.add_argument("--dirreads", type=int, default=4000, help="Option to change the reads within one subfolder (ONT standard is 4000)")
    parser_simulate.add_argument("--workerreads", type=int, default=50, help="Option to change the reads simulated per worker (Standard is 50)")
    parser_simulate.add_argument("--cut_off", type=float, default=1750.0, help="Cut-Off frequency form the low pass filtering")
    parser_simulate.add_argument("--band", type=float, default=40.0, help="Bandwidth frequency for the low pass filtering")
    parser_simulate.add_argument("--reverse", type=bool, default=False, help="Whether choose to get simulation signal of complementary sequence")
    parser_simulate.add_argument("--signal_repeat", type=bool, default=True, help="Whether choose to make signal repeat")
    parser_simulate.add_argument("--uniform_noise", type=bool, default=True, help="Whether choose to create uniform noise")
    parser_simulate.add_argument("--event_repeat", type=bool, default=True, help="Whether choose to make event repeat")
    parser_simulate.add_argument("--low_pass_filter", type=bool, default=True, help="Whether choose to filter by low pass filter")
    parser_simulate.add_argument("--gaussian_noise", type=bool, default=True, help="Whether choose to create gaussian noise")
    parser_simulate.add_argument("--correction", type=bool, default=True, help="Whether choose to correct signal")
    # parser_simulate.add_argument("--filter", type=str, default="None", help="Filter reference genome for specific sizes")
    parser_simulate.set_defaults(func=simulator)

    args = parser.parse_args()
    args.func(parser, args)
