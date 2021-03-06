#!/usr/bin/env python3
import argparse
import subprocess
import logging

from pathlib import Path
import treesearch as ts 

def valid_file(path_str):
    p = Path(path_str)
    if p.is_file():
        return p
    else:
        raise FileNotFoundError(f"could not find a file at {path_str}")

def valid_output(path_str):
    p = Path(path_str)
    if p.parent.exists():
        return p
    else:
        raise NotADirectoryError(f"{p.parent} is not a valid directory")

def valid_start_strategy(strat_str):
    valid_strats = {'random', 'nj', 'upgma'}
    if strat_str in valid_strats:
        return strat_str
    else:
        raise Exception(f"Not a valid starting strategy: {strat_str}")

def is_positive(x):
    ivalue = int(x)
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"Number of processes must be at least 1 (got {x})")
    return ivalue

def make_parser():
    parser = argparse.ArgumentParser(
        description="Wrapper for tree search (serial version)"
    )
    # if i feel like going nuts i can try to validate the formats of these - no thanks
    parser.add_argument("-t", "--template", type=valid_file, help="baseml control file template path", required=True)
    parser.add_argument("-s", "--seq", type=valid_file, help="sequence path (PHYLIP format)", required=True)
    parser.add_argument("-o", "--output", type=valid_output, help="Directory to output results to", required=True)
    parser.add_argument("-z", "--seed", type=int, help="Seed PRNG")
    parser.add_argument("-S", "--start", type=valid_start_strategy, help="Starting tre strategy: 'random', 'nj', or 'upgma'", default='nj') # for now...
    parser.add_argument("-M", "--max_iter", type=int, help="Maximum number of iterations", default=2)
    parser.add_argument("-X", "--mpi_method", type=str, help="Pick parallelization strategy", default="1")
    parser.add_argument("-l", "--logging", choices=["DEBUG", "INFO", "WARNING", "ERROR","CRITICAL"], default="INFO", help="Set the logging level")
    return parser

def main():
    parser = make_parser()
    parser.add_argument("-P", "--MPI", type=is_positive, help="number of MPI processes. If set to 1, program executes serially", default=1)
    try:
        args = parser.parse_args()
    except Exception as e:
        print(e)
        parser.print_help()
        exit(1)
    logging.basicConfig(level=args.logging)
    if args.MPI > 1:
        logging.info(f"Running MPI with {args.MPI} processes")
        # kinda janky, should probably not even do that
        proc_args = [
            "mpiexec", "-n", str(args.MPI), 
            "python", "nhts-mpi.py", 
            "-t", str(args.template.absolute()),
            "-s", str(args.seq.absolute()),
            "-o", str(args.output.absolute()),
            "-M", str(args.max_iter),
            "-X", str(args.mpi_method),
            "-l", str(args.logging)
        ]
        if args.seed:
            proc_args.extend(["-z", str(args.seed)])
        if args.start:
            proc_args.extend(["-S", str(args.start)])
        subprocess.Popen(proc_args, cwd=Path(__file__).parent).wait()
        # raise NotImplementedError("MPI not yet implemented")
    else:
        logging.info(f"Running serial algorithm")
        result = ts.serial_single_shift_search(args)

if __name__ == "__main__":
    main()
