#!/usr/bin/env python3
import argparse
import treesearch as ts
import logging

from pathlib import Path
from Bio import AlignIO
from mpi4py import MPI
from nhts import make_parser
from numpy import array_split
def crude_partition(data, n):
    # i don't really like this, but it works fine enough
    # return [data[i:i+n] for i in range(0, len(data),n)]
    return array_split(data, n)

def main():
    parser = make_parser()
    args = parser.parse_args()
    logging.basicConfig(
        level=args.logging,
        format="%(asctime)s %(levelname)-4s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    strategy_map = {
        '1' : strategy_1,
        '2' : strategy_2
    }
    strategy_map[args.mpi_method](args)

def strategy_1(args):
    """
    This strategy just distributes the work of each neighborhood 
    across all the processes. Rank 0 is in charge of doing that.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    comm_size = comm.Get_size()
    if rank == 0:
        logging.info(f"Using strategy 1 (arg: {args.mpi_method})")
    
    with open(args.template, 'r') as fi:
        ctl_template = "".join(fi.read())
    ctl_template = ctl_template.replace("#SEQFILE", str(args.seq.absolute()))
    
    if rank == 0: # director
        ndigits = len(str(args.max_iter))
        logging.debug(f"{rank}: Choosing starting tree and constructing NNI neighborhood")
        # Get initial trees
        alignment = AlignIO.read(args.seq, format="phylip-relaxed")
        if args.start == 'nj':
            visited_trees = [ts.nj_tree(alignment)]
        elif args.start == 'random':
            visited_trees = [ts.random_tree(alignment)]
        else:
            raise NotImplementedError(f"{args.starting_tree} strategy not implemented")
        del alignment

        stepnum = 0
        best_tree = visited_trees[0]
        prev_tree = best_tree

        neighbors = ts.nearest_neighbors(best_tree)
        visited_trees += neighbors
        neigbors = [best_tree] + neighbors

        l_best_likelihood = float('-inf')
        l_best_info = None
        g_best_likelihood = float('-inf')
        g_best_info = {"log likelihood" : float('-inf')}

        stop_working = False
        
        for iter_num in range(args.max_iter):
            partition = crude_partition(neighbors, comm_size)
            partition = [partition[-1]] + partition[:-1]
            logging.debug(f"{rank}: Sending neighborhood to processes")
            next_trees = comm.scatter(partition, root=0) 
            logging.debug(f"{rank}: Finished sending, Received neighborhood of size {len(next_trees)}")
            for tree in next_trees:
                opath = Path(f"{args.output}/r{rank}_step_{stepnum}")
                result = ts.run_single_shift_baseml(tree, opath, ctl_template)
                if result["log likelihood"] > l_best_likelihood:
                    l_best_likelihood = result["log likelihood"]
                    l_best_info = result
                    l_best_tree = result["tree"]
                stepnum += 1
            # gather the best trees
            local_results = comm.gather(l_best_info, root=0)
            g_best_info = max([g_best_info] + local_results, key=lambda x : x["log likelihood"])
            g_best_likelihood = g_best_info["log likelihood"]
            best_tree = g_best_info["tree"]
            if best_tree == prev_tree:
                # send out an empty list and then stop working
                print(f"{str(iter_num).zfill(ndigits)}:  Did not find a better tree, stopping...")
                next_trees = comm.scatter([None] * comm_size, root=0)
                break
            else:
                prev_tree = best_tree
                unfiltered_neighbors = ts.nearest_neighbors(best_tree)
                neighbors = []
                for t in unfiltered_neighbors:
                    dont_skip = True
                    # this is terrible, really
                    for v in visited_trees:
                        if t.compare(v, unrooted=True)['norm_rf'] == 0:
                            dont_skip = False
                            break
                    if dont_skip:
                        neighbors.append(t)
            print(f"{str(iter_num).zfill(ndigits)}:  {g_best_likelihood}", flush=True)
        
        print(f"Best tree topology:", g_best_info["tree"].write(format=9))
        print("see results in ", g_best_info["Path"])
    else: # worker
        # rank, step num in output name
        # wait for my tree
        stepnum = 0
        
        l_best_likelihood = float('-inf')
        l_best_info = None
        
        for iter_num in range(args.max_iter):
            next_trees = comm.scatter(None, root=0)
            logging.debug(f"{rank}: Received neighborhood of size {len(next_trees)}")
            if next_trees is None:
                # kind of crude - what if one process doesn't have neighboring trees for one iteration ? 
                # it just stops doing work forever?
                break
            for tree in next_trees:
                opath = Path(f"{args.output}/r{rank}_step_{stepnum}")
                result = ts.run_single_shift_baseml(tree, opath, ctl_template)
                if result["log likelihood"] > l_best_likelihood:
                    l_best_likelihood = result["log likelihood"]
                    l_best_info = result
                    l_best_tree = result["tree"]
                stepnum += 1
            comm.gather(l_best_info, root=0)

def strategy_2(args):
    """
    For this strategy, each thread does the same thing
    allgather visited trees
    allreduce max of likelihood 

    visited trees kept track in a set
    trees are encoded by a bipartition of leaves

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    comm_size = comm.Get_size()
    if rank == 0:
        logging.info(f"Using strategy 2 (arg: {args.mpi_method})")
        logging.error("Not yet implemented")


if __name__ == "__main__":
    main()
