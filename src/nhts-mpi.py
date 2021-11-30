#!/usr/bin/env python3
import argparse
import treesearch as ts

from pathlib import Path
from Bio import AlignIO
from mpi4py import MPI
from nhts import make_parser

def crude_partition(data, n):
    # i don't really like this, but it works fine enough
    return [data[i:i+n] for i in range(0, len(data),n)]

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    comm_size = comm.Get_size()
    
    parser = make_parser()
    args = parser.parse_args()

    with open(args.template, 'r') as fi:
        ctl_template = "".join(fi.read())
    ctl_template = ctl_template.replace("#SEQFILE", str(args.seq.absolute()))
    
    if rank == 0: # director
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
            next_trees = comm.scatter(partition, root=0)
        
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
                next_trees = comm.scatter([False] * comm_size, root=0)
                break
            else:
                prev_tree = best_tree
                unfiltered_neighbors = ts.nearest_neighbors(best_tree)
                neighbors = []
                for t in unfiltered_neighbors:
                    dont_skip = True
                    for v in visited_trees:
                        if t.compare(v, unrooted=True)['norm_rf'] == 0:
                            dont_skip = False
                            break
                    if dont_skip:
                        neighbors.append(t)
            print(f"{stepnum}:  {g_best_likelihood}")
        
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
            if not next_trees:
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

if __name__ == "__main__":
    main()
