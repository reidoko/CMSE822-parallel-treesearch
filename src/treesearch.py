from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import NNITreeSearcher, DistanceCalculator, DistanceTreeConstructor
from io import StringIO
import numpy as np

import random
import ete3
import subprocess

from utils import baseml

def fill_baseml_template(ctl_template, replacements):
    result = ctl_template
    for this, with_this in replacements:
        result = result.replace(this, with_this)
    return result

def nearest_neighbors(ete3tree, is_unrooted=True):
    #ete3tree = ete3.Tree(tree_string)
    if is_unrooted:
        ete3tree.resolve_polytomy()
    bio_tree = Phylo.read(StringIO(ete3tree.write(format=9)),format='newick')
    bio_neighbors = NNITreeSearcher._get_neighbors(None, bio_tree)
    all_neighbors = [ete3.Tree(t.format('newick')) for t in bio_neighbors]
    if is_unrooted:
        for t in all_neighbors:
            t.unroot()
    return all_neighbors

def random_tree(alignment, unroot=True):
    tree = ete3.Tree()
    taxa = [x.name for x in alignment]
    tree.populate(len(taxa), names_library=taxa)
    if unroot:
        tree.unroot()
    return tree

def nj_tree(alignment, unrooted=True):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    bio_tree = constructor.build_tree(alignment)
    result = ete3.Tree(bio_tree.format('newick').replace('Inner', ''))
    if unrooted:
        result.unroot()
    return result

def best_result(results):
    best_result = {'log likelihood' : float('-inf')}
    best_log_likelihood = float('-inf')
    for result in results:
        if result:
            if result['log likelihood'] > best_log_likelihood:
                best_log_likelihood = result['log likelihood']
                best_result = result
    return best_result

def get_shift_partitions(tree, k=1):
    # So far I've only implemented this for bipartitions (k=1 shifts)
    # I believe all the ways you can choose k nodes tells you how
    assert k==1
    sps = set()
    nodeset = set(tree.traverse())
    for node in tree.traverse():
        if node.is_root():
            continue
        p1 = set(node.get_descendants())
        p1.add(node)
        p1 = frozenset(p1)
        p2 = frozenset(nodeset - p1)
        sps.add(frozenset((p1,p2)))
    return sps

def bipartition_key(node, all_leaves):
    leaves = {x.name for x in node.get_leaves()}
    notleaves = all_leaves - leaves
    bp = frozenset({frozenset(leaves), frozenset(notleaves)})
    return bp

# kinda cursed but i wanted to be able to use a tree as a key
def bipartition_representation(tree):
    all_leaves = {x.name for x in tree.get_leaves()}
    return frozenset({bipartition_key(x, all_leaves) for x in  tree.traverse()})

def single_shift_assignments(input_tree):
    tree = input_tree.copy()
    shift_partitions = get_shift_partitions(tree)
    for node in tree.traverse():
        node.original_name = node.name

    model_assignments = []
    for sp in shift_partitions:
        for modelnum, nodes in enumerate(sp, start=1):
            for node in nodes:
                if node.is_root():
                    root_name = f"#{modelnum}"
                node.name = f"{node.original_name}#{modelnum}"
        model_assignments.append(tree.write(format=8).replace(";",f"{root_name};"))
    return model_assignments

def run_baseml(tree, output_path, ctl_template, writetree = lambda x : x.write(format=9)):
    control_path = Path(f"{output_path}/baseml.ctl")
    result_path = Path(f"{output_path}/RESULT")
    tree_path = Path(f"{output_path}/tree")
    log_path = Path(f"{output_path}/LOG")
    replacements = [
        ("#TREEFILE",  tree_path.name),
        ("#OUTPUTFILE", result_path.name),
    ]
    baseml_control = fill_baseml_template(ctl_template, replacements)
    with open(tree_path, 'w') as fo:
        fo.write(f"{len(tree)} 1\n")
        fo.write(writetree(tree))
        fo.write('\n')
    with open(control_path, 'w') as fo:
        fo.write(baseml_control)

    proc_args = ["baseml", f"baseml.ctl"]
    with open(log_path, 'w') as fo:
        baseml_proc = subprocess.Popen(proc_args, cwd=output_path, stdout=fo)
        baseml_proc.wait()
    # print(result_path)
    result = best_result(baseml.parse_baseml_result(result_path))
    return result

def run_single_shift_baseml(tree, output_path, ctl_template, cleanup="delete"):
    model_assignments = single_shift_assignments(tree)
    best_info = None
    best_likelihood = float('-inf')
    output_path.mkdir(parents=True, exist_ok=True)
    for ix, model_tree in enumerate(model_assignments, start=1):
        sub_path = Path(f"{output_path}/single_shift_{ix}")
        sub_path.mkdir(exist_ok=True)
        result = run_baseml(tree, sub_path, ctl_template, writetree = lambda x : model_tree)
        if best_likelihood < result["log likelihood"]:
            best_info = result
            best_info["Path"] = sub_path
            best_likelihood = result["log likelihood"]
    if cleanup == "delete":
        for p in output_path.iterdir():
            if p != best_info["Path"]:
                subprocess.run(["rm", "-rf", p.absolute()])
    #elif cleanup == "compress":

    return best_info

def serial_single_shift_search(args):
    if args.seed: # ehhh... untested. a half measure at best
        random.seed(args.seed)

    with open(args.template, 'r') as fi:
        ctl_template = "".join(fi.read())
    ctl_template = ctl_template.replace("#SEQFILE", str(args.seq.absolute()))
    
    alignment = AlignIO.read(args.seq, format='phylip-relaxed')
    if args.start == 'nj':
        visited_trees = [nj_tree(alignment)]
    elif args.start == 'random':
        visited_trees = [random_tree(alignment)]
    else:
        raise NotImplementedError(f"{args.starting_tree} strategy not implemented :(")
    del alignment # was only used to get an initial estimate

    stepnum = 0
    best_tree = visited_trees[0]
    opath = Path(f"{args.output}/step_{stepnum}")
    result = run_single_shift_baseml(best_tree, opath, ctl_template)
    best_info = result
    if result:
        best_likelihood = result["log likelihood"]
    else:
        best_likelihood = float('-inf')
    
    # cursed
    ndigits = len(str(args.max_iter))

    print(f"{str(stepnum).zfill(ndigits)}:\t{best_likelihood}")
    iter_num = 1
    while iter_num <= args.max_iter:
        # strategy is to visit the neighbors of the highest likelihood tree visited
        # to try different strategies probably change this
        prev_best_tree = best_tree
        neighbors = nearest_neighbors(best_tree)
        for neighbor in neighbors:
            # Check if this neighbor has been visited already
            # Maybe there's a faster way of doing this?
            dont_skip = True
            for v_tree in visited_trees:
                if v_tree.compare(neighbor, unrooted=True)["norm_rf"] == 0:
                    dont_skip = False
                    break
            if dont_skip:
                stepnum += 1
                opath = Path(f"{args.output}/step_{stepnum}")
                visited_trees.append(neighbor)
                result = run_single_shift_baseml(neighbor, opath, ctl_template)
                if result and result["log likelihood"] > best_likelihood:
                    best_likelihood = result["log likelihood"]
                    best_tree = result["tree"]
                    best_info = result
                else:
                    # how to handle suboptimal results
                    # in this case, just delete them so they don't pollute the filesystem
                    subprocess.run(["rm", "-rf", opath.absolute()])
        if prev_best_tree.compare(best_tree, unrooted=True)["norm_rf"] == 0:
            print(f"{str(iter_num).zfill(ndigits)}:\tDid not find a better tree, stopping...")
            break
        print(f"{str(iter_num).zfill(ndigits)}:\t{best_likelihood}")
        iter_num += 1

    # cleanup
    for p in args.output.iterdir():
        if not best_info["Path"].is_relative_to(p):
            subprocess.run(["rm", "-rf", p.absolute()])
    print(f"best result: {result['tree'].write(format=9)}")
    return best_info["Path"]
