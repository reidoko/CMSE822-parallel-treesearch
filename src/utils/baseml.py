from ._helpers import *
import numpy as np
import ete3

def parse_baseml_control(file_path):
    file_path = verified_file_path(file_path)

def parse_sub_rates(line):
    return np.array([float(x) for x in line.split()[1:]])

def parse_base_freqs(line):
    return np.array([float(x) for x in line.split()[3:]])

def parse_baseml_result(result_path):
    lines = read_all_lines(result_path)
    
    line_matcher = re.compile(r'(Node)[\s]+\#[\s]*[\d]+[\s]+[\w]+[\s]*\(blength')
    line_matcher = re.compile(r'(Node)[\s]+\#[\s]*[\d]+[\s]+[\w]+[\s]*\(blength')
    node_parameter_lines = [l for l in find_lines(lines, line_matcher, condition=lambda rgx, line : rgx.match(line) is not None)]

    if not node_parameter_lines:
        return []
    all_tree_parameters = []
    last_node_number = 0
    
    lnL_line = first_instance_before_line(lines, node_parameter_lines[0], 'lnL')
    lnL = float(lines[lnL_line].split(':')[3].split()[0])

    tree_line = first_instance_before_line(lines, node_parameter_lines[0], ';')
    
    tree_parameters = {
        'tree' : ete3.Tree(lines[tree_line].strip()),
        'log likelihood' : lnL
    }

    for npl in node_parameter_lines:
        line_nohash = lines[npl].split("#")[1].strip()
        node_number = int(line_nohash.split()[0])
        node_label = lines[npl].split()[1]
        if node_label == '(blength':
            node_label = None
        substitution_rates = parse_sub_rates(lines[npl+1])
        base_freqs = parse_base_freqs(lines[npl+3])
        node_parameters = {
            'node label' : node_label,
            'node #' : node_number,
            'substitution rates' : substitution_rates,
            'base frequencies' : base_freqs
        }
        if last_node_number > node_number:
            all_tree_parameters.append(tree_parameters)
            tree_line = first_instance_before_line(lines, npl, ';')
            lnL_line = first_instance_before_line(lines, npl, 'lnL')
            lnL = float(lines[lnL_line].split(':')[3].split()[0])
            tree_parameters = {'tree' : ete3.Tree(lines[tree_line].strip()), 'log likelihood' : lnL}
        tree_parameters[node_number] = node_parameters
        last_node_number = node_number
    all_tree_parameters.append(tree_parameters)
    return all_tree_parameters

