import re
from pathlib import Path

def verified_file_path(file_path):
    if type(file_path) != Path:
        file_path = Path(file_path)
    assert file_path.is_file()
    return file_path

def read_all_lines(file_path):
    file_path = verified_file_path(file_path)
    with open(file_path, 'r') as fi:
        result = fi.readlines()
    return result

def find_line(lines, content, condition = lambda x, y : x in y):
    for i,line in enumerate(lines):
        if condition(content, line):
            return i

def find_lines(lines, content, condition = lambda x, y : x in y):
    valid_lines = []
    for i,line in enumerate(lines):
        if condition(content, line):
            valid_lines.append(i)
    return valid_lines

def first_instance_before_line(lines, line_number, content, condition = lambda x, y : x in y):
    for i in range(line_number, -1, -1):
        if condition(content, lines[i]):
            return i
    return None

def first_instance_after_line(lines, line_number, content, condition = lambda x, y : x in y):
    for i in range(line_number, len(lines)):
        if condition(content, lines[i]):
            return i
    return None

