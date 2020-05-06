""" A private module to parse multivector expressions """
import re
from typing import List


# The `#` character is included because we generate tokens containing it in
# _parse_paren.
_operand_re = r"([A-Za-z0-9\_\#]+)"

_operator_res = {
    '<>|': re.compile(r'(ARG(\||<|>)ARG)'.replace('ARG', _operand_re)),
    '^': re.compile(r'(ARG[\^]ARG([\^]ARG)*)'.replace('ARG', _operand_re)),
    '*': re.compile(r'(ARG[\*]ARG([\*]ARG)*)'.replace('ARG', _operand_re)),
}


def _contains_interval(interval1, interval2):  # interval1 inside interval2
    if interval1[0] > interval2[0] and interval1[1] < interval2[1]:
        return True
    else:
        return False


# counter to generate unique tokens
_parse_paren_calls = 0


def _parse_paren(line):
    global _parse_paren_calls
    _parse_paren_calls += 1

    if ('(' not in line) or (')' not in line):
        return [[[line]]]
    level = 0
    max_level = 0
    ich = 0
    paren_lst = []
    for ch in line:
        if ch == '(':
            level += 1
            paren_lst.append([level, ich])
        if ch == ')':
            if level < 1:
                raise ValueError('Mismathed Parenthesis in: ' + line + '\n')
            paren_lst.reverse()
            iparen = 0
            for elem in paren_lst:
                if elem[0] == level:
                    paren_lst[iparen].append(ich)
                    break
                iparen += 1
            paren_lst.reverse()
            level -= 1
        max_level = max(max_level, level)
        ich += 1
    if level != 0:
        raise ValueError('Mismatched Parenthesis in: ' + line + '\n')
    if max_level > 0:
        level_lst = []
        for _x in range(max_level + 1):
            level_lst.append([])
        for group in paren_lst:
            level_lst[group[0]].append(group[1:])
        ilevel = max_level
        while ilevel > 1:
            level = level_lst[ilevel]
            level_down = level_lst[ilevel - 1]
            igroup = 0
            for group in level:
                igroup_down = 0
                for group_down in level_down:
                    if _contains_interval(group, group_down):
                        level_lst[ilevel][igroup].append(igroup_down)
                    igroup_down += 1
                igroup += 1
            ilevel -= 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                token = '#' + str(_parse_paren_calls) + '_' + str(ilevel) + '_' + str(igroup) + '#'
                level_lst[ilevel][igroup].append(line[group[0]:group[1] + 1])
                level_lst[ilevel][igroup].append(token)
                igroup += 1
            ilevel += 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                group.append(group[-2])
                level_lst[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        ilevel = max_level
        while ilevel > 1:
            igroup = 0
            for group in level_lst[ilevel]:
                group_down = level_lst[ilevel - 1][group[2]]
                replace_text = group_down[-1].replace(group[-3], group[-2])
                level_lst[ilevel - 1][group[2]][-1] = replace_text
                igroup += 1
            ilevel -= 1
        for group in level_lst[1]:
            line = line.replace(group[2], group[3])
        ilevel = 1
        level_lst[0] = [[line]]
    return level_lst


def _unparse_paren(level_lst):
    line = level_lst[0][0][0]
    for level in level_lst[1:]:
        for group in level:
            new_string = group[-1]
            if new_string[:2] == '((' and new_string[-2:] == '))':
                new_string = new_string[1:-1]
            line = line.replace(group[-2], new_string)
    return line


def _sub_paren(s):
    string = s.group(0)
    return '(%s)' % string


def _add_paren(line, re_exprs):
    paren_flg = False
    if (line[0] == '(') and (line[-1] == ')'):
        paren_flg = True
        line = line[1:-1]
    if ('(' in line) or (')' in line):
        line_levels = _parse_paren(line)
        ilevel = 0
        for level in line_levels:
            igroup = 0
            for group in level:
                group[-1] = re.sub(re_exprs, _sub_paren, group[-1])
                line_levels[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        line = _unparse_paren(line_levels)
    else:
        line = re.sub(re_exprs, _sub_paren, line)
    if paren_flg:
        line = '(' + line + ')'
    return line


def validate_op_order(op_order: List[str]) -> None:
    if not all(op in _operator_res for op in op_order):
        raise ValueError("Illegal operator")


def parse_line(line: str, op_order: List[str]) -> str:
    line = line.replace(' ', '')
    level_lst = _parse_paren(line)
    ilevel = 0
    for level in level_lst:
        igroup = 0
        for group in level:
            string = group[-1]
            for op in op_order:
                string = _add_paren(string, _operator_res[op])
            level_lst[ilevel][igroup][-1] = string
            igroup += 1
        ilevel += 1
    line = _unparse_paren(level_lst)
    return line
