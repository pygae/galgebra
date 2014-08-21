#references.py
import re
book = 'bookGA.tex'
include = re.compile(r'\\include\{[\w]*\}')
label = re.compile(r'\\label\{[a-zA-Z0-9_]*\}')
label_dict = {}
sub_lst = []

book_file = open(book, 'r')
book_str  = book_file.read()
book_file.close()
include_lst= include.findall(book_str)
chapter_lst = []
for chapter in include_lst:
    lbrc = chapter.find('{')
    rbrc = chapter.find('}')
    chapter_lst.append(chapter[lbrc+1:rbrc]+'.tex')
print chapter_lst

chapter_num = 1

def label_equation(eq_str):
    global chapter_num, eq_num, label_dict
    new_label = r'\label{EQ'+str(chapter_num)+'_'+str(eq_num)+'}'
    eq_num += 1
    labels = label.findall(eq_str)

    if len(labels) ==1:
        new_eq_str = eq_str.replace(labels[0],new_label)
        label_dict[labels[0]] = new_label
    else:
        new_eq_str = eq_str.replace(r'\begin{equation}',
                     r'\begin{equation}'+new_label)
    return(new_eq_str)

def label_align(align_str):
    global chapter_num, eq_num, label_dict
    align_lst = align_str.split(r'\\')
    print align_lst
    new_align_str = ''
    for align_line in align_lst[:-1]:
        if r'\nonumber' not in align_line:
            new_label = r'\label{EQ'+str(chapter_num)+'_'+str(eq_num)+'}'
            eq_num += 1
            labels = label.findall(align_line)
            if len(labels) == 1:
                align_line = align_line.replace(labels[0],new_label)
                label_dict[labels[0]] = new_label
            else:
                align_line += new_label
            new_align_str += align_line + r' \\'
        else:
            new_align_str += align_line + r' \\'
    if r'\nonumber' not in align_lst[-1]:
        new_label = r'\label{EQ'+str(chapter_num)+'_'+str(eq_num)+'}'
        eq_num += 1
        labels = label.findall(align_lst[-1])
        if len(labels) == 1:
            align_line = align_lst[-1].replace(labels[0],new_label)
            label_dict[labels[0]] = new_label
        else:
            align_line = align_lst[-1].replace('\n\\end{align}',new_label+'\n\\end{align}')
        new_align_str += align_line
    else:
        new_align_str += align_line
    return(new_align_str)

for chapter in chapter_lst[1:]:
    chapter_file = open(chapter)
    chapter_str  = chapter_file.read()
    chapter_file.close()

    tmp_str = chapter_str.replace('\\be\n','\\begin{equation}\n')
    tmp_str = tmp_str.replace('\\be\\','\\begin{equation}\\')
    tmp_str = tmp_str.replace('\\be ','\\begin{equation} ')
    tmp_str = tmp_str.replace('\\ee\n','\\end{equation}\n')
    tmp_str = tmp_str.replace('\\ee ','\\end{equation} ')

    begin_equation = tmp_str.find(r'\begin{equation}')
    begin_align = tmp_str.find(r'\begin{align}')

    eq_num = 1

    while begin_equation > -1 or begin_align > -1:
        if begin_equation > -1 and begin_align > -1:
            if begin_equation < begin_align:
                end_equation = tmp_str.find(r'\end{equation}',begin_equation)+14
                new_tmp_str = label_equation(tmp_str[begin_equation:end_equation])
                sub_lst.append((new_tmp_str,begin_equation,end_equation))
                begin_equation = tmp_str.find(r'\begin{equation}',end_equation)
            else:
                end_align = tmp_str.find(r'\end{align}',begin_align)+11
                new_tmp_str = label_align(tmp_str[begin_align:end_align])
                sub_lst.append((new_tmp_str,begin_align,end_align))
                begin_align = tmp_str.find(r'\begin{align}',end_align)
        elif begin_equation == -1:
                end_align = tmp_str.find(r'\end{align}',begin_align)+11
                new_tmp_str = label_align(tmp_str[begin_align:end_align])
                sub_lst.append((new_tmp_str,begin_align,end_align))
                begin_align = tmp_str.find(r'\begin{align}',end_align)
        else:
                end_equation = tmp_str.find(r'\end{equation}',begin_equation)+14
                new_tmp_str = label_equation(tmp_str[begin_equation:end_equation])
                sub_lst.append((new_tmp_str,begin_equation,end_equation))
                begin_equation = tmp_str.find(r'\begin{equation}',end_equation)

    chapter_num += 1
print label_dict
for sub in sub_lst:
    tmp_str = tmp_str.replace(tmp_str[sub[1]:sub[2]],sub[0],1)
f = open('tmp.tex','w')
for sub in sub_lst:
    f.write(sub[0])
    f.write('\n'+80*'*'+'\n')
f.close()


