import sys
import numpy as np

def theta(a, b):
    if a == '-' or b == '-' or a != b:   # gap or mismatch
        return -1
    elif a == b:                         # match
        return 1

def score_matrix(seq1, seq2):
    """
    return score matrix and map(each score from which direction)
    0: diagnosis
    1: up
    2: left
    """
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    mat = np.zeros(shape=[len(seq1), len(seq2)], dtype = int)
    matmap = np.zeros(shape=[len(seq1), len(seq2)], dtype = int)

    for i,p in enumerate(seq1):
        for j,q in enumerate(seq2):
            if i == 0:        # first row, gap in seq1
                mat[i,j] = -j
                matmap[i,j] = 2
                if j == 0:
                    matmap[i,j] = -1  # top left
                continue
            if j == 0:        # first column, gap in seq2
                mat[i,j] = -i
                matmap[i,j] = 1
                continue
            ul = mat[i-1, j-1] + theta(p, q)     # from up-left, 0
            l  = mat[i-1, j]   + theta(p, '-')   # from left, 1, gap in seq1
            u  = mat[i, j-1]   + theta('-', q)   # from up, 2, gap in seq2
            picked = max([ul,l,u])
            mat[i,j] = picked
            matmap[i,j] = [ul, l, u].index(picked)
    return mat, matmap

def traceback1(seq1, seq2, mat):
    """
    trace back score matrix to find the best path and return path code
    -!- this can not solve multiple equally optimal alignment -!-
    """
    path_code = ''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1)-1, len(seq2)-1
    while i > 0 or j > 0:
        if (i > 0 and j > 0) and mat[i, j] == mat[i-1, j-1] + theta(seq1[i], seq2[j]):
            path_code = '0' + path_code
            i -= 1
            j -= 1
        elif (i > 0 and mat[i, j] == mat[i-1, j] + theta(seq1[i], '-')):
            path_code = '1' + path_code
            i -= 1
        elif (j > 0 and mat[i, j] == mat[i, j-1] + theta('-' , seq2[j])):
            path_code = '2' + path_code
            j -= 1
    return path_code

def traceback2(seq1, seq2, matmap):
    '''find one optimal traceback path, from map matrix, return path code'''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = matmap[i,j]
        if direction == 0:
            i = i-1
            j = j-1
            path_code = '0' + path_code
        elif direction == 1:
            i = i-1
            path_code = '1' + path_code
        elif direction == 2:
            j = j-1
            path_code = '2' + path_code
    return path_code

def print_m(seq1, seq2, m):
    """ print score matrix or map matrix"""
    seq1 = '-' + seq1; seq2 = '-' + seq2
    print()
    print(' '.join(['%3s' % i for i in ' '+seq2]))
    for i, p in enumerate(seq1):
        line = [p] + [m[i,j] for j in range(len(seq2))]
        print(' '.join(['%3s' % i for i in line]))
    print()
    return

def format_align_result(seq1, seq2, path_code):
    '''
    return pair alignment result string from
    path code: 0 for match, 1 for gap in seq1, 2 for gap in seq2
    '''
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '1':
            align1 = align1 + seq1[0]
            align2 = align2 + '-'
            seq1 = seq1[1:]
            middle = middle + ' '
        elif p == '2':
            align1 = align1 + '-'
            align2 = align2 + seq2[0]
            seq2 = seq2[1:]
            middle = middle + ' '
        elif p == '0':
            align1 = align1 + seq1[0]
            align2 = align2 + seq2[0]
            if seq1[0] == seq2[0]:
                middle = middle + '|'
            else:
                middle = middle + ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]

    return '\n  ' + align1 + '\n  ' + middle + '\n  ' + align2 + '\n'

def main():
    seq1, seq2 = map(str.upper, sys.argv[1:3])
    print('seq1: %s' % seq1)
    print('seq2: %s' % seq2)

    mat, matmap = score_matrix(seq1, seq2)
    print_m(seq1, seq2, mat)
    print_m(seq1, seq2, matmap)

    path_code = traceback1(seq1, seq2, mat)
    print('trace back path code: ', path_code)
    print('Alignment:')
    print(format_align_result(seq1, seq2, path_code))

    print(' ', traceback2(seq1, seq2, matmap))

if __name__ == '__main__':
    main()
