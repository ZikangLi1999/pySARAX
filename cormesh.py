"""
Correct the axial mesh of LAVENDER

@Author: Zikang LI
@Date: 2023-4-9

The current axial layer generation mechanism in pySARAX is incorrect.
This script is used to densify the axial layer in LAVENDER.
It can be used in command line with two arguments:
    --file: the path of lavender.inp
    --height: the maximum tolerated height of mesh cell
like `python cormesh.py --file /path/to/lavender.inp --height 0.1`
"""
import os
import os.path
from math import ceil

def cormesh(lavender_file: str, max_cell_height: float) -> str:
    if not os.path.exists(lavender_file) or not os.path.isfile(lavender_file):
        raise RuntimeError(f'{lavender_file} does not exist or is not a file.')
    
    if max_cell_height <= 0.:
        raise ValueError('max_cell_height should be > 0.')
    
    fa_layers = list()
    with open(lavender_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.startswith('layer'):
                line_splited = line.split()
                # num_layer = int(line_splited[1])
                height_layer = map(float, line_splited[2:])
            elif line.startswith('FA_type'):
                fa_layers.append(
                    list(map(
                        lambda x: list(map(int, x.split('*'))) if '*' in x else [1, int(x)],
                        line.split()[2:]
                    ))
                )
    
    # Densify the mesh cell that is too coarse
    height_layer_dense = list()
    factors_dense = list()
    for height in height_layer:
        if height <= max_cell_height:
            height_layer_dense.append(height)
            factors_dense.append(0)
            continue
        
        factor = ceil(height / max_cell_height)
        factors_dense.append(factor - 1)
        
        height_dense = round(height / factor, 4)
        for _ in range(factor):
            height_layer_dense.append(height_dense)

    trg_path = os.path.join(os.path.split(lavender_file)[0], 'lavender_mesh{:.1f}'.format(max_cell_height))
    if not os.path.exists(trg_path):
        os.mkdir(trg_path)

    trg_file = os.path.join(trg_path, 'lavender.inp')
    with open(lavender_file, 'r', encoding='utf-8') as src:
        with open(trg_file, 'w', encoding='utf-8') as trg:
            fa_idx = -1
            for line in src.readlines():
                if line.startswith('layer'):
                    trg.write('{:<16}'.format('layer'))
                    trg.write('{:<16d}'.format(len(height_layer_dense)))
                    # trg.write(' '.join(map('{:.4f}'.format, height_layer_dense)))
                    trg.write(print_list(height_layer_dense))
                    trg.write('\n')
                    continue

                if not line.startswith('FA_type'):
                    trg.write(line)
                    continue
                
                # Write header of FA_type
                trg.write('{:<16}'.format('FA_type'))
                fa_idx += 1
                trg.write('{:<15d}'.format(fa_idx+1))

                # Write layer_mat_id
                cnt = 0
                for num, mat in fa_layers[fa_idx]:
                    factor = sum(factors_dense[:cnt+num]) - sum(factors_dense[:cnt])
                    cnt += num
                    if factor > 0:
                        num += factor
                    if num > 1:
                        trg.write(' {:d}*{:d}'.format(num, mat))
                    else:
                        trg.write(' {:d}'.format(mat))
                trg.write('\n')
    
    print('Mesh densified from {:d} to {:d}'.format(len(factors_dense), len(factors_dense)+sum(factors_dense)))
    print('Result output to [{}]'.format(trg_file))


def print_list(ls: list) -> str:
    if len(ls) == 0:
        return str()
    
    slist = list()
    prev = ls[0]
    cnt = 0

    for elem in ls:
        if elem == prev:
            cnt += 1
            continue
        if cnt == 1:
            slist.append('{:.4f}'.format(prev))
        else:
            slist.append('{:d}*{:.4f}'.format(cnt, prev))
        prev = elem
        cnt = 1
    
    # Do not forget the last one
    slist.append('{:.4f}'.format(prev))
    
    return ' '.join(slist)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--file', default=os.path.join(os.getcwd(), 'lavender.inp'), type=str)
    parser.add_argument('--height', default=0.1, type=float)
    args = parser.parse_args()
    cormesh(lavender_file=args.file, max_cell_height=args.height)
