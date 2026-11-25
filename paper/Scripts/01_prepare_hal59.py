#!/usr/bin/env python3
"""01_prepare_hal59.py
Builds the small-molecule benchmark set from the HAL59 subset of GMTKN55 (Goerigk et al. 2017).
HAL59 = 45 XB51 dimers (Kozuch & Martin 2013) + 14 X40 dimers (Rezac, Riley & Hobza 2012) with CCSD(T)/CBS
binding energies (kcal/mol, positive = bound).
Input : a clone of https://github.com/grimme-lab/GMTKN55 (folder HAL59)
Output: data/hal59/<name>.xyz and data/hal59/hal59.csv with the halogen-bond atoms X, A and Y (1-based) or the
        acceptor ring, the acceptor element and the reference energy.
Usage : python3 01_prepare_hal59.py /path/to/GMTKN55/HAL59 data/hal59
"""
import os, re, sys, csv
import numpy as np
from bxutil import read_xyz, analyse, write_xyz

def main(src, dst):
    os.makedirs(dst, exist_ok=True); rows = []
    for line in open(os.path.join(src, '.res')):
        if not line.startswith('$tmer'): continue
        p = line.split(); ref = float(p[-1]); spec = p[1]
        m = re.match(r'\{([^,}]+),', spec)
        if m: name = m.group(1); subset = 'XB51'; monoA = re.match(r'\{[^,]+,([^,]+),', spec).group(1)
        else:
            m = re.match(r'(\S+)\{A,B,AB\}', spec); name = m.group(1) + 'AB'; subset = 'X40'; monoA = m.group(1) + 'A'
        el, X = read_xyz(os.path.join(src, name, 'struc.xyz'))
        info = analyse(el, X, read_xyz(os.path.join(src, monoA, 'struc.xyz')))
        short = name[:-2] if subset == 'X40' else name
        out = os.path.join(dst, short + '.xyz')
        write_xyz(out, el, X, 'HAL59 %s %s ref_CCSD(T)/CBS = %.2f kcal/mol' % (subset, short, ref))
        if info['halogen'] == 'F': role = 'control: F-arene stacking, no sigma-hole donor'
        elif info['theta1'] < 90: role = 'control: halogen points away from the acceptor'
        else: role = 'halogen bond'
        rows.append(dict(name=short, subset=subset, ref_kcal=ref, role=role, **info))
    keys = list(rows[0].keys())
    with open(os.path.join(dst, 'hal59.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print('%d dimers written to %s' % (len(rows), dst))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/gmtkn/HAL59', sys.argv[2] if len(sys.argv) > 2 else 'data/hal59')
