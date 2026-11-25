"""bxutil.py - small helpers shared by the BondX paper scripts (geometry of gas-phase dimers, running bondx_cli.js)."""
import json, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'P': 1.07, 'S': 1.05, 'CL': 1.02, 'BR': 1.20, 'I': 1.39, 'LI': 1.28}
HAL = ('CL', 'BR', 'I')

def read_xyz(fn):
    L = open(fn).read().split('\n'); n = int(L[0].split()[0]); el, X = [], []
    for l in L[2:2 + n]:
        p = l.split(); el.append(p[0].upper()); X.append([float(v) for v in p[1:4]])
    return el, np.array(X)

def fragments(el, X, tol=0.45):
    n = len(el); adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(X[i] - X[j]) < COV[el[i]] + COV[el[j]] + tol: adj[i].append(j); adj[j].append(i)
    lab = [-1] * n; k = 0
    for i in range(n):
        if lab[i] >= 0: continue
        st = [i]; lab[i] = k
        while st:
            a = st.pop()
            for b in adj[a]:
                if lab[b] < 0: lab[b] = k; st.append(b)
        k += 1
    return lab, adj

def six_ring(c, adj, el):
    """Return the carbon six-membered ring containing atom c, or None."""
    def dfs(path):
        if len(path) == 6: return path if path[0] in adj[path[-1]] else None
        for b in adj[path[-1]]:
            if el[b] == 'C' and b not in path:
                r = dfs(path + [b])
                if r: return r
        return None
    return dfs([c])

def analyse(el, X, monoA=None):
    # Strong halogen bonds (e.g. H3N...ClF) are short; lower the bond tolerance until the dimer splits into
    # two fragments whose sizes match the monomer files.
    nA = len(monoA[0]) if monoA is not None else None
    for tol in (0.45, 0.35, 0.25, 0.15, 0.05):
        lab, adj = fragments(el, X, tol)
        sizes = sorted(lab.count(k) for k in set(lab))
        if len(sizes) == 2 and (nA is None or nA in sizes): break
    if len(set(lab)) != 2: raise ValueError('expected two fragments, got %d' % len(set(lab)))
    best = None
    for prefer in (HAL, ('F',)):
        for x in range(len(el)):
            if el[x] not in prefer: continue
            for y in range(len(el)):
                if lab[y] == lab[x] or el[y] == 'H': continue
                d = np.linalg.norm(X[x] - X[y])
                if best is None or d < best[0]: best = (d, x, y)
        if best: break
    d, x, y = best
    a = min(adj[x], key=lambda j: np.linalg.norm(X[j] - X[x]))
    ring = six_ring(y, adj, el) if el[y] == 'C' else None
    if ring:
        cen = X[ring].mean(axis=0); d = np.linalg.norm(X[x] - cen)
        v = X[x] - cen; nrm = np.cross(X[ring[1]] - X[ring[0]], X[ring[2]] - X[ring[0]]); nrm /= np.linalg.norm(nrm)
        alpha = np.degrees(np.arccos(abs(v @ nrm) / np.linalg.norm(v)))
        yv = cen
    else:
        alpha = None; yv = X[y]
    u, w = X[a] - X[x], yv - X[x]
    th1 = angle_deg(u, w)
    return dict(X=x + 1, A=a + 1, Y='' if ring else y + 1, ring=','.join(str(i + 1) for i in ring) if ring else '', halogen=el[x], acc_el='pi' if ring else el[y],
                d=round(float(d), 4), theta1=round(float(th1), 2), alpha=round(float(alpha), 1) if alpha is not None else '', n_atoms=len(el),
                frag_donor=sum(1 for i in lab if i == lab[x]), frag_acceptor=sum(1 for i in lab if i != lab[x]))


def angle_deg(u, w):
    c = float(u @ w / np.linalg.norm(u) / np.linalg.norm(w))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, c)))))

def write_xyz(fn, el, X, comment=''):
    with open(fn, 'w') as f:
        f.write('%d\n%s\n' % (len(el), comment))
        for e, c in zip(el, X): f.write('%-3s %14.8f %14.8f %14.8f\n' % (e[0] + e[1:].lower(), *c))

def cli(*args):
    """Run bondx_cli.js and return its JSON output."""
    r = subprocess.run(['node', os.path.join(HERE, 'bondx_cli.js'), *[str(a) for a in args]], capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr)
    return json.loads(r.stdout)
