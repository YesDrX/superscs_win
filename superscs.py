import logging
import nimporter

nimporter.IGNORE_CACHE = True
nimporter.NimCompiler.NIM_CLI_ARGS = [
        '--opt:speed',
        '--parallelBuild:0',
        '--gc:arc',
        '--threads:on',
        '--app:lib',
        '-d:release',
        '-d:strip',
        '-d:lto',
        '-d:ssl',
        '--warning[ProveInit]:off']

# nimporter.NimCompiler.NIM_CLI_ARGS = [
#     '-gc:arc',
#     '--threads:on',
#     '--app:lib',
#     '-tlsEmulation:off'
#     ]

import superscs_win.superscs_wrapper as _superscs
from scipy import sparse
import numpy as np
import pickle 

logging.info(f"Using SuperSCS {_superscs.version()}")

DEFAULT_SETTINGS = {
    "do_super_scs" : 1.0
}

def _np_to_list(data):
    if isinstance(data, np.ndarray):
        return data.tolist()
    else:
        return data
    
def solve(probdata, cone, **kwargs):
    b = probdata['b']
    c = probdata['c']
    if sparse.issparse(b): b = b.todense()
    if sparse.issparse(c): c = c.todense()

    m = len(b)
    n = len(c)

    A = probdata['A']
    if not sparse.isspmatrix_csc(A): A = A.tocsc()

    Ax, Ai, Ap = A.data, A.indices, A.indptr
    nnz = len(Ax)

    stgs_kwargs = {
        k : float(v)
        for k,v in kwargs.items() if k in [
            "normalize",
            "scale",
            "rho_x",
            "max_time_milliseconds",
            "max_iters",
            "previous_max_iters",
            "eps",
            "alpha",
            "cg_rate",
            "verbose",
            "warm_start",
            "do_super_scs",
            "k0",
            "c_bl",
            "k1",
            "k2",
            "c1",
            "sse",
            "ls",
            "beta",
            "sigma",
            "direction",
            "thetabar",
            "memory",
            "tRule",
            "broyden_init_scaling",
            "do_record_progress",
            "do_override_streams",
            "tolerance",
            "anderson_settings"
        ]
    }

    for k,v in DEFAULT_SETTINGS.items():
        if k not in stgs_kwargs:
            stgs_kwargs[k] = v
    
    f = cone['f']
    l = cone['l']
    ep = cone['ep']
    ed = cone.get('ed',0)
    q = cone['q']
    p = cone['p']
    s = cone['s']
    q = q if q else []
    p = p if p else []
    s = s if s else []
    psize = len(p)
    qsize = len(q)
    ssize = len(s)    
    
    Ax = _np_to_list(Ax)
    Ai = _np_to_list(Ai)
    Ap = _np_to_list(Ap)
    b = _np_to_list(b)
    c = _np_to_list(c)
    p = _np_to_list(p)
    q = _np_to_list(q)

    # print("#################################################")
    # for item, val in {
    #     # 'n' : n,
    #     # 'm' : m,
    #     # 'nnz' : nnz,
    #     # 'Ax' : Ax,
    #     # 'Ai' : Ai,
    #     # 'Ap' : Ap,
    #     # 'b' : b,
    #     # 'c' : c,
    #     # 'f' : f,
    #     # 'l' : l,
    #     # 'p' : p,
    #     # 'psize' : psize,
    #     # 'q' : q,
    #     # 'qsize' : qsize,
    #     # 's' : s,
    #     # 'ssize' : ssize,
    #     # 'ep' : ep,
    #     # 'ed' : ed,
    #     'stgs_kwargs' : stgs_kwargs
    # }.items():
    #     print(f"{item} : [{type(val)}] : ")
    #     print(val)
    # print("#################################################")
    
    solution = _superscs.SolutionResult()
    
    _superscs.solve(
        solution,
        n = n,
        m = m,
        nnz = nnz,
        Ax = Ax,
        Ai = Ai,
        Ap = Ap,
        b  = b,
        c  = c,
        cone_f = f,
        cone_l = l,
        q = q,
        qsize = qsize,
        s = s,
        ssize = ssize,
        ep = ep,
        ed = ed,
        p = p,
        psize = psize,
        stgs_kwargs = stgs_kwargs
    )

    x = solution.retrive_solution("x")
    y = solution.retrive_solution("y")
    s = solution.retrive_solution("s")
    status = solution.retrive_solution_status()
    info = solution.retrive_solution_info()

    info['statusVal'] = int(info['statusVal'])
    info['iter'] = int(info['iter'])
    info['status'] = status

    result = {
        "info" : info,
        "x" : np.array(x),
        "y" : np.array(y),
        "s" : np.array(s)
    }

    # print(info)
    # pickle.dump(result, open("c:\\Users\\weixi\\super.pk","wb"))

    return result