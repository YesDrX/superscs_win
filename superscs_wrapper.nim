import superscs_api
import nimpy
import tables
import strformat

type SolutionResult = ref object of PyNimObjectExperimental
    x : seq[float]
    y : seq[float]
    s : seq[float]
    status : string
    info : Table[string, float]

proc version() :string {.exportpy.} =
    $scs_version()

proc update_stgs(data : ptr ScsData, stgs_kwargs : Table[string, float] = Table[string, float]()) =
    for key, value in stgs_kwargs.pairs:
        case key:
            of "normalize":
                data.stgs.normalize = value.cint
            of "scale":
                assert value > 0, "stgs.settings scale should be positive!"
                data.stgs.scale = value
            of "rho_x":
                assert value > 0, "stgs.settings rho_x should be positive!"
                data.stgs.rho_x = value
            of "max_time_milliseconds":
                data.stgs.max_time_milliseconds = value
            of "max_iters":
                data.stgs.max_iters = value.cint
            of "previous_max_iters":
                data.stgs.previous_max_iters = value.cint
            of "eps":
                assert value > 0, "stgs.settings eps should be positive!"
                data.stgs.eps = value
            of "alpha":
                assert value > 0, "stgs.settings alpha should be positive!"
                data.stgs.alpha = value
            of "cg_rate":
                assert value > 0, "stgs.settings cg_rate should be positive!"
                data.stgs.cg_rate = value
            of "verbose":
                data.stgs.verbose = value.cint
            of "warm_start":
                data.stgs.warm_start = value.cint
            of "do_super_scs":
                data.stgs.do_super_scs = value.cint
            of "k0":
                data.stgs.k0 = value.cint
            of "c_bl":
                data.stgs.c_bl = value
            of "k1":
                data.stgs.k1 = value.cint
            of "k2":
                data.stgs.k2 = value.cint
            of "c1":
                data.stgs.c1 = value
            of "sse":
                data.stgs.sse = value
            of "ls":
                data.stgs.ls = value.cint
            of "beta":
                data.stgs.beta = value
            of "sigma":
                data.stgs.sigma = value
            of "direction":
                data.stgs.direction = value.cint.ScsDirectionType
            of "thetabar":
                data.stgs.thetabar = value
            of "memory":
                data.stgs.memory = value.cint
            of "tRule":
                data.stgs.tRule = value.cint
            of "broyden_init_scaling":
                data.stgs.broyden_init_scaling = value.cint
            of "do_record_progress":
                data.stgs.do_record_progress = value.cint
            of "do_override_streams":
                data.stgs.do_override_streams = value.cint
            of "tolerance":
                scs_set_tolerance(data, value)
            of "anderson_settings":
                scs_set_anderson_settings(data, value.cint)
            else:
                raise newException(OSError, fmt"""
            Unknown optimization setting key : {key}
            Check available settings at https://github.com/twvacek/superscs/blob/bfe3e9f968ff6bbb14b5fb3b8e1ebb4f5f233a9a/include/scs.h#L290
            """)

proc get_solution(data : ptr ScsData, sol: ptr ScsSolution, info: ptr ScsInfo): (
    seq[float], # x
    seq[float], # y
    seq[float], # s
    string,
    Table[string, float]
    ) =
    var
        ptr_x = cast[ptr UncheckedArray[float]](sol.x)
        ptr_y = cast[ptr UncheckedArray[float]](sol.y)
        ptr_s = cast[ptr UncheckedArray[float]](sol.s)
        x = newSeq[float](data.n.int)
        y = newSeq[float](data.m.int)
        s = newSeq[float](data.m.int)
        status = $(cast[cstring](info.status.addr))
        rst_info : Table[string, float]
    
    for idx in 0 .. data.n-1:
        x[idx] = ptr_x[idx]
        y[idx] = ptr_y[idx]
    
    for idx in 0 .. data.m-1:
        s[idx] = ptr_s[idx]

    rst_info = {
        "iter" : info.iter.float,
        "statusVal" : info.statusVal.float,
        "pobj" : info.pobj.float,
        "dobj" : info.dobj.float,
        "resPri" : info.resPri.float,
        "resDual" : info.resDual.float,
        "resInfeas" : info.resInfeas.float,
        "resUnbdd" : info.resUnbdd.float,
        "setupTime" : info.setupTime.float,
        "solveTime" : info.solveTime.float,
        "linsys_total_solve_time_ms" : info.linsys_total_solve_time_ms.float,
        "allocated_memory" : info.allocated_memory.float
    }.toTable

    return (x,y,s,status,rst_info)

proc solve(
    rst_container : SolutionResult,
    n  : int,
    m  : int,
    nnz : int,
    Ax : openArray[float],
    Ai : openArray[int],
    Ap : openArray[int],
    b  : openArray[float],
    c  : openArray[float],
    cone_f: int = 0,
    cone_l: int = 0,
    q  : openArray[int] = @[],
    qsize: int = 0,
    s : openArray[int] = @[],
    ssize : int = 0,
    ep : int = 0,
    ed : int = 0,
    p : openArray[float] = @[],
    psize: int = 0,
    stgs_kwargs : Table[string, float] = {
        "eps" : 1e-6,
        "tolerance" : 1e-3,
        "memory" : 10.0
    }.toTable
    ) : void {.exportpy.} = 
    var
        data : ptr ScsData = scs_init_data()
        copy_c = newSeq[float](n)
        copy_b = newSeq[float](m)
        copy_Ax = newSeq[float](nnz)
        copy_Ai = newSeq[cint](nnz)
        copy_Ap = newSeq[cint](n+1)
        copy_p = newSeq[float](p.len)
        copy_q = newSeq[cint](q.len)
        copy_s = newSeq[cint](s.len)
        A = ScsAMatrix()
        cone = ScsCone()
        # copy_xx = newSeq[float](n)
        # copy_yy = newSeq[float](nnz)
        # copy_ss = newSeq[float](nnz)
        sol = scs_init_sol()
        info = scs_init_info()
    
    # sol.x = copy_xx[0].addr
    # sol.y = copy_yy[0].addr
    # sol.s = copy_ss[0].addr   

    # prepare data
    data.m = m.cint
    data.n = n.cint
    data.c = copy_c[0].addr
    data.b = copy_b[0].addr
    for idx in 0 .. n-1:
        copy_c[idx] = c[idx]
    for idx in 0 .. m-1:
        copy_b[idx] = b[idx]
    A.m = m.cint
    A.n = n.cint
    A.p = copy_Ap[0].addr
    A.i = copy_Ai[0].addr
    A.x = copy_Ax[0].addr
    for idx in 0 .. n:
        copy_Ap[idx] = Ap[idx].cint
    for idx in 0 .. nnz-1:
        copy_Ai[idx] = Ai[idx].cint
        copy_Ax[idx] = Ax[idx]
    data.A = A.addr

    scs_set_default_settings(data)
    update_stgs(data, stgs_kwargs)

    cone.ssize = ssize.cint
    cone.ed = ed.cint
    cone.ep = ep.cint
    cone.f = cone_f.cint
    cone.l = cone_l.cint
    cone.psize = psize.cint
    cone.qsize = qsize.cint
    if copy_p.len > 0:
        for idx in 0 .. copy_p.len - 1:
            copy_p[idx] = p[idx]
        cone.p = copy_p[0].addr
    if copy_s.len > 0:
        for idx in 0 .. copy_s.len - 1:
            copy_s[idx] = s[idx].cint
        cone.s = copy_s[0].addr
    if copy_q.len > 0:
        for idx in 0 .. copy_q.len - 1:
            copy_q[idx] = q[idx].cint
        cone.q = copy_q[0].addr
    
    discard scs(data, cone.addr, sol, info)

    var
        x : seq[float]
        y : seq[float]
        s : seq[float]
        rst_status : string
        rst_info : Table[string, float]
    
    (x,y,s,rst_status,rst_info) = get_solution(data, sol, info)

    # deallocShared(data.stgs) # seq in data is collected by nim GC
    # deallocShared(data)
    # scs_free_info(info)
    # scs_free_sol(sol)
    # deallocShared(info)
    # deallocShared(sol)

    rst_container.x = x
    rst_container.y = y
    rst_container.s = s
    rst_container.status = rst_status
    rst_container.info = rst_info

proc retrive_solution_status(self : SolutionResult) : string {.exportpy.} =
    self.status

proc retrive_solution_info(self : SolutionResult) : Table[string, float] {.exportpy.} =
    self.info

proc retrive_solution(self : SolutionResult, key : string) : seq[float] {.exportpy.} =
    case key
        of "x":
            return self.x
        of "y":
            return self.y
        of "s":
            return self.s

# when isMainModule:
#     for idx in 0 .. 10:
#         var
#             rst = SolutionResult()
        
#         solve(
#             rst,
#             n = 11,
#             m = 23,
#             nnz = 122,
#             Ax = @[1.0, -1.0, 0.006891499136395144, 0.21368342853400465, 0.04088481019471927, 0.0005673948788824094, 0.10801876980381059, 0.40650306285578985, 0.13048234478663606, -0.618235115888334, 0.44300774376221985, -0.9415718376174984, 1.0, -1.0, -0.04443933867653382, 0.12117410454967474, -0.0931272548597298, -0.02677465804707761, 0.3824690204589194, -0.12086718107967034, 0.1314821530556831, 0.18340853808104382, 0.647403608088848, 1.1850922010638936, 1.0, -1.0, 0.13391367754176228, 0.042786619022644835, 0.1256997585950353, 0.15127559953123687, 0.020368135523559, 0.2266603515289119, -0.22854859687140927, 0.5305915706531563, 0.3810562449390569, 0.4569722565978448, 1.0, -1.0, -0.057929109472956986, 0.08062784063746409, 0.21934068031758608, 0.03647977354130569, -0.14960159523914482, 0.1742332557865288, 0.6115725923324578, 0.666438308375037, -0.49787308569346606, 0.18047176641580698, 1.0, -1.0, -0.05234768948683697, -0.11548968671244478, 0.08142582560828789, 0.07923173833229143, 0.42870619875969895, 0.5985998864001585, -0.4022293306777772, 0.2309801370375362, -0.18781972141387962, -0.2704434717427616, 1.0, -1.0, -6.408535346501832e-05, -0.12400983053670471, 0.2200502897484133, 0.15488887354235248, 0.11857933791709883, -0.29651090506171857, 0.37480206090644, -0.7342996432377628, 0.6686493551719745, 0.18728512743198145, 1.0, -1.0, 0.036104340292928125, -0.07977534263332724, 0.00800951330603168, -0.46579363042542066, -0.02261400264523733, 0.5217918370452926, 0.34429892254794553, -0.25713302178773806, 0.08550244389936576, 0.4910864396005295, 1.0, -1.0, 0.047197099119830103, 0.013033619766102886, -0.09429127706412471, 0.2228294553764798, 0.18900042756745036, 0.14489546741287165, 0.38127715629153847, -0.5258170307348474, -1.2382889960533618, 0.44226843554657924, 1.0, -1.0, -0.018713551732816636, 0.09208506480476579, 0.18251354525797334, -0.14803342339833175, -0.05531020915979324, -0.15466001858027353, -0.6875615381671514, -0.5487384460380861, -0.6243508155366594, 0.6610162660441602, 1.0, -1.0, 0.04586831998587353, 0.017180564637982097, 0.07207918565982427, -0.2540059914649758, 0.39611041793202684, -0.5600364708417459, 0.1450470242885587, 0.3453333393307771, -0.4418964296232916, -0.6969075212888558, -1.0, 1.0],
#             Ap = @[0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 122],
#             Ai = @[0, 1, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 2, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 3, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 4, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 5, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 6, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 8, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 11, 12],
#             b = @[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#             c = @[-1.6243453636632417, -0.6117564136500754, -0.5281717522634557, -1.0729686221561705, -0.8654076293246785, -2.3015386968802827, -1.74481176421648, -0.7612069008951028, -0.31903909605709857, -0.2493703754774101, 0.21845094725845543],
#             cone_f = 1,
#             cone_l = 10,
#             ep = 0,
#             ed = 0,
#             q = @[12],
#             qsize = 1,
#             s = @[],
#             ssize = 0,
#             p = @[],        
#             psize = 0,
#             stgs_kwargs = {
#                 "verbose" : 0.0
#             }.toTable
#         )
#         echo "idx = ", idx
#         echo "x", rst.x
#         echo "y", rst.y
#         echo "s", rst.s
#         echo "info", rst.info