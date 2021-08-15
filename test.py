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

import superscs_wrapper as _superscs
rst = _superscs.SolutionResult()

for i in range(10):
    print("idx = ", i)
    _superscs.solve(
        rst,
        n = 11,
        m = 23,
        nnz = 122,
        Ax = [1.0, -1.0, 0.006891499136395144, 0.21368342853400465, 0.04088481019471927, 0.0005673948788824094, 0.10801876980381059, 0.40650306285578985, 0.13048234478663606, -0.618235115888334, 0.44300774376221985, -0.9415718376174984, 1.0, -1.0, -0.04443933867653382, 0.12117410454967474, -0.0931272548597298, -0.02677465804707761, 0.3824690204589194, -0.12086718107967034, 0.1314821530556831, 0.18340853808104382, 0.647403608088848, 1.1850922010638936, 1.0, -1.0, 0.13391367754176228, 0.042786619022644835, 0.1256997585950353, 0.15127559953123687, 0.020368135523559, 0.2266603515289119, -0.22854859687140927, 0.5305915706531563, 0.3810562449390569, 0.4569722565978448, 1.0, -1.0, -0.057929109472956986, 0.08062784063746409, 0.21934068031758608, 0.03647977354130569, -0.14960159523914482, 0.1742332557865288, 0.6115725923324578, 0.666438308375037, -0.49787308569346606, 0.18047176641580698, 1.0, -1.0, -0.05234768948683697, -0.11548968671244478, 0.08142582560828789, 0.07923173833229143, 0.42870619875969895, 0.5985998864001585, -0.4022293306777772, 0.2309801370375362, -0.18781972141387962, -0.2704434717427616, 1.0, -1.0, -6.408535346501832e-05, -0.12400983053670471, 0.2200502897484133, 0.15488887354235248, 0.11857933791709883, -0.29651090506171857, 0.37480206090644, -0.7342996432377628, 0.6686493551719745, 0.18728512743198145, 1.0, -1.0, 0.036104340292928125, -0.07977534263332724, 0.00800951330603168, -0.46579363042542066, -0.02261400264523733, 0.5217918370452926, 0.34429892254794553, -0.25713302178773806, 0.08550244389936576, 0.4910864396005295, 1.0, -1.0, 0.047197099119830103, 0.013033619766102886, -0.09429127706412471, 0.2228294553764798, 0.18900042756745036, 0.14489546741287165, 0.38127715629153847, -0.5258170307348474, -1.2382889960533618, 0.44226843554657924, 1.0, -1.0, -0.018713551732816636, 0.09208506480476579, 0.18251354525797334, -0.14803342339833175, -0.05531020915979324, -0.15466001858027353, -0.6875615381671514, -0.5487384460380861, -0.6243508155366594, 0.6610162660441602, 1.0, -1.0, 0.04586831998587353, 0.017180564637982097, 0.07207918565982427, -0.2540059914649758, 0.39611041793202684, -0.5600364708417459, 0.1450470242885587, 0.3453333393307771, -0.4418964296232916, -0.6969075212888558, -1.0, 1.0],
        Ap = [0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 122],
        Ai = [0, 1, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 2, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 3, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 4, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 5, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 6, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 8, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 11, 12],
        b = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        c = [-1.6243453636632417, -0.6117564136500754, -0.5281717522634557, -1.0729686221561705, -0.8654076293246785, -2.3015386968802827, -1.74481176421648, -0.7612069008951028, -0.31903909605709857, -0.2493703754774101, 0.21845094725845543],
        cone_f = 1,
        cone_l = 10,
        ep = 0,
        ed = 0,
        q = [12],
        qsize = 1,
        s = [],
        ssize = 0,
        p = [],
        psize = 0,
        stgs_kwargs = {'verbose' : 0}
    )

    print(rst.retrive_solution_status())