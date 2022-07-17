@sim.register('avg')
def avg(model, modelparams, simparams, miscparams):
    range_a0 = 0.5
    range_a1 = 3.0
    nsteps_a = 26

    range_b0 = 1.0
    range_b1 = 0.1
    nsteps_b = 7

    range_c0 = 0
    range_c1 = -60
    nsteps_c = 61

    range_d0 = 29
    range_d1 = 37
    nsteps_d = 9

    iterlist = [[], [], [], []]
    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    range_c = linspace(range_c0, range_c1, nsteps_c)
    range_d = linspace(range_d0, range_d1, nsteps_d)

    for a in range(nsteps_a):
        for b in range(nsteps_b):
            for i in range_d:
                range_e0 = 0 + i
                range_e1 = 15 + i
                nsteps_e = 61

                range_e = list(linspace(range_e0, range_e1, nsteps_e))

                for j in range(nsteps_c):
                    iterlist[3].append(round(range_e[j], 3))
                    iterlist[2].append(round(range_c[j], 3))
                    iterlist[1].append(round(range_b[b], 3))
                    iterlist[0].append(round(range_a[a], 3))

    models = [model(earheight=ai, msh2=bi, rotorangle=ci, alpha=di, I0=250) for ai, bi, ci, di in
              zip(iterlist[0], iterlist[1], iterlist[2], iterlist[3]])
    with Pool() as pool:
        res = pool.map(execute_model, models)

    result = {'Torque': list(res)}

    with open(ModelDir.DATA / f'avgfullasl_100.json', 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=True)

    return result