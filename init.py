#!/afs/cern.ch/user/b/biwu/python/bin/python3.7
#!/usr/bin/env python3

import sys
work_path = '/afs/cern.ch/user/b/biwu/vn/version9'
sys.path.insert(0, work_path)

if __name__ == '__main__':

    from Code.initcond import *
    from Code.kintran import *
    import numpy as np
    import timeit
    import os
    from Code.getpara import *

    # create the directory for the results if needed
    # dirres = 'results'
    dirres = work_path

    if dirres !='':
        if not os.path.exists(dirres):
            os.mkdir(dirres)
        dirres = dirres + '/'

    finfo = dirres + fn + '/' + fn + '.info.dat'

    if os.path.exists(dirres + fn):
        Fs = [f for f in os.listdir(dirres + fn) if '.kt.' in f]
    else:
        os.mkdir(dirres + fn)
        Fs = []

    # t to save
    tsave = np.concatenate((np.linspace(0.25, 1.5, 6),
                            np.linspace(2.0, tMax, int(2.0 * (tMax - 2.0) + 1)))
                           )

    if len(Fs) == 0:
        # initialize

        # initilaize
        latt = Lattice(('r tilde', 0.0, rMax, nr),
                       ('vz tilde', 0.0, vMax, nv),
                       ('phir', nph),
                       ('theta', nth))

        ob = Observable(t0, latt)
        if ictype == 'gaussian':
            ic = InitCond(ictype, t0, (1, d1, psi1*np.pi/180.0, l1), (2, delta, psi2*np.pi/180.0, l2),
                                    (3, delta3, psi3*np.pi/180.0, l3), (4, d4, psi4*np.pi/180.0, l4),
                                    (5, d5, psi5*np.pi/180.0, l5), (6, d6, psi6*np.pi/180.0, l6))
            if aQ:
                ic.setAncient(aQ)
        elif ictype == 'file':
            if ns != 0:
                ic = InitCond(ictype, t0, (ns, dns))
                if dos != 0.0:
                    ic.set_delta_of_other_modes(dos)
            elif nsMax != 0:
                ic = InitCond(ictype, t0, (nsMax,))
            else:
                ic = InitCond(ictype, t0)
            ic.setFile(dirres + fname)
            ic.set_shiftQ(shiftQ)
            if smoothQ:
                ic.setSmooth(smoothQ)
                ic.setR_smooth(R_smooth)
        C = Kernel(g=g, ob=ob, dr=dr)
        kt = KinTran(t0, ic.initialize(latt), C, dt=dt, err=err)
        if procs > 1:
            kt.kern.ob.setProcs(procs)
        if not adtQ:
            kt.set_adaptive_dtQ(adtQ)
        
        if ictype == 'gaussian':
            parainfo = """The parameters:
              t0 = {}, dt = {}, dr = {}, d1 = {}, d2 = {}, d3 = {}, d4 = {}, d5 = {}, d6 = {}, p1 = {}, p2 = {}, p3 = {}, p4 = {}, p5 = {}, p6 = {},
            g = {}, aQ = {}, err = {}, tMax = {}, rMax = {}, nr = {}, vMax = {}, nv = {}, nph = {}, nth = {}, procs = {}\n
            """.format(ob.t0, kt.dt, C.dr, ic.dn[0][1], ic.dn[1][1], ic.dn[2][1], ic.dn[3][1], ic.dn[4][1], ic.dn[5][1],
                        ic.dn[0][2], ic.dn[1][2], ic.dn[2][2], ic.dn[3][2], ic.dn[4][2], ic.dn[5][2],
                        C.g, ic.ancientQ, kt.err, tMax, ob.latt.lattice[ob.idx_r][-1],
                        len(ob.latt.lattice[ob.idx_r]), ob.latt.vzmax,
                        len(ob.latt.lattice[ob.idx_v]),
                        len(ob.latt.lattice[ob.idx_ph]),
                        len(ob.latt.lattice[ob.idx_th]),
                        ob.procs
                    )
        elif ictype == 'file':
            parainfo = """The parameters:
              t0 = {}, dt = {}, dr = {}, g = {}, err = {}, fname = {}
              tMax = {}, rMax = {}, nr = {}, vMax = {}, nv = {}, nph = {}, nth = {}, procs = {}\n
            """.format(ob.t0, kt.dt, C.dr, C.g, kt.err, fname,
                   tMax, ob.latt.lattice[ob.idx_r][-1],
                   len(ob.latt.lattice[ob.idx_r]), ob.latt.vzmax,
                   len(ob.latt.lattice[ob.idx_v]),
                   len(ob.latt.lattice[ob.idx_ph]),
                   len(ob.latt.lattice[ob.idx_th]),
                   ob.procs
                   )

        print(parainfo)

        f = open(finfo, 'w')
        f.write(parainfo)
        f.close()

        vn = [ob.calcvn(t0, kt.F)]
        fnamevn = dirres + fn + '/' + fn + '.vn'

        filename = dirres + fn + '/' + fn + '.t.{:.2f}.kt'.format(kt.t)
        kt.save(filename)
        np.savez_compressed(fnamevn, vn=np.array(vn))

    else:
        #load F
        tf = 0.0
        fnstart = ''
        for f in Fs:
            #read out the time from the file name
            tt = float(f[f.index('.t.') + 3:f.index('.kt.')])
            if  tt > tf:
                tf = tt
                fnstart = f
        kt = KinTran()
        fnstart = dirres + fn + '/' + fnstart
        kt.load(fnstart)
        print('F is loaded from {} at t = {}.\n'.format(fnstart, kt.t))
        if procs > 1:
            kt.kern.ob.setProcs(procs)
        if not adtQ:
            kt.set_adaptive_dtQ(adtQ)

        parainfo = """The parameters:
                  t0 = {}, t = {}, dt = {}, dr = {}, g = {}, tMax = {},
                  rMax = {}, nr = {}, vMax = {}, nv = {}, nph = {}, nth = {}, procs = {}\n
            """.format(kt.kern.ob.t0, kt.t, kt.dt, kt.kern.dr, kt.kern.g,
                       tMax, kt.kern.ob.latt.lattice[kt.kern.ob.idx_r][-1],
                       len(kt.kern.ob.latt.lattice[kt.kern.ob.idx_r]),
                       kt.kern.ob.latt.vzmax, len(kt.kern.ob.latt.lattice[kt.kern.ob.idx_v]),
                       len(kt.kern.ob.latt.lattice[kt.kern.ob.idx_ph]),
                       len(kt.kern.ob.latt.lattice[kt.kern.ob.idx_th]),
                       kt.kern.ob.procs
                       )

        print(parainfo)

        #load vn
        fvns = [f for f in os.listdir(dirres + fn) if '.vn.' in f]
        if len(fvns) > 0:
            fnamevn = dirres + fn + '/' + fvns[-1]
            vn = np.load(fnamevn)['vn'].tolist()
            print('vn is loaded from {}.\n'.format(fnamevn))
        else:
            vn = [ob.calcvn(t0, kt.F)]
            fnamevn = dirres + fn + '/' + fn + '.vn'

        #shorten the time steps for saving the data
        idx = 0
        for i in range(len(tsave)):
            if kt.t > tsave[i] or tsave[i] - kt.t < 1e-3:
                idx = i
        tsave = tsave[idx + 1:]
        print(tsave)

    # tsave = np.linspace(1.0, tMax, int(tMax))
    start = timeit.default_timer()

    for ts in tsave:
            while ts - kt.t >1e-3:
                kt.nextTimeTo(ts)
            vn.append(kt.kern.ob.calcvn(kt.t, kt.F))
            np.savez_compressed(fnamevn, vn=np.array(vn))
            filename = dirres + fn + '/' + fn + '.t.{:.2f}.kt'.format(kt.t)
            kt.save(filename)
            cmptime = timeit.default_timer() - start
            f = open(finfo, 'a')
            f.write('For t0 = {:.2f} to t = {:.2f}, it takes {:.0f} hours {:.0f} minutes {:.2f} seconds.\n'.format(
                t0, kt.t, cmptime // 3600, (cmptime % 3600) // 60, (cmptime % 3600) % 60))
            f.close()
