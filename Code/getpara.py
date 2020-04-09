#!/usr/bin/env python3

import sys

def decode(s, ss, es=''):
    """
    Read out the value starting from ss to es in s.
    :param s: the string
    :param ss: the starting sting
    :param es: the end string
    :return: the value in between
    """
    if es == '':
        idx_s =  s.index(ss) + len(ss)
        idx_e = idx_s
        for c in s[idx_s:]:
            if c.isalpha():
                break
            else:
                idx_e = idx_e + 1
        res = s[idx_s:idx_e]
        if res[-1] == '.':
            res = res[:-1]
    else:
        res = s[s.index(ss) + len(ss):s.index(es)]

    return res

if len(sys.argv) == 2:
    if sys.argv[1][-1] == 'h' or sys.argv[1][-4:] == 'help':
        print(
            """
            It runs as:
            "./init.py para=value" with para given by
            t0, dt, tMax : initial time, time-step, the maximum time to run
            err : the maximum error used for adaptive time step
            d1, d2, d3, d4, d5, d6, g : the perturbation for nth harmonica and \hat{gamma}
            p1, p2, p3, p4, p5, p6: the nth reaction plane angle in radian
            ln: the radial wavenumber corresponding to n
            aQ : aQ=True, False for old, new initial gaussian condition
            rMax, nr, vMax, nv, nph, nth : parameters for setting up the lattice
            procs: the number of processes for parallel computation
            shiftQ: shiftQ=True or False to determine whether the initial profile shift its origin or not
                    shiftQ=True by default and it is only useful for the case with ic=file.
            or
            "./init.py -f (file name for parameters) (line for the parameters)"
            It also takes the input file and runs like:
            "./init.py ic=file fn=file name para=value"
            Here, one can:
            ns, dns, dos: integer, float, float
                switch on one modes ns=(some integer) and specify the amplitude dns=(some value)
                and the amplitude for other modes dos=(some value), zero by default.
            nsMax : int
                switch off all the harmonic modes larger than nsMax.
            Here, if both (ns, dns) and nsMax are given, (ns, dns) has the priority.    
            smoothQ : bool
                switch on smooth function with R_smooth.
            Rs : float
                the smearing radius (R_smooth) for the smooth function
            continuumQ, nthc : bool, int
                If confinuumQ = True, populate the initial F calculated with a large number of nth = nthc.
            switchboard : a string of integers
                Each digit tells us which fourier mode is kept. Here, I totally rule out the possibility with n >9.
            sba: switch board amplide is a string of floats separated by |
                if not empty, each float gives the amplitude of the corresponding digit in switchboard.
                
            One can also decides on whether to output F by using outputFQ=True or False.
            """
        )
        quit()

#default values
outputFQ = True

nphd = 40 # n_phi
rMaxd = 3.0 # r_max
nrd = 31 # n_r
vMaxd = 10.0 # v_max
nvd = 21 # n_v
procs = 6 # number of processors used
R_smooth = 0.5 # smearing radius
smoothQ = False
nthc = 60
continuumQ = False

switchboard = ''#is empty by default.
switchboard_amps = ''

icd = 'gaussian'
ictype=icd
fname = ''

#default values for nth reaction plan3 angle
psi1d = 0.0
psi2d = 0.0
psi3d = 0.0
psi4d = 0.0
psi5d = 0.0
psi6d = 0.0

#by default all the radial wavenumbers are zero, which means using our default radial profile
l1 = 0
l2 = 0
l3 = 0
l4 = 0
l5 = 0
l6 = 0

adtQ = True
shiftQ = True
t0 = 0.05
dtd = 0.05
dt = 0.05
dr = 0.1
d1 = 0.0
delta = 1.0
delta3 = 0.0
d4 = 0.0
d5 = 0.0
d6 = 0.0
aQ = False
g = 1.0
tMax = 4.0
rMax = rMaxd
nr = nrd
vMax = vMaxd
nv = nvd
nph = nphd
nth = 16
errd = 1.0e-4
err = 1.0e-3

#nth reaction plan3 angle
psi1 = 0.0
psi2 = 0.0
psi3 = 0.0
psi4 = 0.0
psi5 = 0.0
psi6 = 0.0

#initial condition from file, switch on one mode
ns = 0
dns = 1.0
dos = 0.0
nsMax = 0

# Read the parameters
if len(sys.argv) == 4:
    if sys.argv[1] == '-f':
        f = open(sys.argv[2])
        idx = int(sys.argv[3])
        argv = f.readlines()[idx]
        argvs = argv.replace('\n', '').split(' ')
        print(argvs)
    else:
        argvs = sys.argv[1:]
else:
    argvs = sys.argv[1:]

folder = ''
arguments = ''
for para in argvs:
    arguments = arguments + ' ' + para
    cmd = para.split('=')
    if cmd[0] == 'folder':
        folder = cmd[1]
    elif cmd[0] == 'g':
        g = float(cmd[1])
    elif cmd[0] == 'd1':
        d1 = float(cmd[1])
    elif cmd[0] == 'd2':
        delta = float(cmd[1])
    elif cmd[0] == 'd3':
        delta3 = float(cmd[1])
    elif cmd[0] == 'd4':
        d4 = float(cmd[1])
    elif cmd[0] == 'd5':
        d5 = float(cmd[1])
    elif cmd[0] == 'd6':
        d6 = float(cmd[1])
    elif cmd[0] == 'p1':
        psi1 = float(cmd[1])
    elif cmd[0] == 'p2':
        psi2 = float(cmd[1])
    elif cmd[0] == 'p3':
        psi3 = float(cmd[1])
    elif cmd[0] == 'p4':
        psi4 = float(cmd[1])
    elif cmd[0] == 'p5':
        psi5 = float(cmd[1])
    elif cmd[0] == 'p6':
        psi6 = float(cmd[1])
    elif cmd[0] == 'l1':
        l1 = int(cmd[1])
    elif cmd[0] == 'l2':
        l2 = int(cmd[1])
    elif cmd[0] == 'l3':
        l3 = int(cmd[1])
    elif cmd[0] == 'l4':
        l4 = int(cmd[1])
    elif cmd[0] == 'l5':
        l5 = int(cmd[1])
    elif cmd[0] == 'l6':
        l6 = int(cmd[1])
    elif cmd[0] == 'ns':
        ns = int(cmd[1])
    elif cmd[0] == 'dns':
        dns = float(cmd[1])
    elif cmd[0] == 'dos':
        dos = float(cmd[1])
    elif cmd[0] == 'nsMax':
        nsMax = int(cmd[1])
    elif cmd[0] == 'ic':
        ictype = cmd[1]
    elif cmd[0] == 'fn':
        fname = cmd[1]
    elif cmd[0] == 'switchboard':
        switchboard = cmd[1]
        if not switchboard.isalnum():
            print('switchboard has to be numeric.')
            exit()
    elif cmd[0] == 'sba':
        switchboard_amps = cmd[1].split('|')
    elif cmd[0] == 'aQ':
        if cmd[1] == 'True':
            aQ = True
        elif cmd[1] == 'False':
            aQ = False
        else:
            print("aQ can only be True or False.")
            exit()
    elif cmd[0] == 'adtQ':
        if cmd[1] == 'True':
            adtQ = True
        elif cmd[1] == 'False':
            adtQ = False
        else:
            print("adtQ can only be True or False.")
            exit()
    elif cmd[0] == 'smoothQ':
        if cmd[1] == 'True':
            smoothQ = True
        elif cmd[1] == 'False':
            smoothQ = False
        else:
            print("smoothQ can only be True or False.")
            exit()
    elif cmd[0] == 'shiftQ':
        if cmd[1] == 'True':
            shiftQ = True
        elif cmd[1] == 'False':
            shiftQ = False
        else:
            print("shiftQ can only be True or False.")
            exit()
    elif cmd[0] == 'outputFQ':
        if cmd[1] == 'True':
            outputFQ = True
        elif cmd[1] == 'False':
            outputFQ = False
        else:
            print("outputFQ can only be True or False.")
            exit()
    elif cmd[0] == 'continuumQ':
        if cmd[1] == 'True':
            continuumQ = True
        elif cmd[1] == 'False':
            continuumQ = False
        else:
            print("continuumQ can only be True or False.")
            exit()
    elif cmd[0] == 't0':
        t0 = float(cmd[1])
    elif cmd[0] == 'dt':
        dt = float(cmd[1])
    elif cmd[0] == 'dr':
        dr = float(cmd[1])
    elif cmd[0] == 'tMax':
        tMax = float(cmd[1])
    elif cmd[0] == 'Rs':
        R_smooth = float(cmd[1])
    elif cmd[0] == 'rMax':
        rMax = float(cmd[1])
    elif cmd[0] == 'nr':
        nr = int(cmd[1])
    elif cmd[0] == 'vMax':
        vMax = float(cmd[1])
    elif cmd[0] == 'nv':
        nv = int(cmd[1])
    elif cmd[0] == 'nth':
        nth = int(cmd[1])
    elif cmd[0] == 'nph':
        nph = int(cmd[1])
    elif cmd[0] == 'procs':
        procs = int(cmd[1])
    elif cmd[0] == 'err':
        err = float(cmd[1])
    else:
        print("'{}' is not recognized. Please run './init.py -h' to look it up.".format(cmd[0]))
        exit()

    if len(switchboard_amps) > 0:
        if len(switchboard_amps) != len(switchboard):
            print("If you want to provide the amplitudes for the switchboard, please provide enough!\n")
            exit()

# err=err/g

if ictype == icd:
    if folder == '':
        fn = 't0.{:.2f}.d1.{:.2f}.d2.{:.2f}.d3.{:.2f}.d4.{:.2f}.d5.{:.2f}.d6.{:.2f}.g.{:.2f}.nth.{}'.format(
                                                    t0, d1, delta, delta3, d4, d5, d6, g, nth, vMax, nv)
        # the ratational angle
        if psi1 != psi1d:
            fn = fn + '.p1.{:.1f}'.format(psi1)
        if psi2 != psi2d:
            fn = fn + '.p2.{:.1f}'.format(psi2)
        if psi3 != psi3d:
            fn = fn + '.p3.{:.1f}'.format(psi3)
        if psi4 != psi4d:
            fn = fn + '.p4.{:.1f}'.format(psi4)
        if psi5 != psi5d:
            fn = fn + '.p5.{:.1f}'.format(psi5)
        if psi6 != psi6d:
            fn = fn + '.p6.{:.1f}'.format(psi6)

        # add the radial wavenumber
        if l1 != 0:
            fn = fn + '.l1.{}'.format(l1)
        if l2 != 0:
            fn = fn + '.l2.{}'.format(l2)
        if l3 != 0:
            fn = fn + '.l3.{}'.format(l3)
        if l4 != 0:
            fn = fn + '.l4.{}'.format(l4)
        if l5 != 0:
            fn = fn + '.l5.{}'.format(l5)
        if l6 != 0:
            fn = fn + '.l6.{}'.format(l6)

        if vMax != vMaxd or nv != nvd:
            fn = fn + f'.vMax.{vMax}.nv.{nv}'
        if rMax != rMaxd or nr != nrd:
            fn = fn + f'.rMax.{rMax}.nr.{nr}'

        if nph != nphd:
            fn = fn + f'.nph.{nph}'

        if adtQ:
            if err !=errd:
                fn = fn + '.err.{:.4f}'.format(err)
        else:
            if dt != dtd:
                fn = fn + '.dt.{:.4f}'.format(dt)

        if aQ:
            fn = fn + '.old'
    else:
        t0 = float(decode(folder, 't0.', '.d1.'))
        d1 = float(decode(folder, '.d1.', '.d2.'))
        delta = float(decode(folder, '.d2.', '.d3.'))
        delta3 = float(decode(folder, '.d3.', '.d4.'))
        d4 = float(decode(folder, '.d4.', '.d5.'))
        d5 = float(decode(folder, '.d5.'))
        if '.d6.' in folder:
            d6 = float(decode(folder, '.d6.'))
        g = float(decode(folder, '.g.', '.nth.'))
        nth = int(decode(folder, '.nth.'))

        if '.p1.' in folder:
            psi1 = float(decode(folder, '.p1.'))
        if '.p2.' in folder:
            psi2 = float(decode(folder, '.p2.'))
        if '.p3.' in folder:
            psi3 = float(decode(folder, '.p3.'))
        if '.p4.' in folder:
            psi4 = float(decode(folder, '.p4.'))
        if '.p5.' in folder:
            psi5 = float(decode(folder, '.p5.'))
        if '.p6.' in folder:
            psi6 = float(decode(folder, '.p6.'))

        for nl in range(1, 7):
            if 'l{}'.format(nl) in folder:
                exec("l{} = int(decode(folder, '.l{}.'))".format(nl, nl))

        if '.vMax.' in folder:
            vMax = float(decode(folder, '.vMax.'))
        if '.nv.' in folder:
            nv = int(decode(folder, '.nv.'))
        if '.rMax.' in folder:
            rMax = float(decode(folder, '.rMax.'))
        if '.nr.' in folder:
            nr = int(decode(folder, '.nr.'))
        if '.nph.' in folder:
            nph = int(decode(folder, '.nph.'))

        if '.old' in folder:
            aQ = True

        fn = folder

elif ictype == 'file':
    fn = 't0.{:.2f}.fn.{}.g.{:.2f}.nth.{}'.format(t0, fname[:-4], g, nth)
    if continuumQ:
        fn = fn + f'.nthc.{nthc}'
    if ns != 0:
        fn = fn + '.ns.{}.dns.{:.2f}'.format(ns, dns)
        if dos !=0.0:
            fn = fn + '.dos.{:.2f}'.format(dos)
    elif nsMax != 0:
        fn = fn + '.nsMax.{}'.format(nsMax)
        if dos != 0.0:
            fn = fn + '.dos.{:.2f}'.format(dos)
    if vMax != vMaxd or nv != nvd:
        fn = fn + f'.vMax.{vMax}.nv.{nv}'
    if rMax != rMaxd or nr != nrd:
        fn = fn + f'.rMax.{rMax}.nr.{nr}'

    if nph != nphd:
        fn = fn + f'.nph.{nph}'
    if shiftQ:
        fn = fn + '.shifted'
    else:
        fn = fn + '.unshifted'

    if smoothQ:
        fn = fn + '.smoothed'
        fn = fn + '.Rs.{:.2f}'.format(R_smooth)

    if len(switchboard)  > 0:
        if len(switchboard_amps) == 0:
            fn = fn + '.sb.' + switchboard
        else:
            for idx_sb in range(len(switchboard)):
                fn = fn + '.d' + switchboard[idx_sb] + '.' + switchboard_amps[idx_sb]

    if adtQ:
        if err !=errd:
            fn = fn + '.err.{:.4f}'.format(err)
    else:
        if dt != dtd:
            fn = fn + '.dt.{:.4f}'.format(dt)

print(fn)
