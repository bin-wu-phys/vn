This code is used to calculate v_n by numerically solving a kinetic transport as described in  arXiv:1803.02072. It is written in python 3.

This is version 9. In this version, the radial wavenumber is added such that  we can study the dependence of v_n on the radial profile. Here, we works with pairs (dn, ln) with n the theta harmonic number. If ln is zeor or ignore, the initial condition is our previous one.

Run on condor requies the following modification:
--line 1: #!/afs/cern.ch/user/b/biwu/python/bin/python3.7 to your python folder
--line 5: change the value of work_path to your desired folder for the result
--Do not forget create log, output and error folder in the folder containing init.py.
