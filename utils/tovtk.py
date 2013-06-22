from inout import IO
import os
import numpy as np
import numpy.linalg as npl

def postcalc(solution,q_old):
    diff = npl.norm(solution.q[1]-q_old[1],1)
    print diff
    return diff

if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    num_frames = int(sys.argv[2])
    print 'going to path:', path
    print 'number of frames:', num_frames
    for i in range(0,num_frames+1):
        print i
        sol = IO()
        sol.path = path
        sol.frame = i
        sol.read_petsc()
        sol.q_to_vtk()
        if i>0:
            postcalc(sol,q_old)
        q_old = sol.q

    os.remove('petclaw.log')
    os.remove('pyclaw.log')
    os.remove('inout.pyc')