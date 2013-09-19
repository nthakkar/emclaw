#!/usr/bin/env python
# encoding: utf-8
# ddd
import numpy as np
# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
n_frames = 30
# ....... dimensions .............................................
x_lower = 0.0
x_upper = 100.0e-6                    # lenght [m]
y_lower = 0.0
y_upper = 10.0e-6                   # notice that for multilayer this is value will be over-written
# ........ material properties ...................................

# vacuum
eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
co = 1.0/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo = np.sqrt(mo/eo)
# material
mat_shape = 'epsle1'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered

# background refractive index and etas
eta      = np.ones([3])
bkg_n    = np.ones([2])

bkg_n[0] = np.sqrt(eta[0]*eta[2])
bkg_n[1] = np.sqrt(eta[1]*eta[2])


# if interface declare position
x_change = (x_upper-x_lower)/2.0
y_change = (y_upper-y_lower)/2.0

# set moving refractive index or gaussian2D parameters
rip_velocity = np.zeros([2,3])
rip_offset   = np.zeros([2,3])
rip_sigma    = np.zeros([2,3])
delta_eta    = np.zeros([3])

rip_offset[0,:].fill((x_upper-x_lower)/2.0)
rip_offset[1,:].fill((y_upper-y_lower)/2.0)
rip_sigma[0,:].fill((x_upper-x_lower)/25.0)
rip_sigma[1,:].fill((y_upper-y_lower)/25.0)
rip_sigma.fill(10e-6)
rip_sigma.fill(10e-6)
rip_sigma2 = rip_sigma**2

delta_eta = np.zeros([3])
delta_eta = 0.1*eta

# set multilayer parameters

# multilayered definition
n_layers = 2
layers = np.zeros([n_layers,9]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
layers[0,0] = 1.5
layers[0,1] = 1.5
layers[0,2] = 1.5
layers[0,3] = 10
layers[0,4] = 15e-9
layers[1,0] = 2.5
layers[1,1] = 2.5
layers[1,2] = 2.5
layers[1,3] = layers[0,3] - 1
layers[1,4] = 50e-9
N_layers = 5
if mat_shape=='multilayer':
    y_upper = N_layers*np.sum(layers[:,4])+layers[0,4]
    tlp = np.sum(layers[:,4])
    mlp = np.floor(tlp/1e-9)

# set non-linear parameters of the material
chi2 = chi3 = np.zeros( [3], order='F')

# ........ excitation - initial conditoons .......................

# pre-allocate arrays
ex_sigma     = np.ones([3])    # x,y,t
ex_offset    = np.zeros([3])
ex_amplitude = np.ones([3])
ex_kvector   = np.zeros([2])

# fill arrays and set respective values
ex_type   = 'off'
ex_lambda = 1e-6
ex_sigma[0:1] = 1.0*ex_lambda
ex_sigma[2]   = (y_upper-y_lower)/2.0
ex_offset[2]  = (y_upper-y_lower)/2.0

# post calculations
omega    = 2.0*np.pi*co/ex_lambda   # frequency
k        = 2.0*np.pi/ex_lambda      # k vector magnitude
ex_kvector[0] = k                   # propagation along the x-direction


# ........ pre-calculations for wave propagation .................

v = co/bkg_n.min()

# Grid - mesh settings
mx = np.floor(10*(x_upper-x_lower)/ex_lambda)
if mat_shape=='multilayer':
    my = np.floor((y_upper-y_lower)/1e-9)
else:
    my = np.floor(10*(y_upper-y_lower)/ex_lambda)

ddx = (x_upper-x_lower)/mx
ddy = (y_upper-y_lower)/my
ddt = dt=0.50/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)))
max_steps = 1000000
t_final = (x_upper-x_lower)/v
print t_final

# -------- GLOBAL FUNCTION DEFINITIONS --------------

# refractive index map definition function 
def etar(t,X,Y):
    """
    eta = etar(t,X,Y)

    This function returns the refractive index map based on general definitions set earlier,
    Gaussian cases support moving RIPs.
    
    x are the coordinate of the grid centers state.grid.e_j.centers, e_j = x 
         aux holds:
         0: epsilon
         1: mu
         2: epsilon_t
         3: mu_t
    """
    y,x = np.meshgrid(Y,X)
    eta_out = np.zeros( [6,len(X),len(Y)], order='F')

    if mat_shape=='gaussian1dx':

        u_x_eta1 = x - rip_velocity[0,0]*t - rip_offset[0,0]
        u_x_eta2 = x - rip_velocity[0,1]*t - rip_offset[0,1]
        u_x_eta3 = x - rip_velocity[0,2]*t - rip_offset[0,2]
        u_y_eta1 = y - rip_velocity[1,0]*t - rip_offset[1,0]
        u_y_eta2 = y - rip_velocity[1,1]*t - rip_offset[1,1]
        u_y_eta3 = y - rip_velocity[1,2]*t - rip_offset[1,2]

        u_eta1 = (u_x_eta1/rip_sigma[0,0])**2 + (u_y_eta1/rip_sigma[1,0])**2
        u_eta2 = (u_x_eta2/rip_sigma[0,1])**2 + (u_y_eta2/rip_sigma[1,1])**2
        u_eta3 = (u_x_eta3/rip_sigma[0,2])**2 + (u_y_eta3/rip_sigma[1,2])**2

        u_eta1_t = 2*((rip_velocity[0,0]*u_x_eta1)/(rip_sigma[0,0]**2) + (rip_velocity[1,0]*u_y_eta1)/(rip_sigma[1,0]**2))
        u_eta2_t = 2*((rip_velocity[0,1]*u_x_eta2)/(rip_sigma[0,1]**2) + (rip_velocity[1,1]*u_y_eta2)/(rip_sigma[1,1]**2))
        u_eta3_t = 2*((rip_velocity[0,2]*u_x_eta3)/(rip_sigma[0,2]**2) + (rip_velocity[1,2]*u_y_eta3)/(rip_sigma[1,2]**2))

        eta_out[0,:,:] = delta_eta[0]*np.exp(-u_eta1) + eta[0]
        eta_out[1,:,:] = delta_eta[1]*np.exp(-u_eta2) + eta[1]
        eta_out[2,:,:] = delta_eta[2]*np.exp(-u_eta3) + eta[2]
        
        eta_out[3,:,:] = u_eta1_t*delta_eta[0]*np.exp(-u_eta1)
        eta_out[4,:,:] = u_eta2_t*delta_eta[1]*np.exp(-u_eta2)
        eta_out[5,:,:] = u_eta3_t*delta_eta[2]*np.exp(-u_eta3)
    elif mat_shape=='homogeneous':
        eta_out[0,:,:] = eta[0]
        eta_out[1,:,:] = eta[1]
        eta_out[2,:,:] = eta[2]
    elif mat_shape=='interfacex':
        eta_out[0,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta_out[1,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta_out[2,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    elif mat_shape=='interfacey':
        eta_out[0,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
        eta_out[1,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
        eta_out[2,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
    elif mat_shape=='epsle1':
        eta_out[0,:,:] = 1.0*(x<x_change) + 0.5*(x>=x_change)*(x<=(3.0*(x_upper-x_lower)/4.0)) + 1.0*(x>(3.0*(x_upper-x_lower)/4.0))
        eta_out[1,:,:] = 1.0*(x<x_change) + 0.5*(x>=x_change)*(x<=(3.0*(x_upper-x_lower)/4.0)) + 1.0*(x>(3.0*(x_upper-x_lower)/4.0))
        eta_out[2,:,:] = 1.0*(x<x_change) + 0.5*(x>=x_change)*(x<=(3.0*(x_upper-x_lower)/4.0)) + 1.0*(x>(3.0*(x_upper-x_lower)/4.0))
    elif mat_shape=='multilayer':
        for n in range(0,N_layers):
            yi = n*tlp
            for m in range(0,n_layers):
                if m==0:
                    eta_out[0,:,:] = layers[m,0]*(yi<y)*(y<=yi+layers[m,3])
                    eta_out[1,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
                    eta_out[2,:,:] = layers[m,2]*(yi<y)*(y<=yi+layers[m,3])
                else:
                    eta_out[0,:,:] = layers[m,0]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
                    eta_out[1,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
                    eta_out[2,:,:] = layers[m,2]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])


        eta_out[0,:,:] = layers[0,0]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])
        eta_out[1,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])
        eta_out[2,:,:] = layers[0,2]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])  

    return eta_out

def update_aux(solver,state):
    grid = state.grid
    y = state.grid.y.centers
    x = state.grid.x.centers
    t = state.t
    state.aux = setaux(t,x,y)

    return state

#   next function might be redundant since it already exists as deltan  
def setaux(t,x,y):
    aux = np.empty( [6,len(y),len(x)], order='F')
    aux = etar(t,x,y)
    return aux

def setaux_lower(state,dim,t,auxbc,num_ghost):
    grid = state.grid
    X = grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    Y = grid.y.centers_with_ghost(num_ghost)[:num_ghost]
    t = state.t
    auxbc[:,:num_ghost,:num_ghost] = etar(t,X,Y)
    return auxbc

def setaux_upper(state,dim,t,auxbc,num_ghost):
    grid = state.grid
    X = grid.x.centers_with_ghost(num_ghost)[-num_ghost:]
    Y = grid.y.centers_with_ghost(num_ghost)[-num_ghost:]
    t = state.t
    auxbc[:,-num_ghost:,-num_ghost:] = etar(t,X,Y)
    return auxbc

def scattering_bc(state,dim,t,qbc,num_ghost):
    """
    EM scattering boundary conditions with three components in  TM-mode Ex, Ey, Hz.
    """
    grid = state.grid
    X = grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    Y = grid.y.centers_with_ghost(num_ghost)[:num_ghost]
    ts = state.t
    y,x = np.meshgrid(Y,X)
    t0 = 0.0
    aux_left_bc = etar(t,X,Y)
    pulseshape = np.zeros( [len(X),len(Y)], order='F')
    harmonic = np.zeros( [len(X),len(Y)], order='F')

    if ex_type=='plane':
        pulseshape = 1.0
        harmonic = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='gauss-beam':
        pulseshape = np.exp(-(y - ex_yoff)**2/ex_y_sig**2)
        harmonic = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='gauss_pulse':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2 - (y - ex_yoff - ex_vy*(ts-t0))**2/ex_y_sig**2)
        harmonic = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='plane_pulse':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2)
        harmonic = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='simple_pulse2D':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2 - (y - ex_yoff - ex_vy*(ts-t0))**2/ex_y_sig**2)
        harmonic = 1.0
    elif ex_type=='simple_pulse2D_x':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2)
        harmonic = 1.0
    elif ex_type=='off':
        pulseshape = 0.0
        harmonic = 0.0

    qbc[0,:num_ghost,:num_ghost] = 0.0
    qbc[1,:num_ghost,:num_ghost] = zo*ex_amplitude[1]*pulseshape*harmonic
    qbc[2,:num_ghost,:num_ghost] = ex_amplitude[1]*pulseshape*harmonic

    return qbc

def qinit(state):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    
    if ex_type=='off':
        grid = state.grid
        X = grid.x.centers
        Y = grid.y.centers
        y,x = np.meshgrid(Y,X)
        dd1 = (x_upper-x_lower)/5.0
        dd2 = y_upper-y_lower
        sdd = 1e-6
        r2 = (x-dd1/2.0)**2 #+ (y-dd2/2.0)**2
        state.q[0,:,:] = 0.0
        state.q[1,:,:] = zo*np.exp(-r2/(sdd**2))
        state.q[2,:,:] = 1.0*np.exp(-r2/(sdd**2))
    else:
        state.q[0,:,:] = 0.0
        state.q[1,:,:] = 0.0
        state.q[2,:,:] = 0.0
    
    return state

# -------- MAIN SCRIPT --------------

def em2D(kernel_language='Fortran',before_step=False,iplot=False,htmlplot=False,use_petsc=True,save_outdir='./_test_le1',solver_type='sharpclaw'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw


#   Solver settings
    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D()
        solver.dimensional_split=False
        solver.limiters = pyclaw.limiters.tvd.MC
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D()
        solver.num_waves = 2
        solver.lim_type = 4
        solver.interpolation_order = 4


    solver.dt_initial= ddt
    solver.max_steps = max_steps
    import maxwell_2d
    solver.rp = maxwell_2d
    solver.fwave = True
    solver.cfl_max = 1.5
    solver.cfl_desired = 0.4
    solver.dt_variable = True
    print 'setup information:'
    print 'v_wave=',v
    print 'x_lim=',x_upper,' t_f=',t_final 
    print 'mx=',mx,'dx=',ddx,'my=',mx,'dy=',ddy, 'dt=',ddt,'N_max=',max_steps
    print 'lambda=',ex_lambda,'freq=',omega
    if before_step:
        print 'update aux'
        solver.call_before_step_qach_stage = 1
        solver.before_step = update_aux
    


#   define number of waves (eqn) and aux (eps,mu)
    num_eqn = 3
    num_aux = 6

#   abstract domain and state setup
    x_dime  = pyclaw.Dimension('x',x_lower,x_upper,mx)
    y_dime  = pyclaw.Dimension('y',y_lower,y_upper,my)
    domain  = pyclaw.Domain([x_dime,y_dime])
    
    state   = pyclaw.State(domain,num_eqn,num_aux)
    
    grid = state.grid
    X    = grid.x.centers
    Y    = grid.y.centers
    tini = state.t
    state.aux = etar(tini,X,Y)
    state.problem_data['dx']    = x_dime.delta
    state.problem_data['dy']    = y_dime.delta
    state.problem_data['chi2']  = chi2
    state.problem_data['chi3']  = chi3
    state.problem_data['eo']    = eo
    state.problem_data['mo']    = mo
    state.problem_data['co']    = co
    state.problem_data['zo']    = zo
    state.problem_data['vac1']  = eo
    state.problem_data['vac2']  = eo
    state.problem_data['vac3']  = mo

#   Boundary conditions
#   solver.user_bc_lower = scattering_bc
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap

    solver.user_aux_bc_lower = setaux_lower
    solver.user_aux_bc_upper = setaux_upper
    solver.aux_bc_lower[0] = pyclaw.BC.custom
    solver.aux_bc_upper[0] = pyclaw.BC.custom
    solver.aux_bc_lower[1] = pyclaw.BC.custom
    solver.aux_bc_upper[1] = pyclaw.BC.custom

#   Initial solution
    qinit(state)


#   controller
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = t_final
    claw.num_output_times = n_frames
    claw.solver = solver
    claw.solution = pyclaw.Solution(state,domain)
    claw.outdir = save_outdir
    claw.write_aux_always = True
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=save_outdir,file_format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=save_outdir,file_format=claw.output_format)

    return claw


if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em2D)
    


