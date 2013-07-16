from inout import IO
import os
import glob
import numpy as np
import shutil
import scipy as sp
import scipy.special as sps
from scipy.integrate import simps as sq
import numpy.linalg as npl
import matplotlib
matplotlib.use('Agg')
class Postprocess(object):
    def QuadField1D(self,solution,quad):
        #diff = npl.norm(solution.q[1]-q_old[1],1)
        self._x = solution.x.centers.copy()
        quad[0,solution.frame] = solution.t
        quad[1,solution.frame] = sq(solution.q[0]*solution.q[1],solution.x.centers)
        quad[2,solution.frame] = sq(np.sqrt(solution.q[0]**2 + solution.q[1]**2),solution.x.centers)
        return quad

    def QuadField2D(solution,quad):
        #diff = npl.norm(solution.q[1]-q_old[1],1)
        quad[0,solution.frame] = solution.t
        quad[1,solution.frame] = sq(np.sqrt(solution.q[1]*solution.q[2]**2 + (-solution.q[0]*solution.q[2])**2),solution.x.centers)
        quad[2,solution.frame] = sq(np.sqrt(solution.q[0]**2 + solution.q[1]**2 + solution.q[2]**2),solution.x.centers)    
        return quad

    def EMEnergy(self,solution):
        self.u[solution.frame,:] = (solution.q[0]*(solution.q[0]*solution.aux[0])+solution.q[1]*(solution.q[1]*solution.aux[1]))/2.0
        return self

    def EMEnergy_dt(self):
        #self.udt = self.u.copy()
        dt = self.quad[0,1]-self.quad[0,0]
        self.udt = (self.u[1:,:]-self.u[0:-1,:])/dt
        return self

    def Poyinting1D(self,solution):
        self.S = np.zeros([1,len(solution.x.centers)])
        self.S = solution.q[0]*solution.q[1]
        return self

    def Poyinting2D(self,solution):
        self.S = np.zeros([2,len(solution.x.centers),len(solution.y.centers)])
        self.S[0,:,:] = solution.q[1]*solution.q[2]
        self.S[1,:,:] = -solution.q[0]*solution.q[2]
        return self

    def FieldIntensity1D(self,solution):
        self.I = np.zeros([1,len(solution.x.centers)])
        self.I = np.sqrt(solution.q[0]**2 + solution.q[1]**2)
        return self

    def FieldIntensity2D(self,solution):
        self.I = np.zeros([1,len(solution.x.centers),len(solution.y.centers)])
        self.I[0,:,:] = np.sqrt(solution.q[0]**2 + solution.q[1]**2 + solution.q[2]**2)
        return self

    def RefIndex1D(self,solution):
        self.n = np.zeros([1,len(solution.x.centers)])
        self.n = np.sqrt(solution.aux[0]*solution.aux[1])
        return self

    def RefIndex2D(self,solution):
        self.n = np.zeros([2,len(solution.x.centers),len(solution.y.centers)])
        self.n[0,:,:] = np.sqrt(solution.aux[0]*solution.aux[2])
        self.n[1,:,:] = np.sqrt(solution.aux[1]*solution.aux[2])
        return self   

    def get_num_frames(self):
        n = 0
        file_name = self.file_prefix+'.ptc*'
        for dir_name in glob.glob(os.path.join(self.path,file_name)):
            head, tail = os.path.split(dir_name)
            if tail.find('info')==-1:
                n = n+1
        return int(n)

    def PlotQ1D(self,solution):

        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        fig_path = self.outdir
        fig_base_name = 'frame'+str(solution.frame).zfill(4)


        self.Poyinting1D(solution)
        self.FieldIntensity1D(solution)
        self.RefIndex1D(solution)

        from matplotlib import pylab as plt
        plt.close('all')

        plt.figure()
        plt.subplot(211)
        plt.plot(solution.x.centers/1e-6,self.I,'b-')
        title  = '$\mathscr{I}$ and $n$ at $t='+str(solution.t/1e-15)+'fs$'
        ylabel = '$\mathscr{I}$ ($a.u.$)'
        plt.title(title)
        plt.ylabel(ylabel)

        plt.subplot(212)
        plt.plot(solution.x.centers/1e-6,self.n,'g-')
        ylabel = '$n$ ($a.u.$)'
        xlabel = '$x$ ($\mu m$)'
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        fig_name = fig_base_name+'fig0.png'
        plt.draw()
        plt.savefig(os.path.join(fig_path,fig_name),dpi=240)
        plt.close()

        plt.figure()
        plt.subplot(211)
        plt.plot(solution.x.centers/1e-6,self.S,'b-')
        title  = '$\mathscr{S}$ and $n$ at $t='+str(solution.t/1e-15)+'fs$'
        ylabel = '$\mathscr{S}$ ($a.u.$)'
        plt.title(title)
        plt.ylabel(ylabel)

        plt.subplot(212)
        plt.plot(solution.x.centers/1e-6,self.n,'g-')
        ylabel = '$n$ ($a.u.$)'
        xlabel = '$x$ ($\mu m$)'
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        fig_name = fig_base_name+'fig1.png'
        plt.draw()
        plt.savefig(os.path.join(fig_path,fig_name),dpi=240)
        plt.close()

        plt.figure()
        plt.subplot(211)
        plt.plot(solution.x.centers/1e-6,solution.q[0],'b-')
        title  = '$\Psi_1$ and $\eta_1$ at $t='+str(solution.t/1e-15)+'fs$'
        ylabel = '$\Psi_1$ ($a.u.$)'
        plt.title(title)
        plt.ylabel(ylabel)

        plt.subplot(212)
        plt.plot(solution.x.centers/1e-6,solution.aux[0],'g-')
        ylabel = '$\eta_1$ ($a.u.$)'
        xlabel = '$x$ ($\mu m$)'
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        fig_name = fig_base_name+'fig2.png'
        plt.draw()
        plt.savefig(os.path.join(fig_path,fig_name),dpi=240)
        plt.close()

        plt.figure()
        plt.subplot(211)
        plt.plot(solution.x.centers/1e-6,solution.q[1],'b-')
        title  = '$\Psi_2$ and $\eta_e$ at $t='+str(solution.t/1e-15)+'fs$'
        ylabel = '$\Psi_2$ ($a.u.$)'
        plt.title(title)
        plt.ylabel(ylabel)

        plt.subplot(212)
        plt.plot(solution.x.centers/1e-6,solution.aux[1],'g-')
        ylabel = '$\eta_2$ ($a.u.$)'
        xlabel = '$x$ ($\mu m$)'
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        fig_name = fig_base_name+'fig3.png'
        plt.draw()
        plt.savefig(os.path.join(fig_path,fig_name),dpi=240)
        plt.close()

        plt.close('all')

    def PlotQ2D(self,solution):
        pass

    def post_frame_calc(self,frame):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        solution.frame = frame
        solution.read_aux = True
        solution.read_petsc()
        if solution.num_dim==1:
            self.QuadField1D(solution,self.quad)
            self.EMEnergy(solution)
        elif solution.num_dim==2:
            self.QuadField2D(solution,self.quad)

    def PlotFrames(self):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        for frame in self.plot_frames:
            solution.frame = frame
            solution.read_aux = True
            solution.read_petsc()
            if solution.num_dim==1:
                self.PlotQ1D(solution)
            if solution.num_dim==2:
                self.PlotQ2D(solution)

    def _PlotFrame(self,frame):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        solution.frame = frame
        solution.read_aux = True
        solution.read_petsc()
        if solution.num_dim==1:
            self.PlotQ1D(solution)
        if solution.num_dim==2:
            self.PlotQ2D(solution)

    def _SaveCalculation(self):
        if not os.path.isdir(os.path.join(self.outdir,'arrays')):
            os.makedirs(os.path.join(self.outdir,'arrays'))
        
        from array import array
        postcalc_path = os.path.join(self.outdir,'arrays')
        calc_file     = os.path.join(postcalc_path,'quads')
        output_file1  = open(os.path.join(postcalc_path,'quads'),'wb')
        output_file2  = open(os.path.join(postcalc_path,'u'),'wb')
        output_file3  = open(os.path.join(postcalc_path,'udt'),'wb')

        # out_array1 = array('d',self.quad)
        # out_array2 = array('d',self.u)
        # out_array3 = array('d',self.udt)

        self.quad.tofile(output_file1)
        self.u.tofile(output_file2)
        self.udt.tofile(output_file3)

        output_file1.close()
        output_file2.close()
        output_file3.close()
        # np.savetxt(os.path.join(postcalc_path,'postcalc.txt'),self.quad.transpose())
        # np.savetxt(os.path.join(postcalc_path,'u.txt'),self.u.transpose())
        # np.savetxt(os.path.join(postcalc_path,'udt.txt'),self.udt.transpose())

    def PlotFieldQuad(self):
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        fig_path = self.outdir

        from matplotlib import pylab as plt

        plt.close('all')
        plt.figure()
        quad = self.quad.copy()
        plt.plot(quad[0,:]/1e-15,quad[1,:],'g--',label='$\mathscr{S}$')
        plt.plot(quad[0,:]/1e-15,quad[2,:],'r--',label='$\mathscr{I}$')
        xlabel = 'time ($fs$)'
        ylabel = '$\mathscr{S}$, $\mathscr{I}$ ($a.u.$)'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(frameon=False)
        plt.draw()
        plt.savefig(os.path.join(fig_path,'si_vs_t.png'),dpi=240)
        plt.close()
        plt.figure()
        plt.plot(quad[0,:]/1e-15,quad[3,:],'g--',label='$\partial_t \mathscr{S}$')
        plt.plot(quad[0,:]/1e-15,quad[4,:],'r--',label='$\partial_t \mathscr{I}$')
        xlabel = 'time ($fs$)'
        ylabel = '$\partial_t$ $\mathscr{S}$, $\mathscr{I}$ ($fs^{-1}$)'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(frameon=False)
        plt.draw()
        plt.savefig(os.path.join(fig_path,'dtsi_vs_t.png'),dpi=240)
        plt.close()

        plt.figure()
        t = np.linspace(self.quad[0,0],self.quad[0,-1],len(self.quad[0,:]))
        x = self._x.copy()
        #print t.shape, x.shape,self.u.shape
        T,X = np.meshgrid(t,x)
        #print T[0::5,:].shape,X[0::5,:].shape,self.u[0::5,:].transpose().shape
        u = self.u.transpose()
        plt.pcolor(X[::10,:],T[::10,:],u[::10,:],shading='interp')
        xlabel = '$x$ ($\mu m$)'
        ylabel = 'time ($fs$)'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.colorbar()
        plt.draw()
        plt.savefig(os.path.join(fig_path,'u.png'),dpi=240)
        plt.close()
        plt.close('all')

    def get_num_cells(self):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        solution.frame = 0
        solution.read_aux = False
        solution.read_petsc()
        self.__setattr__('num_cells',solution.x.centers.size)
        self.__setattr__('num_dim',solution.num_dim)

    def CalcFieldQuad(self):
        if self.get_all_frames:
            self.__setattr__('num_frames',self.get_num_frames())
            self.u = np.zeros([self.num_frames,self.num_cells])
            self.quad = np.zeros([5,self.num_frames]) 
            for nframe in range(0,self.num_frames):
                self.post_frame_calc(nframe)
            self.quad[3,1:] = (self.quad[1,1:]-self.quad[1,0:-1])/(self.quad[0,1:]-self.quad[0,0:-1])
            self.quad[4,1:] = (self.quad[2,1:]-self.quad[2,0:-1])/(self.quad[0,1:]-self.quad[0,0:-1]) 
            self.EMEnergy_dt()
        else:
            self.__setattr__('num_frames',1)
            self.quad = np.zeros([3,self.num_frames]) 
            self.post_frame_calc(frame)

        if self.save_calc:
            self._SaveCalculation()

        if self.plot_calc:
            self.PlotFieldQuad()

    def _DeleteOldFolder(self):
        if os.path.isdir(os.path.join(self.path,self.outdir)):
            shutil.rmtree(os.path.join(self.path,self.outdir))

    def __init__(self,path='./',file_prefix='claw',get_all_frames=True,frame=0,plot=False,save_calc=True,plot_calc=True,outdir='_postcalc',clean_outdir=False):
        
        self.__setattr__('plot_frames', np.zeros([1]))
        self.__setattr__('path',path)
        self.__setattr__('file_prefix',file_prefix)
        self.__setattr__('outdir',os.path.join(path,outdir))
        self.__setattr__('get_all_frames',get_all_frames)
        self.__setattr__('save_calc',save_calc)
        self.__setattr__('plot_calc',plot_calc)
        self.__setattr__('_cleanold',clean_outdir)
        if self._cleanold:
            self._DeleteOldFolder()
        self.get_num_cells()

if __name__ == "__main__":
    import sys
    import Postprocess
    path = sys.argv[1]
    pst = Postprocess()
    # print 'going to path:', path
    # print 'number of frames:', num_frames
    # quad = np.zeros([3,num_frames])
    # file_name = os.join.path(path,'quad.txt')
    # for i in range(0,num_frames+1):
    #     print i
    #     sol = IO()
    #     pst = Postprocess()
    #     sol.path = path
    #     sol.frame = i
    #     sol.read_petsc()
    #     sol.q_to_vtk()
    #     pst.postcalc(sol,quad)
    #     q_old = sol.q.copy()

    # np.savetxt(file_name,quad)
    # os.remove('petclaw.log')
    # os.remove('pyclaw.log')
    # os.remove('inout.pyc')