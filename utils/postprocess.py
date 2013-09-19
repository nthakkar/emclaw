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
        return self

    def get_num_dim(self):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        solution.frame = 50
        solution.read_aux = True
        solution.read_petsc()
        self.__setattr__('num_dim',solution.num_dim)
        if self.num_dim==1:
            self.RefIndex1D(solution)
        if self.num_dim==2:
            self.RefIndex2D(solution)

        self.__setattr__('cmax',self.co/self.n.min())
        self.__setattr__('cmin',self.co/self.n.max())

        return self

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

        return self

    def _DeleteOldFolder(self):
        if os.path.isdir(os.path.join(self.path,self.outdir)):
            shutil.rmtree(os.path.join(self.path,self.outdir))

    def Sampling(self):
            if not hasattr(self, 'num_dim'):
                self.__setattr__('num_dim',self.get_num_dim())
            if self.sample_all:
                if not hasattr(self, 'num_frames'):
                    self.__setattr__('num_frames',self.get_num_frames())

                self.q_sampled = np.zeros([len(self.sample_mode),self.num_dim+1+2*self.num_dim+1,self.num_frames])
                self.I_sampled = np.zeros([len(self.sample_mode),2*self.num_dim+1,self.num_frames])
                self.S_sampled = np.zeros([len(self.sample_mode),2*self.num_dim+1,self.num_frames])
                for nframe in range(0,self.num_frames):
                        self.PeakSampleFrame(nframe,nframe)
            else:
                samples = len(self.sample_frames)
                self.q_sampled = np.zeros([len(self.sample_mode),self.num_dim+1+2*self.num_dim+1,samples])
                self.I_sampled = np.zeros([len(self.sample_mode),2*self.num_dim+1,self.samples])
                self.S_sampled = np.zeros([len(self.sample_mode),2*self.num_dim+1,self.samples])
                for i,nframe in  enumerate(self.sample_frames):
                    self.PeakSampleFrame(nframe,i)

            if self.plot_sample:
                self.PlotSampling()

    def PlotSampling(self):
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        fig_path = self.outdir

        from matplotlib import pylab as plt

        plt.close('all')
        for k,mode in enumerate(self.sample_mode):
            print 'plotting ', mode
            plt.figure()
            quad = self.q_sampled[k,0:3].copy()
            quad[0] = quad[0]/1e-15
            quad[1] = quad[1]/quad[1].max()
            quad[2] = quad[2]/quad[2].max()
            self.PlotPost1D(quad[0],quad[1:3],xlabel='time',ylabel='$\Psi_{1,2}$',xunits='$fs$',labels=['$\Psi_1$','$\Psi_2$'],save_name='q_sampled_'+mode+'.png')
            del quad
            quad = self.S_sampled[k,0:2].copy()
            quad[0] = quad[0]/1e-15
            self.PlotPost1D(quad[0],quad[1],xlabel='time',ylabel='$S$',xunits='$fs$',labels=['$S$'],save_name='S_sampled_'+mode+'.png')
            del quad
            quad = self.I_sampled[k,0:2].copy()
            quad[0] = quad[0]/1e-15
            self.PlotPost1D(quad[0],quad[1],xlabel='time',ylabel='$S$',xunits='$fs$',labels=['$I$'],save_name='I_sampled_'+mode+'.png')
            del quad
            quad = self.S_sampled[k,0:3].copy()
            quad[0] = quad[0]
            quad[2] = quad[2]
            self.PlotPost1D(quad[0],quad[2],xlabel='time',ylabel='$x$',xunits='$fs$',yunits='$\mu m$',labels=['$S$'],save_name='S_position_'+mode+'.png',lp=True)
            del quad
            quad = self.I_sampled[k,0:3].copy()
            quad[0] = quad[0]
            quad[2] = quad[2]
            self.PlotPost1D(quad[0],quad[2],xlabel='time',ylabel='$x$',xunits='$fs$',yunits='$\mu m$',labels=['$i$'],save_name='I_position_'+mode+'.png',lp=True)

    def PlotPost1D(self,xdata,ydata,xlabel='$x$',ylabel='$y$',xunits='$a.u.$',yunits='$a.u.$',labels=['plot_1'],save_name='plot.png',lp=False):
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        fig_path = self.outdir

        from matplotlib import pylab as plt

        plt.close('all')
        try:
            l,w = ydata.shape
        except ValueError:
            l = w = 1

        plt.figure()

        if lp:
            plt.plot(xdata/1e-15,(self.cmin*xdata)/1e-6,label='$\mathscr{L}_{min}$')
            plt.plot(xdata/1e-15,(self.cmax*xdata)/1e-6,label='$\mathscr{L}_{max}$')
            plt.plot(xdata/1e-15,(self.co*self.vrip*xdata)/1e-6,label='$\mathscr{R}$')
            plt.plot(xdata/1e-15,ydata/1e-6,label=labels[0])
        else:
            if l>1:
                for dim in range(0,l):
                    plt.plot(xdata,ydata[dim],label=labels[dim])
            else:
                plt.plot(xdata,ydata,label=labels[0])

        xlabel = xlabel+'('+xunits+')'
        ylabel = ylabel+'('+yunits+')'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(frameon=False)
        plt.draw()
        plt.savefig(os.path.join(fig_path,save_name),dpi=240)
        plt.close()

    def PeakSampleFrame(self,frame,n):
        solution = IO()
        solution.path = self.path
        solution.file_prefix = self.file_prefix
        solution.frame = frame
        solution.read_petsc()
        x = solution.x.centers.copy()
        # print x.shape
        for k,mode in enumerate(self.sample_mode):
            if self.num_dim==1:
                self.FieldIntensity1D(solution)
                self.Poyinting1D(solution)
                self.q_sampled[k,0,n] = solution.t
                self.S_sampled[k,0,n] = solution.t
                self.I_sampled[k,0,n] = solution.t
                maxif = 0.0
                if mode=='peak':
                    self.q_sampled[k,1,n],self.q_sampled[k,3,n] = self.PeakSample1D(solution.q[0])
                    if self.q_sampled[k,1,n]==-99:
                        self.q_sampled[k,1,n] = self.q_sampled[k,1,n-1]
                        self.q_sampled[k,3,n] = self.q_sampled[k,3,n-1]
                    self.q_sampled[k,2,n],self.q_sampled[k,4,n] = self.PeakSample1D(solution.q[1])
                    if self.q_sampled[k,2,n]==-99:
                        self.q_sampled[k,2,n] = self.q_sampled[k,2,n-1]
                        self.q_sampled[k,4,n] = self.q_sampled[k,4,n-1]
                    self.I_sampled[k,1,n],self.I_sampled[k,2,n] = self.PeakSample1D(self.I)
                    if self.q_sampled[k,1,n]==-99:
                        self.I_sampled[k,1,n] = self.I_sampled[k,1,n-1]
                        self.I_sampled[k,2,n] = self.I_sampled[k,2,n-1]
                    self.S_sampled[k,1,n],self.S_sampled[k,2,n] = self.PeakSample1D(self.S)
                    if self.q_sampled[k,1,n]==-99:
                        self.S_sampled[k,1,n] = self.S_sampled[k,1,n-1]
                        self.S_sampled[k,2,n] = self.S_sampled[k,2,n-1]
                if mode=='width':
                    self.q_sampled[k,1,n],self.q_sampled[k,3,n] = self.PeakWidth1D(solution.q[0],x)
                    if self.q_sampled[k,1,n]==-99:
                        self.q_sampled[k,1,n] = self.q_sampled[k,1,n-1]
                        self.q_sampled[k,3,n] = self.q_sampled[k,3,n-1]
                    self.q_sampled[k,2,n],self.q_sampled[k,4,n] = self.PeakWidth1D(solution.q[1],x)
                    if self.q_sampled[k,2,n]==-99:
                        self.q_sampled[k,2,n] = self.q_sampled[k,2,n-1]
                        self.q_sampled[k,4,n] = self.q_sampled[k,4,n-1]
                    self.I_sampled[k,1,n],self.I_sampled[k,2,n] = self.PeakWidth1D(self.I,x)
                    if self.q_sampled[k,1,n]==-99:
                        self.I_sampled[k,1,n] = self.I_sampled[k,1,n-1]
                        self.I_sampled[k,2,n] = self.I_sampled[k,2,n-1]
                    self.S_sampled[k,1,n],self.S_sampled[k,2,n] = self.PeakWidth1D(self.S,x)
                    if self.q_sampled[k,1,n]==-99:
                        self.S_sampled[k,1,n] = self.S_sampled[k,1,n-1]
                        self.S_sampled[k,2,n] = self.S_sampled[k,2,n-1]

    def PeakWidth1D(self,field,x):
        # print x.size, x.shape
        # print x[np.argwhere(field>=field.max()/2.0).flatten()]
        try:
            xi = x[np.argwhere(field>=field.max()/2.0)].flatten().min()        
            xf = x[np.argwhere(field>=field.max()/2.0)].flatten().max()
            location = np.argwhere(field>=field.max()/2.0).flatten()
            width = np.abs(xf-xi)
        except ValueError:
            xi = -99
            xf = -99
            location = -99
            width = -99


        return width, xi

    def PeakSample1D(self,field):
        try:
            if field.max().size==1:
                maxfield = field.max()
            else:
                maxfield = np.max(field.max())

            location = (field==maxfield).argmax()
        except ValueError:
            xi = -99
            xf = -99
            location = -99
            maxfield = -99

        return maxfield, location


    def __init__(self,path='./',file_prefix='claw',get_all_frames=True,frame=0,plot=False,save_calc=True,plot_calc=True,outdir='_postcalc',clean_outdir=False):
        
        self.__setattr__('plot_frames', np.zeros([1]))
        self.__setattr__('path',path)
        self.__setattr__('file_prefix',file_prefix)
        self.__setattr__('outdir',os.path.join(path,outdir))
        self.__setattr__('get_all_frames',get_all_frames)
        self.__setattr__('save_calc',save_calc)
        self.__setattr__('plot_calc',plot_calc)
        self.__setattr__('_cleanold',clean_outdir)
        self.__setattr__('sample_frames', np.zeros([0]))
        self.__setattr__('sample_all', True)
        self.__setattr__('sample_mode',['peak','width'])
        self.__setattr__('plot_sample',True)
        self.__setattr__('co',299792458.0)
        self.__setattr__('vrip',0.6)
        if self._cleanold:
            self._DeleteOldFolder()
        self.get_num_cells()
        self.get_num_dim()

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