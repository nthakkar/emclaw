import os
import glob
import shutil


src_dir         = '/sim-tools/em-pyclaw/maxwell_1d_source'
src_subdir      = '_postcalc'
src_dir_prefix  = '_moving_'
dst_dir         = '/Users/sanromd/Dropbox/research/fullwave-time/prlx/figures'
cp_name         = ''
respect_copied  = True
load_all        = True
get_frames      = [6, 35, 70, 98]
rename          = ['a','b','c','d','e','f','g','h','i','j','k']
post_calc       = True
transfer_files  = False

# def copyfigs(src,dst,overwrite=False):
#     if os.path.isfile(dst):
#         if overwrite:
#             print 'overwritting ',src_name,' as ',dst_name
#             shutil.copyfile(src_file,dst_file)
#         else: 
#             print 'file already exists'
#     else:
#         shutil.copyfile(src_file,dst_file)
if load_all:
    src_dir_prefix = src_dir_prefix+'*'
for src_path in glob.glob(os.path.join(src_dir,src_dir_prefix)):
    head, tail  =  os.path.split(src_path)
    tail1,tail2 = tail.split('_moving_')
    # tail1, tail2 =  tail.split(src_dir_prefix)
    dst_base_name = cp_name+tail2
    if post_calc:
        from postprocess import Postprocess
        print src_path
        if not os.path.isfile(os.path.join(src_path,'claw_aux.ptc0000')):
            shutil.copy(os.path.join(src_path,'claw_aux.ptc0001'),os.path.join(src_path,'claw_aux.ptc0000'))
        pst = Postprocess(path=src_path,clean_outdir=True)
        if src_path.find('v0')==-1:
            pst.vrip = int(src_path[src_path.find('vm0')+3:src_path.find('vm0')+5])/100.0
        else:
            pst.vrip = int(src_path[src_path.find('v0')+2:src_path.find('v0')+5])/100.0

        if not src_path.find('1um')==-1:
            pst.plot_frames = [6,35,70,105,140,175,198]
            print pst.plot_frames
            pst.CalcFieldQuad()
            pst.PlotFrames()
            pst.Sampling()
        else:
            pst.plot_frames = get_frames
            print pst.plot_frames
            pst.CalcFieldQuad()
            pst.Sampling()
            pst.PlotFrames()
    if transfer_files:
        fig_path = os.path.join(src_path,src_subdir)
        for figfile in glob.glob(os.path.join(fig_path,'*.png')):
            tr = figfile.split('/')
            if not tr[-1].find('frame')==-1:
                tail1,tail2 = tr[-1].split('frame')
                dst_file = dst_base_name+'_'+tail2
            else:
                dst_file = dst_base_name+'_'+tr[-1]

            dst = os.path.join(dst_dir,dst_file)
            print dst
            shutil.copy(figfile,dst)



# from postprocess import Postprocess
# src_path = '../maxwell_1d_source/_moving_vm061_gpulse_1'
# pst = Postprocess(path=src_path)
