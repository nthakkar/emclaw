import os
import glob
import shutil

# define the src_dir and dst_dir
src_dir         = '/Volumes/sim-tools/em-pyclaw/maxwell_1d_source'
src_subdir      = '_plots'
src_dir_prefix  = '_moving'
dst_dir         = '/Users/sanromd/Dropbox/research/fullwave-time/prlx/figures'
cp_name         = 'moving'
respect_copied  = True
overwrite       = False
load_all        = True

# number of frames to be copied and new name assigment if any
get_frames = [6,35,70,98]
frame_s = ['a','b','c','d']


# ................... main routine .................................................
if load_all:
    src_dir_prefix = src_dir_prefix+'*'

for dir_name in glob.glob(os.path.join(src_dir, src_dir_prefix)):
    if os.path.isdir(os.path.join(dir_name,src_subdir)):
        head, tail =  os.path.split(dir_name)
        tail1, tail2 =  tail.split(src_dir_prefix)
        tail2 = cp_name+tail2
        if respect_copied:
            if not os.path.isfile(os.path.join(dir_name,'.copied')):
                open(os.path.join(dir_name,'.copied'),'a').close()
                for frame in get_frames:
                    src_name = 'frame'+str(frame).zfill(4)+'fig0.png'
                    dst_name = tail2+'_'+frame_s[n]+'.png'
                    print src_name            
                    src_file = os.path.join(dir_name,'_plots',src_name)
                    dst_file = os.path.join(dst_dir,dst_name)
                    if os.path.isfile(src_file):
                        if os.path.isfile(dst_file):
                            if overwrite:
                                print 'overwritting ',src_name,' as ',dst_name
                                shutil.copyfile(src_file,dst_file)
                            else: 
                                print 'file already exists'
                        else:
                            print 'copying ',src_name,' as ',dst_name
                            shutil.copyfile(src_file,dst_file)
        else:
            open(os.path.join(dir_name,'.copied'),'a').close()
            for n in range(0,4):
                src_name = 'frame'+frame_n[n]+'fig0.png'
                dst_name = tail2+'_'+frame_s[n]+'.png'
                print src_name            
                src_file = os.path.join(dir_name,'_plots',src_name)
                dst_file = os.path.join(dst_dir,dst_name)
                if os.path.isfile(src_file):
                    if os.path.isfile(dst_file):
                        if overwrite:
                            print 'overwritting ',src_name,' as ',dst_name
                            shutil.copyfile(src_file,dst_file)
                        else: 
                            print 'file already exists'
                    else:
                        print 'copying ',src_name,' as ',dst_name
                        shutil.copyfile(src_file,dst_file)