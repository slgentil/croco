import os, sys
from os.path import join
import shutil
from glob import glob
import yaml
import time

class run(object):
    """ Object to automate run launchings
    """
    
    def __init__(self, rundir, nbchains=1, elaptim=None, config=None,
                 jobname='job', workdir=None, restart=False,
                 launch=False, qsub=True,
                 **params):
        #
        self.startdir = os.getcwd()
        self.user = os.getenv('USER')
        if workdir is None:
            self.workdir = os.getenv('SCRATCH')
        elif not os.path.isdir(workdir):
            self.workdir = os.getenv(workdir)
        else:
            self.workdir = workdir
        # main directory where the simulation will be run:
        self.rpath = join(self.workdir, rundir)
        print('Run will be stored in {}'.format(self.rpath))
        self._create_rpath()
        # individual run directories
        self.nbchains = nbchains
        self._create_tdirs()
        # backup code files
        self._backup_code()
        # change input parameters in croco.in
        self.restart = restart
        self.update_input_file(**params)    
        # guess config if necessary and create job files
        self.jobname = jobname
        self._load_config(config)
        self._create_job_files()
        # create commands to be executed
        self.qsub = qsub
        self._create_commands()
        # launch runs
        if launch:
            self.launch()
        os.chdir(self.startdir)

    def _create_rpath(self):
        if os.path.exists(self.rpath) :    
            os.system('rm -Rf '+self.rpath)
        os.mkdir(self.rpath)
        #os.chdir(self.rpath)
    
    def _create_tdirs(self):
        #
        os.mkdir(join(self.rpath, 't0'))
        #
        self.tdirs = [join(self.rpath, 't'+str(t)) 
                        for t in range(1, self.nbchains+1)]
        for t, tdir in zip(range(1,self.nbchains+1), self.tdirs):
            os.mkdir(tdir)
            if t != 0 :
                # copy executables and input files
                shutil.copy(self.startdir+'/croco',tdir)
                shutil.copy(self.startdir+'/croco.in',tdir)
                shutil.copy(self.startdir+'/config.yaml',tdir)
                shutil.copy(self.startdir+'/iodef.xml',tdir)
                shutil.copy(self.startdir+'/domain_def.xml',tdir)
                shutil.copy(self.startdir+'/axis_def.xml',tdir)
                shutil.copy(self.startdir+'/grid_def.xml',tdir)
                shutil.copy(self.startdir+'/field_def.xml',tdir)
                shutil.copy(self.startdir+'/floats.in',tdir)
                shutil.copy(self.startdir+'/xios_server.exe',tdir)
                
    def _backup_code(self):
        _bpath = join(self.rpath,'backup/')
        if ~os.path.exists(_bpath):
            os.mkdir(_bpath)
            code_files = glob(join(self.startdir,'*.[Fh]'))
            #cmd='find '+startdir+' -name "*.[Fh]" -exec grep -l aponte {} \;'
            #listfiles = os.popen(cmd).readlines()  
            for f in code_files:
                shutil.copy(f,_bpath)

    def update_input_file(self, **params):
        """ Update parameters in croco.in
        
        Parameters
        ----------
            **param: dict
                Dictionary of parameters to change, values are 
                given as lists or single values
        """

        for t, tdir in zip(range(1,self.nbchains+1), self.tdirs):
            def wrap(p):
                if isinstance(p, list):
                    return '         '+'  '.join(map(str,p))+'\n'
                else:
                    return '         '+str(p)+'\n'
            # move to tdir
            os.chdir(tdir)
            with open('croco.in','r') as f:
                lines = f.readlines()
                for index, line in enumerate(lines):
                    key = line.split(':')[0]
                    if key=='initial' and t==1:
                        if not self.restart:
                            lines[index+1]=wrap(0)
                        else:
                            lines[index+1]=wrap(1)
                    if key=='initial' and t>1:
                        lines[index+1]=wrap(1)
                    if key=='nrpfflt' and t>1:
                        lines[index+1]=lines[index+1].rstrip()[:-2]+'1\n'
                    for p, v in params.items():
                        if key==p:
                            lines[index+1]=wrap(v)
            with open('croco.in','w') as f:
                f.writelines(lines)

    def _load_config(self, config):
        with open(join(self.startdir,'config.yaml')) as f:
            configs = yaml.full_load(f)
        if config is None:
            self.config = list(configs.keys())[0]
            _required = ['elapse_time', 'nbproc_roms', 'nbproc_xios']
            assert all([r in configs[self.config] for r in _required])
            for k, v in configs[self.config].items():
                setattr(self, k, v)
        self.nb_cores = self.nbproc_roms+self.nbproc_xios
        self.nb_nodes = int((self.nb_cores)/28)+1

    def _create_job_files(self):
        
        for t, tdir in zip(range(1,self.nbchains+1), self.tdirs):
            with open(join(tdir,'job_datarmor'),'w') as f:
                f.write('#!/bin/csh\n')
                f.write('#PBS -N '+self.jobname+str(t)+'\n')
                f.write('#PBS -q mpi\n')
                f.write('#PBS -l select='+str(self.nb_nodes)+
                        ':ncpus=28:mpiprocs=28:mem=120G\n')
                f.write('#PBS -l walltime='+self.elapse_time+'\n')
                f.write('\n')
                f.write('# get the path fr command module,mpirun\n')
                f.write('source /usr/share/Modules/3.2.10/init/csh\n') 
                f.write('module purge\n')  
                f.write('module load   NETCDF/4.3.3.1-mpt-intel2016'+
                        '  intel-cmkl-16/16.0.4.258\n')  
                f.write('\n')
                f.write('# cd to the directory you submitted your job\n')
                f.write('cd $PBS_O_WORKDIR\n')
                f.write('\n')
                f.write('date\n') 
                if t!=1 or self.restart:
                   f.write('cp ../t'+str(t-1)+'/jetn_rst*.nc .\n')  
                   f.write('cp ../t'+str(t-1)+'/float.rst.* .\n')   
                f.write('time $MPI_LAUNCH -n '+str(self.nbproc_roms)
                        +' croco : -n '+str(self.nbproc_xios)
                        +' xios_server.exe >& output.mpi\n')
                    
    def _create_commands(self):

        # all runs
        for t, tdir in zip(range(1,self.nbchains+1), self.tdirs):
            os.chdir(tdir)
            if t == 1:
                command='qsub -h job_datarmor'  
            else:
                f=os.popen('qselect -N '+self.jobname+str(t-1)+' -u '
                           +self.user)
                numjob=f.read()
                command='qsub -W depend=afterany:'+str(numjob[:-1]) \
                        +' job_datarmor'
            if self.qsub:
                os.system(command)
        
        # first run
        os.chdir(join(self.rpath,'t1'))
        f=os.popen('qselect -N '+self.jobname+str(1)+' -u '
                    +self.user)
        numjob=f.read()
        self.command = 'qrls '+numjob[:-1]
        #print(self.command)

        if self.restart and isinstance(self.restart,str):
            for f in glob(self.restart):
                shutil.copy(f, join(self.rpath,'t0'))        
        elif self.restart:
            print('Put the restart files in '+join(self.restart,'t0'))
                
    def launch(self):
        """ Launch simulations
        """
        time.sleep(1)
        os.chdir(join(self.rpath,'t1'))
        os.system(self.command)
        os.chdir(self.startdir)

