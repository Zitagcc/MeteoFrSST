#!/usr/bin/env python

from os.path import abspath,dirname,basename,isfile,isdir
from os import mkdir,remove,system,chmod
from time import sleep
from datetime import datetime
from socket import gethostname

def gen_platform(exe_name, topo, ncpu, cores_per_node, bandwidth, latency, params, amm="amm1", model=dict(congestion='pisces', node='simple', memory='simple', nic='pisces', switch='pisces'), allocation="first_available", indexing="block", routing='minimal', stat=0, app_argv=''):
    platform = './%08d/parameters.ini' % ncpu

    if topo == 'fattree':
        routing = 'fattree'

    if stat == 0:
        nic_stat = ''
        xbar_stat = ''
        buff_stat = ''
    else:
        nic_stat = """
    traffic_matrix {
      fileroot = %08d/pattern
    }
    packet_allocator = delay_stats
    ejection {
      stats = multi
      callbacks = congestion_spyplot delay_histogram
      congestion_spyplot {
        fileroot =  %08d/congestion
      }
      delay_histogram {
        fileroot=  %08d/delay
        bin_size = %fns
      }
    }
""" % (ncpu, ncpu, ncpu, latency['nic'])
        xbar_stat = '    stats = congestion_delay'
        buff_stat = """
  output_buffer {
    stats = congestion_delay
  }
"""

    with open(platform, 'wt') as f:
        f.write("""
# Network parameters
amm_model = %s
congestion_model = %s


# Node parameters
node {
  app1 {
    # Launch parameters
    allocation = %s
    indexing = %s
    launch_cmd = aprun -n %d -N %d
    name = %s
    argv = %s
  }

  model = %s

  ncores = %d
  proc {
    frequency = 5.0GHz
  }

  memory {
    model = %s
    bandwidth = %fGB/s
    latency = %fns
  }

  nic {
    model = %s
    injection {
      bandwidth = %fGB/s
      latency = %fns
    }
%s
  }
}

# Switch
switch {
  model = %s

  xbar {
    bandwidth = %fGB/s
    latency = %fns
%s
  }
%s
  ejection {
    bandwidth = %fGB/s
    latency = %fns
  }

  link {
    bandwidth = %fGB/s
    latency = %fns
  }

  router {
    name = %s
  }
}
""" % (amm, model['congestion'], allocation, indexing, ncpu, cores_per_node, exe_name.split('.')[0], app_argv, model['node'], cores_per_node, model['memory'], bandwidth['memory'], latency['memory'], model['nic'], bandwidth['nic'], latency['nic'], nic_stat, model['switch'], bandwidth['switch_xbar'], latency['switch_xbar'], xbar_stat, buff_stat, bandwidth['switch_ejection'], latency['switch_ejection'], bandwidth['switch_link'], latency['switch_link'], routing))

        words = params.split(':')
        if topo == 'torus' or topo == "hypercube" or topo == "fattree":
            f.write("""
# Topology
topology {
  name = %s
  geometry = %s
  concentration = %s
}
""" % (topo, ' '.join(words[:-1]), words[-1]))
        elif topo == 'dragonfly':
            f.write("""
# Topology
topology {
  name = %s
  geometry = %s
  concentration = %s
  group_connections = %s
}
""" % (topo, ' '.join(words[:-2]), words[-2], words[-1]))

    return platform

def gen_hostfile(pfname, ncpu, cores_per_node):
    nnode = int(ncpu / cores_per_node)
    if ncpu % cores_per_node > 0:
        nnode += 1

    hostfile = './%08d/hostfile' % ncpu

    with open(hostfile, 'wt') as f:
        for n in range(1, nnode+1):
            f.write('%s%d:%d\n' % (pfname, n, cores_per_node))

    return hostfile

# partition
def tile_setup_2d(ncpu, nx, ny,):
    # aspect ratio of the grid
    aspect_grid = float(nx) / ny
    # find the prime factors of ncpu
    factors = prime_factors(ncpu)
    nfactors = len(factors)

    aspect_tile = 1.0E20
    aspect_tile_opt = 1.0E20
    ntile_x_opt = 0
    ntile_y_opt = 0;
    # partition so that aspect ratio of tile is closed to that of grid as far as possible
    for k in range(1, nfactors+1):
        for j in range(nfactors+1-k):
            mx = 1
            for i in range(j, j+k):
                mx *= factors[i]
            my = ncpu / mx
            aspect = float(mx) / my
            # best closed to the aspect ratio of grid but mx <= nx and my <= ny is not guarantteed
            if abs(aspect - aspect_grid) < abs(aspect_tile - aspect_grid):
                ntile_x = mx
                ntile_y = my
                aspect_tile = aspect
            # closed to the aspect ratio of grid also mx <= nx and my <= ny is satified
            if abs(aspect - aspect_grid) < abs(aspect_tile_opt - aspect_grid) and  mx <= nx and my <= ny:
                ntile_x_opt = mx
                ntile_y_opt = my
                aspect_tile_opt = aspect
    # use optimal partition if found one
    if ntile_x_opt > 0 and ntile_y_opt > 0:
        return (ntile_x_opt, ntile_y_opt)
    else:
        return (ntile_x, ntile_y)

# find all prime factors
def prime_factors(num):
    factors = []
    n = num

    # first check whether 2 is its prime factor
    while n % 2 == 0:
        factors.append(2)
        n /= 2

    # find all odd prime factors
    j = 3
    while n > 1:
        while n % j == 0:
            factors.append(j)
            n /= j
        j += 2

    # if num is a prime number
    if len(factors) == 0:
        factors.append(num)

    return factors;

# find a number whose prime factors are not greater than max_prime
def find_num(num, max_prime):
    for i in range(0,1001,2):
        # after
        fmax = 0
        for f in prime_factors(num+i):
            if f > fmax:
                fmax = f
        if fmax < max_prime:
            return num+i
        # before
        fmax = 0
        for f in prime_factors(num-i):
            if f > fmax:
                fmax = f
        if fmax < max_prime:
            return num-i
    print('Can not find a number between %d and %d' % (num-1000, num+1000))
    return num

def fprint_out(line):
    from sys import stdout

    stdout.write(line+'\n')
    stdout.flush()

def submit_mpi(work_dir,exe_file,cores_per_node,ntasks,wait=None,max_wait=1440,command_args='',dependency=None,sid=''):
    if dependency is None:
        dependency = ''
    else:
        dependency = """
until [[ -e '%s' ]]
do
    sleep %d
done
""" % (dependency['chk_file'], dependency['chk_freq'])

    exe_name = basename(exe_file) + sid
    cmd_file = '%s/%08d/%s.cmd' % (work_dir, ntasks, exe_name)
    chk_file = '%s/%08d/%s.chk' % (work_dir, ntasks, exe_name)
    job_file = '%s/%08d/%s.job' % (work_dir, ntasks, exe_name)

    if ntasks % cores_per_node > 0:
        nodes = int(ntasks / cores_per_node) + 1
    else:
        nodes = ntasks / cores_per_node

    # prepare the command file for the job
    with open(cmd_file, 'wt') as f:
        f.write("""#!/bin/bash
#SBATCH --job-name="%s"
#SBATCH --nodes=%d
#SBATCH --partition=normal64,huge256
#SBATCH --exclusive
#SBATCH --output=%08d/slurm-%%j.out

if [ -f ${HOME}/.bashrc ]; then
    ${HOME}/.bashrc
fi

%s

# put the following exports in profile and use intelmpi/2017.1.132 if "UCM connect: REQ RETRIES EXHAUSTED"
export I_MPI_FABRICS=shm:dapl
export I_MPI_FALLBACK=disable
export I_MPI_SLURM_EXT=1
export I_MPI_LARGE_SCALE_THRESHOLD=%d # a value larger than MPI tasks
export I_MPI_DYNAMIC_CONNECTION=1
export I_MPI_CHECK_DAPL_PROVIDER_COMPATIBILITY=0
export DAPL_UCM_REP_TIME=8000
export DAPL_UCM_RTU_TIME=8000
export DAPL_UCM_CQ_SIZE=2000
export DAPL_UCM_QP_SIZE=2000
export DAPL_UCM_DREQ_RETRY=4
export DAPL_UCM_WAIT_TIME=10000

mpirun -np %d %s/%s %s
wait
touch %s
""" % (exe_name, nodes, ntasks, dependency, ntasks, ntasks, work_dir, exe_file, command_args, chk_file))

    # submit the job
    fprint_out('submit %s (%d nodes and %d tasks) at %s from %s' % (cmd_file, nodes, ntasks, datetime.now().strftime('%Y-%m-%d %H:%M:%S'), gethostname()))
    if isfile(chk_file):
        remove(chk_file)
    system('sbatch --no-requeue %s > %s' % (cmd_file, job_file))

    if wait is not None:
        # waiting for completion
        t = 0
        while t < max_wait and not isfile(chk_file):
            sleep(wait)
            t = t + wait / 60.0

        # task finished or max wait time is reached
        if isfile(chk_file):
            return chk_file
        else:
            with open(job_file, 'rt') as f:
                line = f.readline()
                words = line.split()
                system("scancel %s" % words[-1])
            return None
    else:
        return chk_file

def submit_smpi(work_dir,exe_file,ntasks,nnode=1,cores_per_node=1,wait=None,max_wait=1440,command_args='',dependency=None,sid='',app_argv=None):
    if dependency is None:
        dependency = ''
    else:
        dependency = """
until [[ -e '%s' ]]
do
    sleep %d
done
""" % (dependency['chk_file'], dependency['chk_freq'])

    if len(ntasks) % cores_per_node > 0:
        nodes = int(len(ntasks) / cores_per_node) + 1
    else:
        nodes = len(ntasks) / cores_per_node
    if len(ntasks) == 1:
        nodes = nnode

    exe_name = basename(exe_file) + sid

    for ntask in ntasks:
        cmd_file = '%s/%08d/%s.cmd' % (work_dir, ntask, exe_name)
        with open(cmd_file, 'wt') as fcmd:
            fcmd.write("#!/bin/sh\n")
            fcmd.write("%s/%s -f %s/%08d/parameters.ini %s\n" % (work_dir, exe_file, work_dir, ntask, command_args))
        chmod(cmd_file, 448) # 0700 = 448: -rwx------

    bat_file = '%s/%08d/%s.bat' % (work_dir, ntasks[0], exe_name)
    chk_file = '%s/%08d/%s.chk' % (work_dir, ntasks[0], exe_name)
    job_file = '%s/%08d/%s.job' % (work_dir, ntasks[0], exe_name)
    mac_file = '%s/%08d/%s.mac' % (work_dir, ntasks[0], exe_name)

    if len(ntasks) > 4 or ntasks[-1] < 70000:
        reservation = ''
        partition = 'normal64'
    else:
        reservation = '#SBATCH --reservation=test'
        reservation = ''
        partition = 'huge256'

    with open(bat_file, 'wt') as fbat:
        fbat.write("""#!/bin/bash
#SBATCH --job-name="%s"
#SBATCH --nodes=%d
#SBATCH --ntasks-per-node=%d
#SBATCH --partition=%s
%s
#SBATCH --exclusive
#SBATCH --output=%08d/slurm-%%j.out

if [ -f ${HOME}/.bashrc ]; then
    ${HOME}/.bashrc
fi

%s

> %s
nn=0
for h in `scontrol show hostnames $SLURM_NODELIST`
do
  echo "$h:1" >> %s
  nn=$((nn+1))
done

""" % (exe_name, nodes, cores_per_node, partition, reservation, ntasks[0], dependency, mac_file, mac_file))
        if len(ntasks) == 1:
            ntask = ntasks[0];
            if app_argv is None or app_argv[-1] > 0:
                fbat.write("mpirun -machinefile %s -np %d %s/%s -f %s/%08d/parameters.ini %s\n" % (mac_file, nnode*cores_per_node, work_dir, exe_file, work_dir, ntask, command_args))
            else:
                radixes = []
                for k in range(2,abs(app_argv[-1])):
                    por = k**int(log(ntask)/log(k))
                    rem = ntask - por
                    if rem <= por:
                        radixes.append(k)
                for k in radixes:
                    fbat.write("sed -i 's/argv = .*/argv = %d %d %d/g' %s/%08d/parameters.ini\n" % (app_argv[0], app_argv[1], k, work_dir, ntask))
                    fbat.write("mpirun -machinefile %s -np %d %s/%s -f %s/%08d/parameters.ini %s > %s/%08d/%03d.log 2>&1\n" % (mac_file, nnode*cores_per_node, work_dir, exe_file, work_dir, ntask, command_args, work_dir, ntask, k))
        else:
            for ntask in ntasks:
                fbat.write('%s/%08d/%s.cmd > %s/%08d/%s.log 2>&1 &\n' % (work_dir, ntask, exe_name, work_dir, ntask, exe_name))
            fbat.write('wait\n')
        fbat.write('touch %s\n' % (chk_file))
    # submit the job
    fprint_out('submit %s (%d nodes and %d tasks) at %s from %s' % (bat_file, nodes, len(ntasks), datetime.now().strftime('%Y-%m-%d %H:%M:%S'), gethostname()))
    if isfile(chk_file):
        remove(chk_file)
    system('sbatch --no-requeue %s > %s' % (bat_file, job_file))

    if wait is not None:
        # waiting for completion
        t = 0
        while t < max_wait and not isfile(chk_file):
            sleep(wait)
            t = t + wait / 60.0

        # task finished or max wait time is reached
        if isfile(chk_file):
            return chk_file
        else:
            with open(job_file, 'rt') as f:
                line = f.readline()
                words = line.split()
                system("scancel %s" % words[-1])
            return None
    else:
        return chk_file

def run_cmd(work_dir,exe_file,ntasks,wait=None,max_wait=1440,command_args='',dependency=None,sid=''):
    exe_name = basename(exe_file) + sid
    cmd_file = '%s/%08d/%s.cmd' % (work_dir, ntasks[0], exe_name)
    chk_file = '%s/%08d/%s.chk' % (work_dir, ntasks[0], exe_name)
    job_file = '%s/%08d/%s.job' % (work_dir, ntasks[0], exe_name)

    # prepare the command file for the job
    with open(cmd_file, 'wt') as f:
        f.write("#!/bin/bash\n\n")
        f.write("export GO_MBLOCK_SHMID=`%s/shm_mgr.x %08d %d`\n" % (work_dir, ntasks[0], ntasks[0]))

        for ntask in ntasks:
            f.write("/usr/bin/time --format='Elapsed time: %%e seconds.' %s/%s -s $GO_MBLOCK_SHMID %s\n" % (work_dir, exe_file, command_args))

        f.write("touch %s\n" % (chk_file))
        f.write("%s/shm_mgr.x $GO_MBLOCK_SHMID\n" % (work_dir))

    # submit the job
    fprint_out('run %s at %s from %s' % (cmd_file, datetime.now().strftime('%Y-%m-%d %H:%M:%S'), gethostname()))
    if isfile(chk_file):
        remove(chk_file)
    system('bash %s >& %s' % (cmd_file, job_file))

    if wait is not None:
        # waiting for completion
        t = 0
        while t < max_wait and not isfile(chk_file):
            sleep(wait)
            t = t + wait / 60.0
    else:
        return chk_file

def wait_all(chk_files, wait, max_wait=1440):
    # waiting for completion
    t = 0
    while t < max_wait:
        n = 0
        for chk_file in chk_files:
            if isfile(chk_file):
                n += 1

        if n == len(chk_files):
            break
        else:
            sleep(wait)
            t = t + wait / 60.0

def run_cpus(work_dir, exe_file, nx, ny, nz, upper, nstep, njob, wait, restart, allocation='first_available', indexing='block', topo=None, bandwidth=None, latency=None, params=None, nnode=1, cores_per_node=1, routing='minimal', stat=0):
    if (njob < 0):
        njobs = [1,2,8,16]
        ntasks = []
    else:
        chk_files = []

    MINNCPU = 100
    min_ncpu = log10(MINNCPU);
    max_ncpu = log10(upper)
    del_ncpu = (max_ncpu - min_ncpu) / nstep
    old_ncpu = 0
    mjob = 0

    dcpu = MINNCPU
    mcpu = 0
    for k in range(nstep + 1):
        #ncpu = int(pow(10, min_ncpu + del_ncpu * k))
        mcpu += dcpu
        if (mcpu / dcpu == 10):
            dcpu *= 10
        if (mcpu > upper):
            break

        ncpu = mcpu
        if (ncpu < MINNCPU):
            ncpu = MINNCPU
        if (ncpu % 2 > 0):
            ncpu += 1
        if (ncpu > upper):
            ncpu = upper

        ncpu = find_num(ncpu, 100)
        if (ncpu > old_ncpu):
            old_ncpu = ncpu
            mjob += 1

            if (ncpu >= restart):
                # create directory to store the output
                odir = '%08d' % ncpu
                if not isdir(odir):
                    mkdir(odir,493)

                if njob == 0:
                    chk_file = submit_mpi(work_dir, exe_file, 40, ncpu, wait=wait, max_wait=180)
                elif njob > 0:
                    # generate the platform description file and hostfile
                    gen_platform(exe_file, topo, ncpu, 1, bandwidth, latency, params, routing=routing, stat=stat, app_argv='%d %d %d' % (nx,ny,nz))
                    # submit
                    chk_file = submit_smpi(work_dir, exe_file, [ncpu], max_wait=7200, nnode=nnode, cores_per_node=cores_per_node, app_argv=(nx,ny,nz))
                    # prevent to submit too quickly
                    sleep(8)
                    # wait
                    chk_files.append(chk_file)
                    if len(chk_files) == njob or ncpu == upper or k == nstep:
                        wait_all(chk_files, wait, max_wait=7200)
                        chk_files = []
                elif njob < 0:
                    ntasks.append(ncpu)
                    # generate the platform description file and hostfile
                    gen_platform(exe_file, topo, ncpu, 1, bandwidth, latency, params, routing=routing, stat=stat)
                    if mjob == njobs[-1] or ncpu == upper or k == nstep:
                        # submit
                        chk_file = submit_smpi(work_dir, exe_file, [ntask for ntask in ntasks if ntask >= restart], wait=wait, max_wait=7200)
                        # next round
                        mjob = 0
                        ntasks = []
                        if len(njobs) > 0:
                            njobs.pop()
    print('****** SUCCEED ******')

if __name__ == '__main__':
    from math import log, log10, pow
    from sys import argv, exit
    from glob import glob

    if len(argv) < 4:
        print("%s nx:ny:ncpu_upper nstep wait [restart | njob allocation:indexing topo:routing nnode:cores_per_node bandwidth latency params [restart]]" % (argv[0]))
        exit()
    else:
        work_dir=dirname(abspath(argv[0]))
        nx, ny, nz, upper = map(lambda x: int(x), argv[1].split(':'))
        nstep = int(argv[2])
        wait = int(argv[3])
        if len(argv) > 10:
            njob, stat = map(lambda x: int(x), argv[4].split(':'))
            allocation, indexing = argv[5].split(':')
            topo, routing = argv[6].split(':')
            nnode, cores_per_node = map(lambda x: int(x), argv[7].split(':'))
            A = list(map(lambda x: float(x), argv[8].split(':')))
            bandwidth = dict(memory=A[0], nic=A[1], switch_xbar=A[2], switch_ejection=A[3], switch_link=A[4])
            A = list(map(lambda x: float(x), argv[9].split(':')))
            latency = dict(memory=A[0], nic=A[1], switch_xbar=A[2], switch_ejection=A[3], switch_link=A[4])
            params = argv[10]
            # restart
            if len(argv) == 12:
                restart = int(argv[11])
            else:
                restart = 2
        else:
            njob = 0

            # restart
            if len(argv) == 5:
                restart = int(argv[4])
            else:
                restart = 2

    exe_file = glob('./*.exe')
    if len(exe_file) == 0:
        print('Can not find the executable file!')
        exit()
    else:
        exe_file = basename(exe_file[0])

    if len(argv) > 10:
        run_cpus(work_dir, exe_file, nx, ny, nz, upper, nstep, njob, wait, restart, allocation=allocation, indexing=indexing, topo=topo, bandwidth=bandwidth, latency=latency, params=params, nnode = nnode, cores_per_node=cores_per_node, routing=routing, stat=stat)
    else:
        run_cpus(work_dir, exe_file, nx, ny, upper, nstep, njob, wait, restart)
