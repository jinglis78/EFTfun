import os,multiprocessing,subprocess,time,stat
import logging

def prepareScripts(outname,options,pois,njobs,nsteps):
    if len(nsteps)==1:
        total_nsteps=nsteps*len(pois)
    assert len(nsteps)==len(pois)
    from operator import mul
    total_nsteps=reduce(mul,nsteps)
    nsteps_job=total_nsteps//njobs
    remainder=total_nsteps%njobs
    first=0
    if not '/' in outname:
        outname='/'.join(['results',outname])
    if not os.path.exists(outname):
        os.makedirs(outname)
    jobs=[]
    for ijob in range(njobs):
        jobname=outname+'/'+"job_{}".format(ijob)
        last=first+nsteps_job
        if remainder>0:
            last+=1
            remainder-=1
        options_dict = vars(options)            
        new_options={'method':'scan','first':first,'last':last,'njobs':1,'outname':jobname}
        f=open(jobname+'.sh','w')
        f.write('#!/usr/bin/env bash\n')
        f.write('cd '+os.getcwd()+'\n')
        f.write('./bin/eftfun.py')
        for o in options_dict:
            if not o in new_options:
                if options_dict[o] is None:
                    continue
                if isinstance(options_dict[o],bool):
                    if options_dict[o]:
                        f.write(' --{}'.format(o))
                else:
                    f.write(' --{} {}'.format(o,options_dict[o]))
        for o in new_options:
            f.write(' --{} {}'.format(o,new_options[o]))
        f.write(' --parallel_job')
        f.write(' --quiet')
        f.write('\n')
        f.write('touch {}.done\n'.format(jobname))
        chmodX(jobname+'.sh')
        jobs.append(jobname)
        first=last            
        f.close()
    os.makedirs(outname+'/combined')
    f=open(outname+'/plot.sh','w')
    f.write('#!/usr/bin/env bash\n')
    f.write('cd '+os.getcwd()+'\n')
    f.write('./bin/eftfun.py -i "{}" -m collect -p {}  -o {}'.format(' '.join([j+'.root' for j in jobs]),','.join(pois),outname+'/combined'))
    f.close()
    chmodX(outname+'/plot.sh')
    return jobs

def chmodX(fn):
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)

def runParallel(parallel_jobs,cluster=None):
    outfiles=[j.replace('.sh','.root') for j in parallel_jobs]
    donefiles=[j.replace('.sh','.done') for j in parallel_jobs]
    for l in donefiles:
        if os.path.exists(l):
            raise Exception('This job was already run (maybe try a clean directory?)',l.replace('.done',''))

    submitJobs(parallel_jobs,donefiles,cluster)
    for d in donefiles:
        if not os.path.exists(d):
            raise Exception('This job did not run?',d.replace('.done',''))
    for o in outfiles:
        if not os.path.exists(o):
            raise Exception('Missing outfile',o)
    return outfiles

def submitJobs(scripts,outputs,cluster):
    ncores=multiprocessing.cpu_count()
    scripts_running=[]
    expected_outputs=[]
    scripts_to_submit=list(scripts)
    outputs_to_expect=list(outputs)
    while len(scripts_to_submit)>0 or len(scripts_running)>0 :
        for script,output in zip(scripts_running,expected_outputs):
            if os.path.exists(output):
                scripts_running.remove(script)
                expected_outputs.remove(output)
                logging.info('finished '+script)

        nsubmit=ncores-len(scripts_running)
        nsubmit=min(len(scripts_to_submit),nsubmit)
        for s,o in zip(scripts_to_submit[:nsubmit],outputs_to_expect[:nsubmit]):
            scripts_to_submit.remove(s)
            scripts_running.append(s)
            expected_outputs.append(o)
            outputs_to_expect.remove(o)
            if cluster is None:
                runLocal(s)
            elif cluster=='condor':
                runCondor(s)
            else:
                raise Exception('Cluster type '+c+' not implemented')
        time.sleep(3)

def runLocal(script):
    logging.info('running '+script)
    run = subprocess.Popen([script],stdout=None,stderr=None,stdin=None)        

def runConodor(jobs):
    raise Exception('Not yet implemented')
