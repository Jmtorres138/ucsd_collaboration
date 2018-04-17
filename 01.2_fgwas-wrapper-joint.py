#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run fGWAS across annotations
Usage: python fgwas_wrapper_genome-wide.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import operator
import time
from select import select
from math import ceil
from math import floor
import moniter_rescomp_jobs

# globals
fgwas = "LD_LIBRARY_PATH=/apps/well/gsl/2.2.1-gcc4.9.3/lib /users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
home_dir = "/well/mccarthy/users/jason/projects/ucsd_collaboration/"
in_dir=home_dir+"fgwas_input/"
out_dir = home_dir + "fgwas_output/"
input_file=in_dir+"ukbb_diamante-euro.fgwas.gz" # Optional: Run 01.1 script and use this for abbreviated annotations: in_dir+"fgwas_input_file.renamed.fgwas.gz"
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)
if os.path.isdir(job_dir)==False:
    os.mkdir(job_dir)
if os.path.isdir(log_dir)==False:
    os.mkdir(log_dir)
start_index = 9 #0-based index of column in fgwas input file where annotations start
job_prefix = "bing_" # change this if you have concurrent runs of the fgwas_wrapper (for example, running wrapper seperately across tissues)


def step1():
    '''
    Start index is the first index for an annotation in the file
    Here, the start index is 10 (column 11)
    '''
    fin = gzip.open(input_file,'rb')
    annot_list = fin.readline().strip().split()[start_index:]
    fin.close()
    for annot in annot_list:
        job_file = job_dir+job_prefix+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",  "-dists",
                            annot+":"+home_dir+"dist_model", "-o", out_dir+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            annot, "-o", out_dir+annot]
        command = " ".join(command_list)
        # removed #$ -V from script
        script='''
#$ -N %s%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (job_prefix,annot, log_dir,job_prefix+annot,log_dir,job_prefix+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        out_path = out_dir+annot+".llk"
        if os.path.exists(out_path) == False:
            sp.check_call(call)
        if os.path.exists(out_path) == True and os.stat(out_path).st_size == 0:
            sp.check_call(call)
    job_list = moniter_rescomp_jobs.get_job_ids(job_prefix)
    moniter_rescomp_jobs.wait_for_jobs(job_list)

def sig_annot_list():
    fin1 = gzip.open(input_file,'rb')
    annot_list = fin1.readline().strip().split()[start_index:]
    fin1.close()
    sig_list = []
    for annot in annot_list:
        f = out_dir+annot+".params"
        fin=open(f,'r')
        fin.readline() # header
        fin.readline() # pi_region
        l = fin.readline().strip().split() # annotation param val list
        aname, ci_low, est, ci_hi = l[0],l[1],l[2],l[3]
        try:
            if (float(ci_hi.replace(">","")) > 0 and float(ci_low.replace("<","")) < 0) != True:
                sig_list.append(annot)
        except:
            print ("Annotation Failed: %s" % annot)
        fin.close()
    print sig_list
    return(sig_list)

def step2(sig_list):
    print("Finding annotation with highest model likelihood..")
    annot_list = sig_list
    track_dic = {}
    for annot in annot_list:
        f = out_dir+annot+".llk"
        fin=open(f,'r')
        ## Use next line for ln(lk) directly
        l = fin.readline().strip().split()
        ## Use next 3 line for AIC
        ##fin.readline()
        ##fin.readline()
        ##l = fin.readline().strip().split()
        fin.close()
        try:
            if l[0]=="ln(lk):":
        #if l[0]=="AIC:":
                track_dic[annot] = float(l[1])
        except:
            print f
    sorted_annot = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_annot.reverse() # use for ln(lk)
    top_annot = sorted_annot[0][0]
    top_val = sorted_annot[0][1]
    #print "Top annotation: %s ; ln(lk): %f" % (top_annot, top_val)
    ##print "Top annotation: %s ; AIC: %f" % (top_annot, top_val)
    for annot in sorted_annot[1:]:
        annot = annot[0]
        #print annot
        job_file = job_dir+job_prefix+top_annot+"-"+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", annot+":"+home_dir+"dist_model",
                            "-w", annot, "-o", out_dir+top_annot+"+"+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            top_annot+"+"+annot, "-o", out_dir+top_annot+"+"+annot]
        command = " ".join(command_list)
        script='''
#$ -N %s%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (job_prefix,top_annot+"-"+annot, log_dir,job_prefix+top_annot+"-"+annot,
        log_dir,job_prefix+top_annot+"-"+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        out_path = out_dir+top_annot+"+"+annot+".llk"
        if os.path.exists(out_path) == False:
            sp.check_call(call)
        if os.path.exists(out_path) == True and os.stat(out_path).st_size == 0:
            sp.check_call(call)
    job_list = moniter_rescomp_jobs.get_job_ids(job_prefix)
    moniter_rescomp_jobs.wait_for_jobs(job_list)
    top = [top_annot,top_val]
    return(top)


def run_models(fixed_list,eval_list):
    #print "Fixed annotations: " + ", ".join(fixed_list)
    fixed_name = "-".join(fixed_list)
    fixed = "+".join(fixed_list)
    for annot in eval_list:
        job_file = job_dir+job_prefix+fixed_name+"-"+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", annot+":"+home_dir+"dist_model",
                            "-w", fixed, "-o", out_dir+fixed+"+"+annot]
        elif "distance_tss" in fixed_list:
            temp_list = list(fixed_list)
            temp_list.remove("distance_tss")
            fixed_sub = "+".join(temp_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", fixed_sub+"+"+annot, "-o", out_dir+fixed+"+"+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc", "-w",
                            fixed+"+"+annot, "-o", out_dir+fixed+"+"+annot]
        command = " ".join(command_list)
        script='''
#$ -N %s%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (job_prefix,fixed+"-"+annot, log_dir,job_prefix+fixed+"-"+annot,
        log_dir,job_prefix+fixed+"-"+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        out_path = out_dir+fixed+"+"+annot+".llk"
        if os.path.exists(out_path) == False:
            sp.check_call(call)
        if os.path.exists(out_path) == True and os.stat(out_path).st_size == 0:
            sp.check_call(call)


def step3(top_annot, top_val,sig_list):

    annot_list = sig_list
    top_annot_list = top_annot.split("+")
    annot_list = [x for x in annot_list if x not in top_annot_list]
    out_list = [top_annot+"+"+ x for x in annot_list]
    run_models(top_annot_list,annot_list)
    job_list = moniter_rescomp_jobs.get_job_ids(job_prefix)
    moniter_rescomp_jobs.wait_for_jobs(job_list)
    track_dic = {}
    for name in out_list:
        f = out_dir+name+".llk"
        fin=open(f,'r')
        #fin.readline() #use for AIC
        #fin.readline() #use for AIC
        l = fin.readline().strip().split() # get 3rd line
        fin.close()
        if l[0]=="ln(lk):":
        #if l[0]=="AIC:":
            track_dic[name] = float(l[1])
    sorted_annot = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_annot.reverse() #used when comparing ln(lk)
    #print str(top_val) #+ " : " + str(floor(top_val))
    ## Next line(s) used for comparing ln(lk)
    sig_annot_list = [x for x in sorted_annot if float(x[1]) > float(top_val)]
    #sig_list = [x for x in sorted_annot if float(x[1]) > ceil(float(top_val))]
    ## Next line used for comparing AIC, comment out if unecessary
    #sig_list = [x for x in sorted_annot if float(x[1]) < floor(float(top_val))]
    #sig_list = [x for x in sorted_annot if float(x[1]) < float(top_val)]

    #print sig_annot_list
    #print len(sig_annot_list)
    if len(sig_annot_list) > 0:
        return [sig_annot_list[0], len(sig_annot_list)]
    else:
        return [top_annot, len(sig_annot_list)]



def step4(model_list):
    print "Finding the penalty with the best cross-validation likelihood..."
    p_list = ["0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50",
              "0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.0"]
    model_name = "-".join(model_list)
    model = "+".join(model_list)
    for p in p_list:
        job_file = job_dir+job_prefix+model_name+"-"+p+".sh"
        fout=open(job_file,'w')
        if "distance_tss" in model_list:
            temp_list = list(model_list)
            temp_list.remove("distance_tss")
            model_sub = "+".join(temp_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", model_sub, "-p", p, "-xv", "-print",
                            "-o", out_dir+model+"-p"+p]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            model, "-p", p, "-xv", "-print",
                            "-o", out_dir+model+"-p"+p]
        command = " ".join(command_list)
        script='''
#$ -N %s%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q long.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (job_prefix,model_name+"-p"+p, log_dir,job_prefix+model_name+"-p"+p,
        log_dir,job_prefix+model_name+"-p"+p, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        out_path = out_dir+model+"-p"+p+".ridgeparams" # try onlyp for comparison
        if os.path.exists(out_path) == False:
            sp.check_call(call)
        if os.path.exists(out_path) == True and os.stat(out_path).st_size == 0:
            sp.check_call(call)
    job_list = moniter_rescomp_jobs.get_job_ids(job_prefix)
    moniter_rescomp_jobs.wait_for_jobs(job_list)
    print "Finding best parameter value..."
    track_dic = {}
    for p in p_list:
        fin = open(out_dir+model+"-p"+p+".ridgeparams",'r')
        line_list = fin.readlines()
        fin.close()
        line = line_list[-1]
        llk = line.strip().split()[-1]
        track_dic[p]=llk
        print (p + ": " + str(llk))
    sorted_p = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_p.reverse() #used when comparing ln(lk)
    best = sorted_p[0]
    print "Optimal parameter value evaluated: %s"   % best[0]
    return best

def step5(model_list,best_p,best_llk,best_dropped_mod="NA",previously_dropped=[]):
    print "Test dropping each annotation from the model, using cross-validation likelihood"
    print "Keep dropping annotations as long as the cross-validation likelihood keeps increasing"
    if len(previously_dropped) > 0:
        dropped = "+".join(previously_dropped) + "+"
    else:
        dropped = ""
    for mod in model_list:
        keep_list = list(model_list)
        dropped_mod = mod
        keep_list.remove(mod)
        if len(keep_list) <= 9:
            qc = "short.qc"
        else:
            qc = "long.qc"
        keep_mods = "+".join(keep_list)
        job_file = job_dir+job_prefix+"drop-"+dropped+mod+".sh"
        fout=open(job_file,'w')
        if "distance_tss" in keep_list:
            keep_list.remove("distance_tss")
            model_sub = "+".join(keep_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", model_sub, #keep_mods,
                            "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+dropped+mod]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-w", keep_mods, "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+dropped+mod]
        command = " ".join(command_list)
        script='''
#$ -N %sdrop-%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q %s
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (job_prefix,dropped+mod, qc, log_dir,job_prefix+"drop-"+dropped+mod,
        log_dir,job_prefix+"drop-"+dropped+mod, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        out_path = out_dir+"drop-"+dropped+mod+".ridgeparams"
        if os.path.exists(out_path) == False:
            sp.check_call(call)
        if os.path.exists(out_path) == True and os.stat(out_path).st_size == 0:
            sp.check_call(call)
    job_list = moniter_rescomp_jobs.get_job_ids(job_prefix)
    moniter_rescomp_jobs.wait_for_jobs(job_list)
    print "The best likelihood value to beat: %s" % str(best_llk)
    track_dic = {}
    for mod in model_list:
        fin = open(out_dir+"drop-"+dropped+mod+".ridgeparams",'r')
        line_list = fin.readlines()
        fin.close()
        line = line_list[-1]
        llk = line.strip().split()[-1]
        track_dic[mod]=llk
        print ("dropped " + mod + ": " + str(llk))
    sorted_mods = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_mods.reverse() #used when comparing ln(lk)
    check_list = [x for x in sorted_mods if float(x[1]) > float(best_llk)]
    try:
        best = check_list[0]
        best_dropped_mod = best[0]
        best_dropped_llk = best[1]
        report_list = list(model_list)
        report_list.remove(best_dropped_mod)
        print ("Best dropped model: %s" % best_dropped_mod)
        print ("Best dropped llk: %s" % best_dropped_llk)
        print ("Annotations to keep: %s" % ",".join(report_list))
        status_complete = False
        best_llk = best_dropped_llk
        return best_dropped_mod,best_dropped_llk, report_list, best_llk, status_complete
    except:
        print ("Dropping models didn't improve cross-validated likelihood")
        print ("Keep the current model!")
        status_complete = True
        return best_dropped_mod,False,model_list, best_llk, status_complete


def step1_2():
    sys.stdout.write("Step 1: Running each annotation separately and identifying signficant annotations\n")
    step1()
    sig_list = sig_annot_list() # limiting to only annotations that didn't overlap zero (log2FE) from single analysis

    sys.stdout.write("Step 2: Finding the single best annotation to seed the model\n")
    top = step2(sig_list)
    top_annot, top_val = top[0], top[1]
    print "Top annotation: " + str(top_annot) + " Value: " + str(top_val)


def wrapper():
    sys.stdout.write("Step 1: Running each annotation separately and identifying signficant annotations\n")
    step1()
    sig_list = sig_annot_list() # limiting to only annotations that didn't overlap zero (log2FE) from single analysis

    sys.stdout.write("Step 2: Finding the single best annotation to seed the model\n")
    top = step2(sig_list)
    top_annot, top_val = top[0], top[1]
    print "Top annotation: " + str(top_annot) + " Value: " + str(top_val)

    sys.stdout.write("Step 3: Iteratively grow the model based on annotations that improve model likelihood\n")
    iter1 = step3(top_annot, top_val,sig_list)
    top_annot,top_val,remaining  = iter1[0][0], iter1[0][1], iter1[1]
    while remaining > 0:
        iteration = step3(top_annot, top_val,sig_list)
        print len(iteration[0])
        try:
            assert not isinstance(iteration[0], basestring)
            top_annot, top_val = iteration[0][0], iteration[0][1]
        except:
            top_annot = iteration[0]
        remaining = iteration[1]
        print "REMAINING: %d" % remaining

    sys.stdout.write("Step 4: Find penalty with best cross-valitated likelihood\n")
    model_list = top_annot.split("+")
    print top_annot
    best_p, best_llk = step4(model_list)

    sys.stdout.write("Step 5: Test dropping annotations from the model and evaluating cross-valitated likelihood\n")
    dropped_mods = []
    mod,llk,keep,best_llk,status = step5(model_list,best_p,best_llk,previously_dropped=[])
    dropped_mods.append(mod)
    print ("Dropped Models:")
    print dropped_mods
    while status == False:
        model_list.remove(mod)
        mod,llk,keep,best_llk,status = step5(model_list,best_p,best_llk,previously_dropped=dropped_mods)
        dropped_mods.append(mod)

    sys.stdout.write("Step 6: Determine the best cross-validated model\n")
    print "Here are the annotations in the best model:"
    print model_list
    pre_list = [x for x in dropped_mods if x != 'NA']
    if (len(pre_list)>0) == False:
        prefix = out_dir + "+".join(model_list)+"-p"+str(best_p)
    else:
        prefix = out_dir+"drop-"+"+".join(pre_list)
    print "Prefix of files for best model: %s" % (prefix)
    print "Copying best model input files with prefix: 'best-joint-model'"
    command = ["cp",prefix+".llk", out_dir+"best-joint-model.llk"]
    sp.check_call(" ".join(command),shell=True)
    command = ["cp",prefix+".params", out_dir+"best-joint-model.params"]
    sp.check_call(" ".join(command),shell=True)
    command = ["cp",prefix+".ridgeparams", out_dir+"best-joint-model.ridgeparams"]
    sp.check_call(" ".join(command),shell=True)
    command = ["cp",prefix+".segbfs.gz", out_dir+"best-joint-model.segbfs.gz"]
    sp.check_call(" ".join(command),shell=True)
    command = ["cp",prefix+".bfs.gz", out_dir+"best-joint-model.bfs.gz"]
    sp.check_call(" ".join(command),shell=True)
    print ("Removing all intermediate files...")
    all_files = os.listdir(out_dir)
    intermed_files1 = [x for x in all_files if "+" in x]
    intermed_files2 = [x for x in all_files if "drop" in x]
    intermed_files = list(set(intermed_files1 + intermed_files2))
    print len(all_files)
    print len(intermed_files)
    keep_files = list(set(all_files).difference(intermed_files))
    print keep_files
    print len(keep_files)
    for f in intermed_files:
        command = ["rm",out_dir+f]
        #sp.check_call(" ".join(command),shell=True)


def main():
    step1_2()

if (__name__=="__main__"): main()
