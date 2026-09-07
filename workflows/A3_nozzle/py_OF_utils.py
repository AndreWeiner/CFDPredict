'''
Created on 23.05.2023

@author: martin
'''
import random 
import os 
import shutil
import subprocess
import sys
import threading 
import ctypes
import time 
import math

verbose = True

# ----------------------------------------------------------------------
# nProcs auto-detection
# ----------------------------------------------------------------------
# Auto-detected physical-core count is capped at MAX_AUTO_NPROCS before
# being used as the solver's nProcs default. Avoids accidentally swamping
# huge servers; for explicit higher counts set STREAMLIT_NPROCS=<N> in the
# environment (systemd EnvironmentFile or shell). Cap matches the
# our production HPC server's user-policy upper bound (max 128 phys cores
# per process).
#
# OpenMPI 5.x on that server handles parallel mpirun runs cleanly --
# two jobs each with -np 128 oversubscribe via OS scheduling, they don't
# silently kill each other. (The "two jobs both crash" symptom that
# originally drove a cap=64 attempt turned out to be a German decimal
# comma typo "36,0" in an angle-list input: parsed as a 4th case at
# angle=0 -> degenerate sector geometry -> simpleFoam silently dies
# during init. MPI was never the root cause.)
MAX_AUTO_NPROCS = 128


def _detect_phys_cores(default: int = 4) -> int:
    """Return the number of physical CPU cores available on this host.

    Strategy:
      1. lscpu output (Sockets x Cores per socket) -- structured.
      2. /proc/cpuinfo fallback (parses 'physical id' + 'core id').
      3. Hard fallback: caller-supplied default.

    Counts PHYSICAL cores, not logical (no HT). Works inside VMs and
    containers since both data sources reflect the topology presented
    to the OS, which is what the OpenFOAM solver should run on.
    """
    try:
        r = subprocess.run(["lscpu"], capture_output=True, text=True, check=False)
        sockets, cores_per_socket = 0, 0
        for line in r.stdout.split("\n"):
            s = line.strip()
            if "Socket(s):" in s:
                sockets = int(s.split(":")[1].strip())
            elif "Core(s) per socket:" in s:
                cores_per_socket = int(s.split(":")[1].strip())
        if sockets > 0 and cores_per_socket > 0:
            return sockets * cores_per_socket
    except Exception:
        pass
    try:
        with open("/proc/cpuinfo") as f:
            content = f.read()
        physical_ids = set()
        core_ids = {}
        cur = None
        for line in content.split("\n"):
            if line.startswith("physical id"):
                cur = line.split(":")[1].strip()
                physical_ids.add(cur)
            elif line.startswith("core id") and cur is not None:
                core_ids.setdefault(cur, set()).add(line.split(":")[1].strip())
        if physical_ids and core_ids:
            return len(physical_ids) * max(len(c) for c in core_ids.values())
    except Exception:
        pass
    return default


def auto_nprocs(default: int = 4) -> int:
    """Return the nProcs value the solver should run with.

    Resolution order:
      1. STREAMLIT_NPROCS env var (admin override, e.g. via systemd
         EnvironmentFile=) -- if set and parseable as int, used verbatim.
      2. Auto-detected physical cores, capped at MAX_AUTO_NPROCS.
      3. `default` (caller fallback for hosts where detection fails).
    """
    env = os.environ.get("STREAMLIT_NPROCS")
    if env:
        try:
            return int(env)
        except ValueError:
            pass
    return min(_detect_phys_cores(default=default), MAX_AUTO_NPROCS)


def case_getValue(keyword, case, location, foam_dictionary):
    '''Read an OpenFOAM dictionary text file, search for the keyword and return the value.'''
    infile = open(os.path.join(case, location, foam_dictionary), "r")
    for line in infile:
        if line.strip().startswith(keyword):
            #value = list(filter(None, line.replace("\t", b' ').strip().split(";")[0].split(" ")))[1]
            value = list(filter(None, line.replace("\t", " ").strip().split(";")[0].split(" ")))[1]
            return value
    print("Warning: looking for keyword %s in %s, but nothing found!" %(keyword, os.path.join(case, location, foam_dictionary)))
    return None


def get_nproc(case_path):
    '''Read the number of processors to be used for parallel processing.'''
    return case_getValue("numberOfSubdomains", case_path, "system", "decomposeParDict")


def clean_processor_directories(case_dir):
    '''
    Remove all processor* directories to save storage.
    '''
    try:
        proc_dirs = os.listdir(case_dir)
        for each_dir in proc_dirs:
            if each_dir.startswith("processor"):
                shutil.rmtree(os.path.join(case_dir, each_dir))
    except:
        pass

def run_OF_utility(case_path, utility_name, params=[]):
    if verbose:
        print("Running %s" %utility_name)
    if not os.path.isdir(os.path.join(case_path, "logs")):
        os.makedirs(os.path.join(case_path, "logs"))
    subprocess.call([utility_name] + params + ["-case", case_path],
        stdout=open("%s/logs/%s.log" %(case_path, utility_name), "w"))

def source_OF(source_version="/opt/OpenFOAM/OpenFOAM-v2312/etc/bashrc"):
    command = ['bash', '-c', 'source %s && env' %source_version]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    
    for line in proc.stdout:
        line = line.decode("utf-8")
        try:
            if not "BASH_FUNC_mc" in line:
                line = line.split("=")
                if len(line) == 2:
                    key = line[0].strip()
                    value = line[1].strip()
                    os.environ[key] = value
                    #print("using %s" %line)
        except:
            #print("skipping %s" %line)
            pass
            
        #print line
    proc.communicate()  
    
def run_pvbatch_utility(case_path, params=[]):
    old_dir = os.getcwd()
    try:
        os.chdir(case_path)
        if verbose:
            print("Running %s" %"pvbatch")
        subprocess.call(["pvbatch"] + params)    
    except:
        print("Failure while calling pvbatch!")
    os.chdir(old_dir)
    
def execute_get_stdout(command):
    if verbose:
        print("Running %s" %command)
    p = subprocess.Popen(command.split(" "), shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = p.stdout
    stderr = p.stderr
    return stdout

def get_latestTime(case_path=os.getcwd()):
    return execute_get_stdout("foamListTimes -case %s -latestTime -withZero" %case_path).read().decode("utf-8").splitlines()[-1]

def run_solver_copy_progress(series_path, png_list, case_path, solver_name, log_file="logs/solver.log", n_proc=16, start_observer=False):
    """
    Runs the OpenFOAM solver with name solver_name using standard run_solver.
    Trys to start the logfile_observer within a thread.
    Starts a thread for copying a list of *.png files from current case dir to the series_path/progress directory.
    Useful for the cloud simulation project.
    """
    class thread_copy_pngs(threading.Thread):
        def __init__(self, case_path, progress_path, png_list):
            threading.Thread.__init__(self)
            self.png_list = png_list 
            self.png_path = os.path.join(case_path, "logs", "to_plot")
            self.case_path = case_path 
            self.progress_path = progress_path 
            self.do_continue = True 
        
        def run(self):
            try:  
                while(self.do_continue):
                    time.sleep(10)
                    print("...copy pngs...")
                    for i, each_png in enumerate(self.png_list):
                        try:
                            shutil.copy(os.path.join(self.png_path, each_png), os.path.join(self.progress_path, "%s_%s" %(i, each_png)))
                            print("... ... target file should be %s ... ..." %os.path.join(self.progress_path, "%s_%s" %(i, each_png)))
                        except:
                            pass
                print("...copy_pngs_thread ends now...")
            finally:
                print("...copy_pngs_thread is finishing now...")
            
        def get_id(self):
            if hasattr(self, '_thread_id'):
                return self._thread_id 
            for id, thread in threading._active.items():
                if thread is self:
                    return id 
        
        def end(self):
            self.do_continue = False
        
        def raise_exception(self):
            thread_id = self.get_id()
            res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, ctypes.py_object(SystemExit))
            if res > 1:
                ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
                print("Exception raise failure")
                
    progress_path = os.path.join(series_path, "progress")
    print("Starting copy_png-thread now...")
    t = thread_copy_pngs(case_path, progress_path, png_list)
    t.start()
        
    run_solver(case_path, solver_name, log_file, n_proc, start_observer)
    
    print("Ending copy_png-thread now...")
    #t.raise_exception()
    t.end()
    t.join()
        
def run_solver(case_path, solver_name, log_file="logs/solver.log", n_proc=16, start_observer=False):
    """
    Runs the OpenFOAM solver with name solver_name.
    Trys to start the logfile_observer within a thread.
    Trys to create the ParaView file suffixes.
    """  
    def which(program):
        """Find `program` on $PATH or in well-known rungui install dirs.

        is_readable instead of is_executable: rungui ships
        logfile_observer.py with the +x bit unset on many installs, but
        we invoke it via `python3 <path>` anyway so executability is
        irrelevant. The original is_exe check was the primary reason
        the live residuals chart never appeared on our production server --
        /opt/rungui/logfile_observer.py is mode 644.
        """
        def is_readable(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.R_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_readable(program):
                return program
            return None
        # 1) standard PATH
        for path in os.environ.get("PATH", "").split(os.pathsep):
            cand = os.path.join(path, program)
            if is_readable(cand):
                return cand
        # 2) RUNGUI_PATH env-var (admin override) then well-known dirs.
        # /opt/rungui-portal is the production server's stable symlink (points
        # at the rungui snapshot the portal-service is allowed to read);
        # /opt/rungui is the legacy ubuntu-streamlit / pre-migration
        # path. The env-var takes priority for non-standard installs.
        rungui_candidates = [os.environ.get("RUNGUI_PATH", ""),
                             "/opt/rungui-portal",
                             "/opt/rungui",
                             os.path.expanduser("~/rungui"),
                             "/usr/local/share/rungui"]
        for prefix in rungui_candidates:
            if not prefix:
                continue
            cand = os.path.join(prefix, program)
            if is_readable(cand):
                return cand
        return None

    def start_logfile_observer_thread(case_path):
        """
        Search for the logfile_observer (legacy single-file *or*
        rungui-next package layout) and start it as a subprocess.
        Sleep for 10 seconds so that the OpenFOAM solver can be started.
        Log clearly on miss so the next debugger doesn't have to grep
        through worker logs to find out why no live residuals chart
        shows up in the UI.
        """
        try:
            time.sleep(10)
            solver_logfile = os.path.join(case_path, "logs", "solver.log")
            observer_file = which("logfile_observer.py")
            if observer_file is not None:
                print("logfile_observer found at %s, starting observer "
                      "subprocess for %s" % (observer_file, solver_logfile))
                subprocess.call(["python3", observer_file,
                                 solver_logfile, "10"])
                return
            # rungui-next ships logfile_observer as a Python *package*
            # (with __main__.py), not a single .py file -- probe the
            # rungui-next install prefixes and invoke as `python3 -m`.
            for pkg_prefix in (os.environ.get("RUNGUI_PATH", ""),
                               "/opt/dhcae/rungui-next",
                               "/opt/rungui-next",
                               os.path.expanduser("~/dhcae/rungui-next"),
                               os.path.expanduser("~/rungui-next")):
                if not pkg_prefix:
                    continue
                if os.path.isfile(os.path.join(
                        pkg_prefix, "logfile_observer", "__main__.py")):
                    print("logfile_observer (rungui-next package) at %s, "
                          "starting `python3 -m logfile_observer` for %s"
                          % (pkg_prefix, solver_logfile))
                    env = dict(os.environ)
                    env["PYTHONPATH"] = (pkg_prefix + os.pathsep
                                         + env.get("PYTHONPATH", ""))
                    subprocess.call(["python3", "-m", "logfile_observer",
                                     solver_logfile, "10"], env=env)
                    return
            print("logfile_observer not found in $PATH or /opt/rungui / "
                  "~/rungui / /usr/local/share/rungui (legacy single-"
                  "file layout) or /opt/dhcae/rungui-next / ~/dhcae/"
                  "rungui-next (rungui-next package layout) -- live "
                  "residual / iterations charts will not appear in the "
                  "streamlit UI for this run.")
        except Exception as e:
            print("Failed launching the logfile_observer process: %s" % e)
    old_path = os.getcwd()
    os.chdir(case_path)
    if verbose:
        print("Running %s" %solver_name)
    try:
        run_OF_utility(case_path, "paraFoam", params=["-touchAll"])
    except:
        print("Failed to run paraFoam -touchAll")
    
    threads = []
    if start_observer:
        t = threading.Thread(target=start_logfile_observer_thread, args=[case_path])
        # t.daemon = True # <-- if set to True: Thread can run withouth the main program
        t.daemon = False # <-- if set to False: Thread is killed if main program exits
        threads.append(t)
        t.start()

    # Capture mpirun's exit code so silent failures (MPI init errors,
    # solver SIGTERM, prterun killed because too many ranks requested)
    # are visible in the worker.log instead of being swallowed. Logs a
    # clear marker that the next debugger will spot in the worker.log
    # tail; doesn't raise so the rest of the worker (post-processing,
    # zip, sentinel) still completes.
    rc = subprocess.call(["mpirun", "-np",
                          n_proc,
                          solver_name, "-parallel",
                          "-case", case_path],
                         stdout=open(os.path.join(case_path, log_file), "w"))
    if rc != 0:
        print("WARNING: mpirun -np %s %s -parallel returned exit code %d "
              "for case %s. Check %s and the system load -- typical cause "
              "on a shared production server is two concurrent jobs each demanding 128 "
              "ranks (solver.log stops mid-init). Solver may have produced "
              "no time directories; downstream post-processing will likely "
              "report 'Converged in 0 iterations'."
              % (n_proc, solver_name, rc, case_path, log_file))
    for each_thread in threads:
        each_thread.join()
    os.chdir(old_path)

def touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()

if __name__ == '__main__':
    pass