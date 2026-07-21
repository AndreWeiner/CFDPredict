'''
Created on Oct 28, 2012

@author: martin
'''
from __future__ import absolute_import
from __future__ import print_function
try:
    from future import standard_library
    standard_library.install_aliases()
except:
    pass # probably python 2.7
#from builtins import str
#from past.builtins import basestring
try:
    from builtins import str
except:
    pass # python2
    #from __builtin__ import str
import os
#from . import globaldata
#import modules.globaldata as globaldata
import sys
import subprocess
import re
#from Crypto.Cipher import AES
#import base64 

def get_application_path():
    '''
    Return the path of the running application.
    Works with a frozen executable or a python script running.
    '''
    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)
    return application_path

def get_active_OpenFOAM_version():
    '''
    Check, which OpenFOAM version is currently active.
    '''
    if "WM_PROJECT_VERSION" in os.environ:
        return os.environ["WM_PROJECT_VERSION"]
    else:
        return None
    
def get_path_to_rungui():
    '''
    Return the path to the location of rungui.
    '''
    own_app_path = get_application_path()
    if own_app_path.endswith("modules"):
        own_app_path = os.path.abspath(os.path.join(own_app_path, os.pardir)) 
    return own_app_path

def get_path_to_modules():
    '''
    Return the path to the rungui/modules location.
    '''
    return os.path.join(get_path_to_rungui(), "modules")       

def get_path_to_rungui_app():
    '''
    Return the path to the rungui application.
    '''
    return os.path.join(get_path_to_rungui(), "rungui.py")       

def is_linux():
    return sys.platform.startswith("linux")

def linux_is_pid_accessing_filename(pid, full_checkfilename):
    '''
    For linux only: check, if a process has a file handle to a specific file.
    '''
    if not is_linux():
        return False
    dir = '/proc/'+str(pid)+'/fd'
    if not os.access(dir,os.R_OK|os.X_OK): return False

    for fds in os.listdir(dir):
        for fd in fds:
            full_name = os.path.join(dir, fd)
            try:
                the_file = os.readlink(full_name)
                if the_file == '/dev/null' or \
                  re.match(r'pipe:\[\d+\]',the_file) or \
                  re.match(r'socket:\[\d+\]',the_file):
                    the_file = None
                    continue
            except OSError as err:
                if err.errno == 2:     
                    the_file = None
                    continue
                else:
                    raise(err)
                    continue
            
            if the_file == full_checkfilename:
                return True
    return False
            #yield (fd,file)

def win_get_pid_dict(key_is_name=True):
    '''
    Return the process table as dictionary.
    '''
    if not is_linux():
        try:
            task_dict = dict()
            #task_list = os.popen(r'tasklist /FO "CSV" 2>&1' , 'r')
            si = subprocess.STARTUPINFO()
            si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            task_list = subprocess.check_output(r'tasklist /FO CSV', startupinfo=si).splitlines()                   
            for line in task_list: 
                process_name, process_pid = line.split(",")[0:2]
                if "\"" in process_name:
                    process_name = process_name.strip("\"")
                if "\"" in process_pid:
                    process_pid = process_pid.strip("\"")
                    
                if key_is_name:
                    if process_name in task_dict:
                        task_dict[process_name].append(process_pid)
                    else:
                        task_dict[process_name] = [process_pid]
                else:
                    task_dict[process_pid] = process_name
            #task_list.close()
            return task_dict
        except:
            pass
    return {}

def is_PID_running(pid):
    if is_linux():
        return os.path.exists("/proc/%s" % pid)
    else:
        task_dict = win_get_pid_dict(key_is_name=False)
        return pid in task_dict

def check_PID_and_process_name(process_name, pid):
    '''
    Check, if process_name is running with PID.
    '''
    if is_linux():
        path_to_process = "/proc/%s/cmdline" % pid
        if os.path.isfile(path_to_process):
            p_name_file = open(path_to_process, "r")
            return process_name in p_name_file.read()
    else:
        pid_dict = win_get_pid_dict(key_is_name=True)
        if process_name in pid_dict:
            return pid in pid_dict[process_name]
    return False
           
def get_PIDs_for_process_name(process_name):
    '''
    Return a list of all PIDs that contain the process_name.
    '''
    if is_linux():
        all_pids = os.listdir("/proc")
        pid_list = []
        for each_pid in all_pids:
            try:
                path_to_process = "/proc/%s/comm" % each_pid
                p_name_file = open(path_to_process, "r")
                if process_name in p_name_file.read():
                    pid_list.append(each_pid)
            except:
                pass
        return pid_list
    else:
        pid_dict = win_get_pid_dict(key_is_name=False)
        pid_list = []
        for key, value in pid_dict.items():
            if process_name in value:
                pid_list.append(key)
        return pid_list
            
           
def decode(e):
    try:
        keylist = [
                   "2013-12-31: S4SPP9eOKSUrZGM5sLU3JQ==",
                   "2014-04-01: LcdS0dDUq1FmngVHOoWB4g==",
                   "2014-10-01: vEZYreweWfTUbnUnFort0w==",
                   "2015-04-01: 8FcIwxgM3RJ87Emm+CFp5Q==",
                   "2015-10-01: D5bB3VHcX9F2onVJfo4hyw==",
                   "2016-04-01: jcLZ4kc3LXcty8AGJFiSoQ==",
                   "2016-10-01: CGI2YrAkynNzaYlNhOYqdg==",
                   "2017-04-01: FxWw3EDUbML7VdUrfk2/dw==",
                   "2017-10-01: bbne033R73XKkH0zNzOpkA==",
                   "2018-04-01: GwBBZFnF2zE5CgM+LSac4Q==",
                   "2018-10-01: peig81yHq3GDj9rBwD9v1w==",
                   "2019-04-01: PqQwldxWOmnSxDAZP1xFRg==",
                   "2019-10-01: gGoIf/4wr79oUxuOZ4KjCw==",
                   "2020-04-01: XGSocSmZVj4UZokEpC6LrQ==",
                   "2020-10-01: hbaLyNlmp/U8Ej01LZZzzQ==",
                   "2021-04-01: /FDpOMl0qOCxzL1nPMxHew==",
                   "2021-10-01: atxaKmqtIB8MFh/ubzInBg==",
                   "2022-04-01: +Mwtw2nZsdhNf5WOGQ/EDw==",
                   "2022-10-01: fBU3TJVheEH2Lb8BtzYZ8Q==",
                   "2023-04-01: kNCw53iuqMML00xnkHZ9mQ==",
                   "2023-10-01: 0d5gmsA9NZc3ToWGBYFtJw==",
                   "2024-04-01: RPrv8CIXg/wkZvY14qy+9w==",
                   "2024-10-01: Vy1QhZ+MNY3x987Y9fdazA==",
                   "2025-04-01: B6sIRX+qX4945eEUyboZQw==",
                   "2025-10-01: yknLM5V2DDhOuqeVQOSZ3w==",
                   "2026-04-01: e0ykBJo9r5kYFVolis0wOQ==",
                   "2026-10-01: DTFRTXrpFODHqtb3CV6w/Q==",
                   ]
        for entry in keylist:
            d,k = [s.strip() for s in entry.split(":")]
            if e == k:
                if len(d.split("-")) == 3:
                    return [int(i) for i in d.split("-")[0:2]] 
    except:
        pass
    return [0,0]
                
def get_castnet_command():    
    path = os.environ['PATH'].split(os.pathsep)
    path_to_simapps = None
    
    for item in path:
        try:
            if os.path.exists(os.path.join(item,"simapps")):
                path_to_simapps = item
                break
        except:
            pass
    if path_to_simapps:
        path_to_castnet = os.path.join(path_to_simapps, "apps", "castNet", "MANIFEST")        
        if os.path.isfile(path_to_castnet):
            path_to_simapps = os.path.join(path_to_simapps, "simapps")
            return "%s %s %s" %(path_to_simapps, path_to_castnet, "--cae ")
    return None

def find_castnet_on_windows():
    '''
    Search for the CastNet entry in the Start Menu and get the installation directory.
    '''     
    import winreg
    import win32com.client
    import tempfile
    
    def get_castnet_in_StartPrograms():
        '''Lookup the Start->Programs folder and search for the CastNet installation. Return path to CastNet in Start->Programs.'''
        objShell = win32com.client.Dispatch("WScript.Shell")
        allUserProgramsMenu = objShell.SpecialFolders("AllUsersPrograms")
        all_entries = os.listdir(allUserProgramsMenu)
        castnet_entries = [n for n in all_entries if n.startswith("CastNet")]
        castnet_entries.sort()
        if len(castnet_entries) > 0:
            return os.path.join(allUserProgramsMenu, castnet_entries[-1])
        
        userMenu = objShell.SpecialFolders("StartMenu")
        all_entries = os.listdir(userMenu)
        castnet_entries = [n for n in all_entries if n.startswith("CastNet")]
        castnet_entries.sort()
        if len(castnet_entries) > 0:
            return os.path.join(userMenu, castnet_entries[-1])
        return None
    
    def get_castnet_link():
        path_to_castnet = get_castnet_in_StartPrograms()
        if path_to_castnet:
            path_to_castnet_link = os.path.join(path_to_castnet, "CastNet.lnk")
            if os.path.isfile(path_to_castnet_link):
                return path_to_castnet_link
        return None
    
    def get_path_to_simapps():
        path_to_castnet_link = get_castnet_link()
        if path_to_castnet_link:
            shell = win32com.client.Dispatch("WScript.Shell")
            shortcut = shell.CreateShortCut(path_to_castnet_link)
            target_file = shortcut.Targetpath
            if os.path.isfile(target_file):
                return os.path.dirname(target_file)
        return None
            
    def get_castnet_command():
        exec_path = get_path_to_simapps()
        if not exec_path:
            return
        manifest_file = os.path.join(os.path.abspath(os.path.join(exec_path, os.pardir)), "apps", "castNet", "MANIFEST")
        path_to_simapps = os.path.join(exec_path, "simapps.bat")
        castnet_command = "\"%s\" \"%s\" %s" %(path_to_simapps, manifest_file, "--cae")
        return castnet_command, exec_path
        

    return get_castnet_command()

def expand_path(the_path):
    if the_path == ".":
        the_path = os.getcwd()
    elif the_path == "~":
        the_path = os.path.expanduser('~')
    return the_path

def win_to_unix_path(path):
    r'''
    checks if path is in Windows format, e.g. "C:\workwith\test.txt", and transfers it back to
    UNIX format, i.e. "C:/workwith/test.txt". This is necessary for the blueCFD msys environment.
    '''
    def repair_drive_letter(path):
        if ":\\\\" in path:
            path = path.replace(":\\\\", ":/")
        if ":\\" in path:
            path = path.replace(":\\", ":/")
        return path
    from . import globaldata
    if globaldata.settings["transform_Windows_path_to_UNIX_path"] and type(path) is str: #isinstance(path, basestring):
        if "\\" in path or ":" in path:
            path = repair_drive_letter(path)            
            path = path.replace('\\', '/')
    return path

def win_slash_path(path):
    '''
    Take a path like C:/Program Files (x86)/RunGui/logfile_observer.exe and return
    /c/Program Files (x86)/RunGui/logfile_obersver.exe
    '''
    if ":/" in path:
        drive = path.split(":")[0]
        path = path.split(":")[1]
        path = path.replace(" ", r"\ ")
        return "/%s%s" %(drive.lower(), path)
    elif ":\\" in path:
        drive = path.split(":")[0]
        path = path.split(":")[1]
        path = path.replace(" ", r"\ ")
        return "%s:%s" %(drive, path)
    else:
        return path

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def swap_dict(the_dict):
    '''swap key and entry[0]'''
    swapped_dict = dict()
    for k in list(the_dict.keys()):
        entry = the_dict[k]
        if entry in swapped_dict:
            swapped_dict[entry].append(k)
        else:
            swapped_dict[entry] = [k]
    return swapped_dict

def folderToZip(path):
    '''Make an archive of folder at path.'''
    import zipfile
    if os.path.isdir(path):
        try:
            dirname1 = os.path.basename(path) 
            archive_name = os.path.abspath(os.path.join(os.path.dirname(path), dirname1+".zip"))
            relroot = os.path.abspath(os.path.join(path, os.pardir))
            with zipfile.ZipFile(archive_name, "w", zipfile.ZIP_DEFLATED) as zip:
                for root, dirs, files in os.walk(path):
                    # add directory (needed for empty dirs)
                    zip.write(root, os.path.relpath(root, relroot))
                    for file in files:
                        filename = os.path.join(root, file)
                        if os.path.isfile(filename): # regular files only
                            arcname = os.path.join(os.path.relpath(root, relroot), file)
                            zip.write(filename, arcname)
        except:
            print("Warning: can't archive %s" %path)
    else:
        print("Warning: %s is not a directory" %path)


def folderToTarGz(path):
    '''Make an archive of folder at path.'''
    if not is_linux():
        folderToZip(path)
        return
    import tarfile
    if os.path.isdir(path):
        try:
            dirname1 = os.path.basename(path) 
            archive_name = os.path.abspath(os.path.join(os.path.dirname(path), dirname1+".tar.gz"))
            tar = tarfile.open(archive_name, 'w:gz')            
            tar.add(path, arcname=dirname1)
            tar.close()
        except:
            print("Warning: can't archive %s" %path)
    else:
        print("Warning: %s is not a directory" %path)

def packArchive(dirname):
    '''Create a tar.gz archive.'''    
    try:
        import tarfile
        with tarfile.open(dirname + ".tar.gz", "w:gz") as tar:
            tar.add(dirname, arcname=os.path.basename(dirname))        
    except:
        pass
    
def unpackArchive(filename):
    '''Unpack the selected archive.'''
    try:
        import tarfile
        if filename.endswith(".gz"):
            tar = tarfile.open(filename)
            tar.extractall()
            tar.close()
    except:
        pass

def get_terminal():
    '''Test different terminals and return the standard terminal application, pe. xterm.'''
    term_list = [
                 "xterm",
                 "gnome-terminal",
                 "konsole"
                 ]
    for each_terminal in term_list:
        if which(each_terminal):
            return each_terminal
    return None

if __name__ == '__main__':
    pass
    