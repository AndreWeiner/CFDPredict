'''
Created on Oct 13, 2012

@author: martin
'''
from __future__ import print_function
from __future__ import absolute_import
import utils

def switch_name_and_position(the_dict):
    '''
    Switch the key and the preferred position of the dictionary. For example:
        "axisPt" : ["0", ...]
        "axis" : ["1", ...]
    becomes
        "0": ["axisPt", ...]
        "1": ["axis", ... ]
    '''
    new_dict = dict()
    for key in list(the_dict.keys()):
        if the_dict[key][4]:
            subdict_selector = the_dict[key][4]
            for k, e in subdict_selector.items():
                subdict = switch_name_and_position(e[1])
                subdict_selector[k][1] = subdict
        first_entry = the_dict[key][0]
        new_dict[first_entry] = [key]
        for entry in the_dict[key][1:]:
            new_dict[first_entry].append(entry)
    return new_dict

def get_line_len(the_dict):
    '''Return the number of entries in the data part of the dict.'''
    key = list(the_dict.keys())[0]
    first_line = the_dict[key]
    return len(first_line)

def print_dict(the_dict):
    keylist = list(the_dict.keys())
    keylist.sort()
    for k in keylist:
        if the_dict[k][4]:
            subdict_selector = the_dict[k][4]
            subdict_selection = the_dict[k][1]
            subdict_name = subdict_selector[subdict_selection][0]
            subdict = subdict_selector[subdict_selection][1]
            if isinstance(subdict, dict):
                print("%s : %s" % (k, the_dict[k]))
                print(".... subdict %s ...." % subdict_name)
                print_dict(subdict)
                print(".... subdict over ....")
        else:
            print("%s : %s" % (k, the_dict[k]))
            
def print_complete_dict(the_dict):
    keylist = list(the_dict.keys())
    keylist.sort()
    for k in keylist:
        if the_dict[k][4]:
            subdict_selector = the_dict[k][4]
            for k, e in subdict_selector.items():
                print("---")
                print("%s : %s" % (k, e[0]))
                print_complete_dict(e[1])
                print("***")
        else:
            print("%s : %s" % (k, the_dict[k]))

def unroll_dict(the_dict):
    unrolled_dict = dict()
    keylist = list(the_dict.keys())
    keylist.sort()
    for k in keylist:
        if the_dict[k][4]:
            subdict_selector = the_dict[k][4]
            subdict_selection = the_dict[k][1]
            subdict_name = subdict_selector[subdict_selection][0]
            subdict = subdict_selector[subdict_selection][1]
            if isinstance(subdict, dict):
                unrolled_dict[k] = the_dict[k][1]
                unrolled_dict.update(unroll_dict(subdict))
                
        else:
            unrolled_dict[k] = the_dict[k][1]
    return unrolled_dict

def extract_subdict(dictionary_line):
    '''
    Pick the selected subdict and return it.
    '''
    subdict_selector = dictionary_line[4]
    subdict_selection = dictionary_line[1]
    if subdict_selector:
        subdict_name = subdict_selector[subdict_selection][0]
        subdict = subdict_selector[subdict_selection][1]
        if isinstance(subdict, dict):
            return subdict
    return None
        
def get_data_to_write(the_dict):
    '''
    Extract the data that will be written to the OpenFOAM dictionary.
    '''    
    data_to_write = []
    if not isinstance(list(the_dict.keys())[0], int):
        the_dict = switch_name_and_position(the_dict)
    keys = list(the_dict.keys())
    keys.sort()
    for k in keys:
        if isinstance(the_dict[k][4], dict):
            data_to_write.append("%s %s;" % (the_dict[k][0], the_dict[k][1]))
            sub_dict = the_dict[k][4]
            choice = the_dict[k][1]
            sub_template = sub_dict[choice]
            data_to_write.append(sub_template[0])
            if isinstance(sub_template[1], dict):
                data_to_write.append("{")
            else:
                print("Warning: unrecognized sub_template %s" % sub_template[1])
            insertion = get_data_to_write(sub_template[1])
            for line in insertion:
                data_to_write.append(line)
            if isinstance(sub_template[1], dict):
                data_to_write.append("}")            
        elif isinstance(the_dict[k][1], list):
            line = "%s (" % the_dict[k][0]
            for item in the_dict[k][1]:
                line += " %s" %item
            line += " );"
            data_to_write.append(line)
        else:
            entry = utils.win_to_unix_path(the_dict[k][1])            
            data_to_write.append("%s %s;" % (the_dict[k][0], entry))
    return data_to_write

def get_options_string(options_dict):
    '''
    Extract the options from the options_dict and return them to be appended to 
    the command string.
    '''
    def append_option(the_option):
        '''If necessary appends the option to the command string'''
        if isinstance(options_dict[the_option][1], bool):
            if options_dict[the_option][1]:
                return " -" + the_option 
        else:
            entry = options_dict[the_option][1]
            if entry != "":
                entry = utils.win_to_unix_path(entry)
                return " -" + the_option + " " + entry 
        return ""
    option_string = ""
    if isinstance(list(options_dict.keys())[0], int):
        options_dict = switch_name_and_position(options_dict)        
    for the_option in list(options_dict.keys()):
        option_string += append_option(the_option)
    return option_string

def get_position_version(the_dict):
    '''Calls switch_name_and_position if necessary.'''
    if not isinstance(list(the_dict.keys())[0], int):
        the_dict = switch_name_and_position(the_dict)
    return the_dict

def force_position_version(the_dict):
    '''Traverse through the subdicts and guarantee the position version.'''
    new_dict = dict()
    for key in list(the_dict.keys()):
        if the_dict[key][4]:
            subdict_selector = the_dict[key][4]
            for k, e in subdict_selector.items():
                subdict = force_position_version(e[1])
                subdict_selector[k][1] = subdict            
        if not isinstance(key, int):
            first_entry = the_dict[key][0]
            new_dict[first_entry] = [key]
            for entry in the_dict[key][1:]:
                new_dict[first_entry].append(entry)
        else:
            new_dict[key] = [first_entry]
            for entry in the_dict[key][1:]:
                new_dict[key].append(entry)
    return new_dict

def get_key_version(the_dict):
    '''Calls switch_name_and_position if necessary.'''
    if len(the_dict) > 0:
        if isinstance(list(the_dict.keys())[0], int):
            the_dict = switch_name_and_position(the_dict)
    return the_dict

def get_arguments_string(arguments_dict):
    '''
    Extract the arguments from the arguments_dict and return them to be appended to 
    the command string.
    If more than one argument is necessary, the sequence is determined by the 
    integer key of the arguments_dict.
    '''
    argument_string = ""
    argument_list = list()
    if not isinstance(list(arguments_dict.keys())[0], int):
        arguments_dict = switch_name_and_position(arguments_dict)
    keys = list(arguments_dict.keys())
    keys.sort()
    for key in keys:
        argument_list.append(arguments_dict[key][0])
    
    arguments_dict = switch_name_and_position(arguments_dict)
    
    for each_arg in argument_list:
        entry = arguments_dict[each_arg][1]
        entry = utils.win_to_unix_path(entry)
        argument_string += " " + entry
    return argument_string

if __name__ == '__main__':
    from modules import foam_extrudeMesh
    the_dict = foam_extrudeMesh.dictionary_dict
    #print_dict(the_dict)
    #the_dict = force_position_version(the_dict)
    print_complete_dict(the_dict)
    print() 
    the_dict = get_position_version(the_dict)
    print_complete_dict(the_dict)
    print()
    the_dict = get_key_version(the_dict)
    print_complete_dict(the_dict)
    print() 
    the_dict = get_position_version(the_dict)
    print_complete_dict(the_dict)
    print()