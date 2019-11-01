"""Change reference to paper in all files."""

import os
import shutil
import checkMatlabHelp as chM
import numpy as np

class checker():
    def __init__(self, pattern):
        self.pattern = pattern
        self.no_lines = len(self.pattern)

    def check(self, window):
        if len(window) != self.no_lines:
            raise ValueError('Pattern has %d number of lines, you gave %d'%(self.no_lines, len(window)))
    
        for k in range(self.no_lines):
            if self.pattern[k].strip('\n').strip() != window[k].strip('\n').strip():
                return False
    
        # if we arrive here, pattern coincides
        return True

def find_line_to_replace(lines, pattern):
    line_idx = np.nan
    ch = checker(pattern)
    for l in range(len(lines) - ch.no_lines + 1):
        window = lines[l : l + ch.no_lines]
        if ch.check(window):
            line_idx = l
            break
        
    return line_idx

def rewrite_file(file_name, pattern_new, line_idx, len_pattern_old):
    print('rewrite: ', file_name)
    tmp = file_name.rsplit('.', 1)
    cache_name = tmp[0]+'_cache' + '.' + tmp[1]
    shutil.copyfile(file_name, cache_name)
    if not np.isnan(line_idx):
        with open(cache_name, 'r') as f:
            lines = f.readlines()
            with open(file_name, 'w') as f_new:
                for k in range(line_idx):
                    f_new.write(lines[k])
                for k in range(len(pattern_new)):
                    f_new.write(pattern_new[k] + '\n')
                for k in range(line_idx + len_pattern_old, len(lines)):
                    f_new.write(lines[k])
    os.remove(cache_name)
        
def get_list_of_matlab_python_files(all_files):
    """Get list of matlab and python files."""
    all_m_py_files = []
    for f in all_files:
        if f.split('/')[-1].endswith('.m') or f.split('/')[-1].endswith('.py'):
            all_m_py_files.append(f)
    
    return all_m_py_files     


if __name__ == '__main__':
    pattern = ['%', '%#$&replaceMe&$#%', '%']
    pattern_new  = ['%', '%See The second-order formulation of the P_N equations with Marshak boundary conditions ', '%Matthias Andres, Florian Schneider',  '%']
    
#    k = 0
#    k_max = 3
    all_files = chM.getListOfFiles('./')
    all_files.remove('./changeRef.py')
    all_m_py_files = get_list_of_matlab_python_files(all_files)
    for file_name in all_m_py_files:
#        k  += 1
#        if k == k_max:
#            break
        with open(file_name, 'r') as f:
            lines = f.readlines()
            line_idx = find_line_to_replace(lines, pattern)
        rewrite_file(file_name, pattern_new, line_idx, len(pattern))
        

