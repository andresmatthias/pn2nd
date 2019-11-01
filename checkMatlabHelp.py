"""Check if references in Matlab functions to its unit tests are correct."""

import os
import pprint
import numpy as np

def getListOfFiles(dirName):
    # from https://thispointer.com/python-how-to-get-list-of-files-in-directory-and-sub-directories/
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles

def getListOfMatlabTestFiles(allMFiles):
    """Get list of matlab unit test files."""
    allTestFiles = []
    for f in allFiles:
        if f.split('/')[-1].endswith('Test.m'):
            allTestFiles.append(f)
    
    return allTestFiles


def getTestsPerFunction(allScriptNames, allTestFiles):
    """ Assign to each m-file a list of tests which call this file.
        As we look for exact correspondences and not for the Upper-case 
        version, we avoid that within a test file its own name is recognized 
        at the beginning of its help section.
    """
    
    file2test = dict()
    for m in allScriptNames:
        file2test[m] = []
        for t in allTestFiles:
            if checkIfFileIsTested(m, t):
                file2test[m].append(t)
    
    for f in file2test:
        file2test[f] = list(set(file2test[f]))
    return file2test 


def checkMatlabRef(file2test, allMFiles, allTestFiles):
    """Check if all files called by a unit test give a reference to this."""
    report = dict()
    for fName in allMFiles:
            nameShort = fName.split('/')[-1].split('.')[0]
            idx = np.ones(len(file2test[nameShort]))
            with open(fName) as f:
                lines = f.readlines()
            
            for l in lines:
                for k, t in enumerate(file2test[nameShort]):
                    tShort = t.split('/')[-1].split('.')[0]
                    if tShort.upper() in l:
                        idx[k] = 0
            
            if np.sum(idx) > 0:
                report[fName] = [file2test[nameShort][k] for k in range(len(idx)) if idx[k] > 0]
        
    return report

def checkIfFileIsTested(scriptName, testFileName):
    """Is a file tested within a given test file?"""
    with open(testFileName) as f:
        lines = f.readlines()
    for l in lines:
        if scriptName in l:
            return True
    return False


def checkReferencedUnitTests(allMFiles, allScriptNames, allTestFiles):
    """Check if a function references to a unit test, that doesn't call this function."""
    report = dict()
    for fName in allMFiles:
        nameShort = fName.split('/')[-1].split('.')[0]
        with open(fName) as f:
                lines = f.readlines()
        lines = lines[2:] # skip first and second row, containing the name of the function
      
        for l in lines:
            for t in allTestFiles:
                 tShort = t.split('/')[-1].split('.')[0]
                 if tShort.upper() in l:
                     if not checkIfFileIsTested(nameShort, t):
                         try:
                             report[fName].append(t)
                         except KeyError:
                             report[fName] = [t]
    return report
        
    

if __name__ == '__main__': 
    allFiles = getListOfFiles('./')
    allMFiles = [f for f in allFiles if f.endswith('.m')]
    allScriptNames = [f.split('/')[-1].split('.')[0] for f in allMFiles]
    allTestFiles = getListOfMatlabTestFiles(allMFiles)
    
    # check if unit tests calling a function are also referenced in it
    file2test = getTestsPerFunction(allScriptNames, allTestFiles)
    report1 = checkMatlabRef(file2test, allMFiles, allTestFiles)
    if len(report1) == 0:
        print('\n------\nReferences on unit tests in Matlab all correct!\n------\n')
    else:
        print('\n------\nFunctions appear in the following unit tests:\n------\n')
        pprint.pprint(report1)
        
    #  check if unit test referenced in a function indeed calls that function
    report2 = checkReferencedUnitTests(allMFiles, allScriptNames, allTestFiles)
    if len(report2) == 0:
        print('\n------\nNo function references to a unit test, which does not call that function :) \n------\n')
    else:
        print('\n\n\n------\nFollowing unit tests are referenced, but not used::\n------\n')
        pprint.pprint(report2)
