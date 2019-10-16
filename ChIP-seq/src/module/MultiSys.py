#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     handle the file related  operation #
#          2017.12.06                    #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.1"

import sys, os
import itertools
from functools import reduce


############ commands operation ############
# accepted parameters for functions
def AcceptArgs(arg, *args):
    if arg not in args:
        args
        sys.stderr.write('Unkown parameters for "' + arg + '"\n')
        sys.exit()

# show argparse help when no aruguments
def ArgparseHelp(parser):
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

############ system operation ############
# retrive the file name and remove sufix
def DeSufix(file):
    return os.path.splitext(os.path.basename(file))[0]

# return file status
def FileExist(file, error=None):
    if file is None:
        if error is not None:
            Error('The file is NoneType\n')
    else:
        if os.path.exists(file)==False:
            if error is not None:
                Error('The file', file, 'is not existed!\n')

def Error(*messageList):
    sys.stderr.write(' '.join(messageList))
    sys.exit()

def find(directory):
    #make sure directory argument should be a directory
    directory = directory.encode('utf-8')
    assert os.path.isdir(directory)
    fileList = []
    for root,dirs,files in os.walk(directory, topdown=True):
        for fl in files:
            fileList.append(os.path.join(root,fl).decode('utf-8'))
    return fileList

############ list operation ############
def Nest2List(listA):
    return list(itertools.chain.from_iterable(listA))

# corresponding elements addition in lists
def ListsAdd(*lists):
    return list(map(sum, zip(*lists)))

# map list to dict structure
def List2Dict(listA, step=2):
    return dict(zip(listA[::step], listA[step - 1::step]))

# map list to str list
def List2Str(listA):
    return list(map(str, listA))

# map list to int list
def List2Int(listA):
    return list(map(int, listA))

# retrive the intersected elements in lists
# nestList = [list1, list2, list3 ...]
# return a set
def ListIntersect(*nestList):
    if len(nestList) == 1:
        return set(nestList[0])
    else:
        return reduce(lambda x,y: set(x) & set(y), nestList)



############ dict operation ############
# return a set
def DictKeysSet(*listDict):
    newList = []
    for i in listDict:
        newList.append(i.keys())
    if len(newList) == 1:
        return set(newList[0])
    else:
        return reduce(lambda x,y: set(x) & set(y), newList)

if __name__ == '__main__':
    print(find("/public/zhoukr/softwares/python3/lib/python3.6/multiBioPro"))
