import argparse
import sys
import os

def readArguments():
    parser = argparse.ArgumentParser(description='Create exercise sheet.')
    parser.add_argument('sourceFolder', type=str, help='Source folder of .geo file')
    parser.add_argument('meshNamePrefix', type=str,
                        help='n/y, how often should characteristic length be halved?')
    parser.add_argument('-l', '--lcHalvesteps', type=int, default=0,
                        help='How often should characteristic length be halved?')
    parser.add_argument('-d', '--destination', type=str, default='../../build/',
                        help='Destination folder')
    args = parser.parse_args()
    

    return args

if __name__ == '__main__':
    args = readArguments()
    for k in range(args.lcHalvesteps + 1):
        # create meshs
        geoName = args.meshNamePrefix + '.geo'
        mshName = '%s_%d.msh'%(args.meshNamePrefix, k) 
        xmlName = '%s_%d.xml'%(args.meshNamePrefix, k) 
        os.system('gmsh %s -2 -clscale %1.2f -o %s'%(args.sourceFolder + geoName, 0.5**k,args.sourceFolder + mshName))
        # move meshs
        os.system('mv %s %s'%(args.sourceFolder + mshName, args.destination + mshName))
        # convert
        os.system('dolfin-convert %s %s'%(args.destination + mshName, args.destination + xmlName))



