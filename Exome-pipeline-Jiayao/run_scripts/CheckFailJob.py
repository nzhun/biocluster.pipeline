#!/home/yufengshen/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# CheckFailJob.py
#========================================================================================================

import argparse
import os

class CheckFailJob:
    def __init__(self,args):
        self.dir = args.dir
        self.prefix = args.prefix
        self.ErrorArray = self.prefix + '.ErrorArray'
        self.ErrorInput = self.prefix + '.ErrorInput'
        self.InpArrayFil = args.input
    def run(self):
        self.GetInputArray()
        ArrayFil = open(self.ErrorArray, 'wb')
        InpFil = open(self.ErrorInput, 'wb')
        Array = []
        for _file in os.listdir(self.dir):
            if _file.startswith(self.prefix):
                res = self.ProcessOneFil(_file)
                if res != False:
                    Array.append(int(res))
        Array.sort()
        ArrayFil.write("\t".join(map(str, Array)))
        for idx in Array:
            InpFil.write(self.InpArray[idx])
        ArrayFil.close()
        InpFil.close()
    def GetInputArray(self):
        self.InpArray = []
        fin = open(self.InpArrayFil, 'rb')
        for l in fin:
            self.InpArray.append(l)
        fin.close()
    def ProcessOneFil_2(self, _file):
        fin = open(_file, 'rb')
        FlagError = False
        print _file
        for l in fin:
            if "ERROR" in l:
                print l.strip()
                FlagError = True
        fin.close()
        if FlagError:
            return _file.split(".")[-1]
        else:
            return False
    def ProcessOneFil(self, _file):
        fin = open(_file, 'rb')
        FlagError = False
        print _file
        lines = fin.readlines()
        last_line = lines[-1]
        if not last_line.startswith("Done:"):
            FlagError = False
        fin.close()
        if FlagError:
            return _file.split(".")[-1]
        else:
            return False
        

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--dir',required=True, type=str, help = 'Directory that place the log file for check')
    parser.add_argument('-p','--prefix',required=True, type=str, help = 'Prefix of the log file')
    parser.add_argument('-i','--input',required=True, type=str, help = 'InputFil of the Job')
    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    ins = CheckFailJob(args)
    ins.run()
    return

if __name__=='__main__':
    main()
