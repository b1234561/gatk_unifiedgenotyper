#!/usr/bin/env python

import dxpy
import subprocess, logging

def main():
    print job['input']['sam']
    inputFileName = dxpy.download_dxfile(job['input']['sam'], "input.bam")
    subprocess.check_call("java -jar AddOrReplaceReadGroups.jar I=input.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1")
    job['output']['rg'] = dxpy.upload_local_file("output.rg.bam")
    print job['output']['rg'].describe()