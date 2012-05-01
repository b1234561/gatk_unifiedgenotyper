#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys

def main():
    #os.environ["CLASSPATH"] = "/opt/jar"

    
    if "CLASSPATH" in os.environ: classpath = os.environ["CLASSPATH"].split(":")
    else: classpath = []
        
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'
    
    referenceFileName = dxpy.download_dxfile(job['input']['reference'], "ref.fa")
    
    
    print job['input']['sam']
    print "Downloading"
    inputFileName = dxpy.download_dxfile(job['input']['sam'], "input.sam")
    print "Converting to BAM"
    subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
    print "Sorting"
    subprocess.check_call("samtools sort input.bam input.sorted", shell=True)
    print "Adding Read Groups"
    subprocess.call("java net.sf.picard.sam.AddOrReplaceReadGroups I=input.sorted.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1", shell=True)
    print "Indexing"
    subprocess.check_call("samtools index input.rg.bam", shell=True)
    
    
    

    outputFile = dxpy.upload_local_file("input.rg.bam")
    print outputFile.describe()
    
    indexFile = dxpy.upload_local_file("input.rg.bam.bai")
    print indexFile.describe()

    referenceFile = dxpy.upload_local_file("ref.fa")
    print referenceFile.describe()


    subprocess.call("java org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.rg.bam -o output.vcf -out_mode EMIT_ALL_SITES", shell=True)

    outputFile = dxpy.upload_local_file("output.vcf")
    job['output']['rg'] = outputFile.get_id()
    print job['output']['rg'].describe()
    
    