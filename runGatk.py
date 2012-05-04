#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator

def main():
    #os.environ["CLASSPATH"] = "/opt/jar"
    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar'

    mappings_schema = [
            {"name": "chr", "type": "string"}, 
            {"name": "lo", "type": "int32"},
            {"name": "hi", "type": "int32"},
            {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
            {"name": "ref", "type": "string"},
            {"name": "alt", "type": "string"},
            {"name": "qual", "type": "int32"},
            {"name": "coverage", "type": "int32"},
            {"name": "genotypeQuality", "type": "int32"},    
        ]

    simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi",'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    
    referenceFileName = dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    print "Indexing Dictionary"
    subprocess.check_call("java net.sf.picard.sam.CreateSequenceDictionary R=ref.fa O=ref.dict", shell=True)
    subprocess.check_call("samtools faidx ref.fa", shell=True)
    
    
    referenceDictionary = dxpy.dxlink(dxpy.upload_local_file("ref.dict"))
    referenceIndex = dxpy.dxlink(dxpy.upload_local_file("ref.fa.fai"))
    
    maxLength = 16000
    chunkSize = 16000
    
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
    
    bam = dxpy.dxlink(dxpy.upload_local_file("input.rg.bam"))
    bamIndex = dxpy.dxlink(dxpy.upload_local_file("input.rg.bam.bai"))
    
    reduceInput = {}

    for i in range(1, maxLength, chunkSize):
        intervalStart = i
        intervalEnd = min(intervalStart + chunkSize - 1, maxLength - 1)
        mapInput = {
            'bam': bam,
            'bam_index': bamIndex,
            'reference_sequence': job['input']['reference_sequence']['$dnanexus_link'],
            'reference_dict': referenceDictionary,
            'reference_index': referenceIndex,
            'from': intervalStart,
            'to': intervalEnd,
            'tableId': tableId,
            'command': buildCommand(job)
        }
        # Run a "map" job for each chunk
        mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapGatk").get_id()
        reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}


    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()
    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})
    job['output'] = {'simplevar': {'job': reduceJobId, 'field': 'simplevar'}}
    
    
def mapGatk():
    print "In GATK"
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'
    
    referenceFileName = dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    dictFileName = dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
    indexFileName = dxpy.download_dxfile(job['input']['reference_index'], "ref.fa.fai")
    
    dxpy.download_dxfile(job['input']['bam'], "input.rg.bam")
    dxpy.download_dxfile(job['input']['bam_index'], "input.rg.bam.bai")
    
    simpleVar = dxpy.open_dxgtable(job['input']['tableId'])
    
    command = job['input']['command'] + " -L chrM:" + str(job['input']['from'])+"-"+str(job['input']['to'])
    print command
    
    subprocess.call(command, shell=True)
    parseVcf(open("output.vcf", 'r'), simpleVar)


def buildCommand(job):
    
    command = "java org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.rg.bam -o output.vcf "
    command += " -stand_call_conf " +str(job['input']['call_confidence'])
    command += " -stand_emit_conf " +str(job['input']['emit_confidence'])
    command += " -pcr_error " +str(job['input']['pcr_error_rate'])
    command += " -hets " + str(job['input']['heterozygosity'])
    command += " -indelHeterozygosity " + str(job['input']['indel_heterozygosity'])
    command += " -glm " + job['input']['genotype_likelihood_model']
    command += " -mbq " + str(job['input']['minimum_base_quality'])
    command += " -maxAlleles " + str(job['input']['max_alternate_alleles'])
    command += " -deletions " + str(job['input']['max_deletion_fraction'])
    command += " -minIndelCnt " + str(job['input']['min_indel_count'])
    command += " -pnrm " + str(job['input']['non_reference_probability_model'])
    
    if job['input']['downsample_to_coverage'] != 50000:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    command += " -dt " + job['input']['downsampling_type']
    if job['input']['nondeterministic']:
        command += " -ndrs "
    print command
    return command



def reduceGatk():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    t.close(block=True)
    print "Closing Table"
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())
    

#input will require vcfFile, simpleVar, 
def parseVcf(vcfFile, simpleVar):
    #
    #compressNoCall = job['input']['compressNoCall']
    #compressReference = job['input']['compressReference']
    #storeFullVcf = job['input']['storeFullVcf']
    
    compressNoCall = False
    compressReference = False
    storeFullVcf = False
    
    header = ''

    #These prior variables are used for keeping track of contiguous reference/no-call
    #   in the event that compressReference or compressNoCall is True
    priorType = "None"
    priorPosition = -1

    mappings_schema = [
            {"name": "chr", "type": "string"}, 
            {"name": "lo", "type": "int32"},
            {"name": "hi", "type": "int32"},
            {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
            {"name": "ref", "type": "string"},
            {"name": "alt", "type": "string"},
            {"name": "qual", "type": "int32"},
            {"name": "coverage", "type": "int32"},
            {"name": "genotypeQuality", "type": "int32"},    
        ]

    #simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi",'gri')])
    #tableId = simpleVar.get_id()
    #simpleVar = dxpy.open_dxgtable(tableId)

    
    fileIter = vcfFile.__iter__()
    count = 1

    #Additional data will contain the extra format and info columns that are optional in VCF and may not be
    #   present in the VCF file. These are stored in an extended table 
    additionalData = []
    
    while 1:
        try:
            input = fileIter.next()
            if count%100000 == 0:
                print "Processed count %i variants " % count
            count += 1
            
            if input[0] == "#":
                header += input
                #extract additional column header data
                if(input[1] != "#"):
                    tabSplit = input.split("\t")
                    additionalColumns = tabSplit[7:]
                    
                        
            else:
                tabSplit = input.split("\t")
                chr = tabSplit[0]
                lo = int(tabSplit[1])
                hi = lo + len(tabSplit[3])
                ref = tabSplit[3]
                
                #In VCF format, the ALT column holds possible candidate alleles. The actual call as to the
                #   variant and its zygosity is a combination of ALT and the genotype specified in the info field.
                #   We store all of the options (including ref) and calculated the actual calls later
                altOptions = [ref.upper()]
                altOptions.extend(tabSplit[4].upper().split(","))
                qual = tabSplit[5]
                type = "Unknown"
                if qual == ".":
                    type = "No-call"
                else:
                    qual = int(float(tabSplit[5]))

                formatColumn = tabSplit[7]
                infoColumn = tabSplit[8]
                
                genotypeQuality = 0
                
                coverage = re.findall("DP=(\d+);", formatColumn)
                if(len(coverage) > 0):
                    coverage = int(coverage[0])
                else:
                    coverage = 0

                if altOptions == [ref, '.']:
                    if type == "No-call":
                        if compressNoCall == False:
                            simpleVar.add_rows([[chr, lo, hi, type, "", "", 0, 0, 0]])
                            additionalData.append(tabSplit[7:])
                    else:
                        type = "Ref"
                        if compressReference == False:
                            simpleVar.add_rows([[chr, lo, hi, type, "", "", 0, 0, 0]])
                            additionalData.append(tabSplit[7:])
                            #print [chr, lo, hi, type, "", "", 0, 0, 0]
                else:
                    #Find all of the genotypes 
                    genotypePossibilities = {}
                    for x in tabSplit[9:]:
                        genotype = getInfoField("GT", infoColumn, x)
                        genotypeQuality = float(getInfoField("GQ", infoColumn, x))
                        if genotype != False and genotypeQuality != False:
                            if genotypePossibilities.get(genotype) == None:
                                genotypePossibilities[genotype] = float(genotypeQuality)
                            else:
                                genotypePossibilities[genotype] += float(genotypeQuality)
                        else:
                            genotypeQuality = 0
                    genotypePossibilities = sorted(genotypePossibilities.iteritems(), key=operator.itemgetter(1), reverse=True)
                    genotype = genotypePossibilities[0][0]
                    genotypeQuality = genotypePossibilities[0][1]
                    if len(genotypePossibilities) > 1:
                        genotypeQuality -= genotypePossibilities[1][1]
                    alt = ""
                    if genotype == "0/0" or genotype == "0|0" or genotype == False:
                        if(len(genotypePossibilities) > 1):
                            genotype = genotypePossibilities[1][0]
                            genotypeQuality = 0
                            
                    genotypeSplit = re.split("[\|\/]", genotype)
                    for i in range(len(genotypeSplit)):
                        
                    #This is done to ensure the convention of placing the ref allele first
                    #   in practice, it seems that all VCFs already place the ref first
                        genotypeSplit[i] = int(genotypeSplit[i])
                    genotypeSplit.sort()

                    #In VCF format, the prior character to a sequence change is given in some cases (Ins, Del)
                    #   we are removing this in our format, and so need to figure out which characters to filter   
                    overlap = findMatchingSequence(ref, altOptions)

                    for x in genotypeSplit:
                        if len(alt) > 0:
                            alt += "/"
                        alt += altOptions[x][overlap:]
                        if len(altOptions[x][overlap:]) == 0:
                            alt += "-"
                        typeList = []                
                        #These rules determine how to characterize the type of change that has occurred
                        for x in altOptions:
                            if len(x) == len(ref) and len(ref) == 1:
                                type = "SNP"
                            elif ref in x:
                                type = "Ins"
                            elif x in ref:
                                type = "Del"
                            else:
                                type = "Complex"
                            for x in typeList[1::]:
                                if typeList[0] != x:
                                    type = "Mixed"
                        
                    ref = ref[overlap:]
                    if len(ref) == 0:
                        ref = "-"
                    simpleVar.add_rows([[chr, lo-overlap, lo+len(ref[overlap:]), type, ref[overlap:], alt, qual, coverage, int(genotypeQuality)]])
                    additionalData.append(tabSplit[7:])
                if compressReference:
                    if priorType == "Ref" and type != priorType:
                        simpleVar.add_rows([[chr, priorPosition, hi, type, "", "", 0, 0, 0]])
                        additionalData.append(generateEmptyList(len(additionalColumns)))
                if compressNoCall:
                    if priorType == "No-call" and type != priorType:
                        simpleVar.add_rows([[chr, priorPosition, hi, type, "", "", 0, 0, 0]])
                        additionalData.append(generateEmptyList(len(additionalColumns)))
                if type != priorType:
                    priorType = type
                    priorPosition = lo 
        except StopIteration:
            break
        
    simpleVar.set_details({"header":header})    
    #simpleVar.close(block=True)
    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})    
    #job['output']['simplevar'] = dxpy.dxlink(simpleVar.get_id())
    #
    #if storeFullVcf:
    #    extension = []
    #    for x in additionalColumns:
    #        extension.append({"name":x, "type":"string"})
    #            
    #    vcfTable = simpleVar.extend(extension)
    #    vcfTable.add_rows(additionalData)
    #    vcfTable.set_details({"header":header})    
    #    vcfTable.close(block=True)
    #    
    #    print "Full VCF table" + json.dumps({'table_id':vcfTable.get_id()})
    #    
    #    job['output']['extendedvar'] = dxpy.dxlink(vcfTable.get_id())
    
    

def findMatchingSequence(ref, altOptions):
    position = 0
    minLength = len(ref)
    for x in altOptions:
        if len(x) < minLength:
            minLength = len(x)
    for i in range(minLength):
        for x in altOptions:
            if ref[i] != x[i]:
                return i
    return minLength

def getInfoField(fieldName, infoColumn, infoContents):
    if infoColumn.count(fieldName) > 0:
        entrySplitColumn = infoColumn.split(":")
        position = -1
        for i in range(len(entrySplitColumn)):
            if entrySplitColumn[i] == fieldName:
                position = i
                entrySplitInfo = infoContents.split(":")
                if len(entrySplitInfo) == len(entrySplitColumn):
                    return entrySplitInfo[position]
    return False
    
def generateEmptyList(columns):
    result = []
    for i in range(columns):
        result.append('')
    return result
                    
