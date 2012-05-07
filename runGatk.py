#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator

def main():
    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar'

    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['reference_contig_set']['$dnanexus_link']), shell=True)
    reference_sequence = dxpy.dxlink(dxpy.upload_local_file("ref.fa"))
    
    #referenceFileName = dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    print "Indexing Dictionary"
    subprocess.check_call("java net.sf.picard.sam.CreateSequenceDictionary R=ref.fa O=ref.dict", shell=True)
    subprocess.check_call("samtools faidx ref.fa", shell=True)
    
    referenceDictionary = dxpy.dxlink(dxpy.upload_local_file("ref.dict"))
    referenceIndex = dxpy.dxlink(dxpy.upload_local_file("ref.fa.fai"))
    
    
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
    
    ##Run the trivial case to find the columns
    trivialOutput = runTrivialTest(job['input']['reference_contig_set'], buildCommand(job))
    
    for x in trivialOutput['additionalColumns']:
        mappings_schema.append({"name": x.strip(), "type": "string"})
    #simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri'), dxpy.DXGTable.substring_index("type", "typeIndex")])
    simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details({'header':trivialOutput['header']})
    
    reduceInput = {}

    commandList = splitGenomeLength(job['input']['reference_contig_set'], job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
        
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:         
            mapInput = {
                'bam': bam,
                'bam_index': bamIndex,
                'reference_sequence': reference_sequence,
                'reference_dict': referenceDictionary,
                'reference_index': referenceIndex,
                'interval': commandList[i],
                'tableId': tableId,
                'command': buildCommand(job),
                'compress_reference': job['input']['compress_reference'],
                'compress_no_call' : job['input']['compress_no_call'],
                'store_full_vcf' : job['input']['store_full_vcf']
                
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
    os.environ['CLASSPATH'] = '/opt/jar/GenomeAnalysisTK.jar'
    
    referenceFileName = dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    dictFileName = dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
    indexFileName = dxpy.download_dxfile(job['input']['reference_index'], "ref.fa.fai")
    
    dxpy.download_dxfile(job['input']['bam'], "input.rg.bam")
    dxpy.download_dxfile(job['input']['bam_index'], "input.rg.bam.bai")
    
    simpleVar = dxpy.open_dxgtable(job['input']['tableId'])
    
    command = job['input']['command'] + job['input']['interval']
    print command
    
    subprocess.call(command, shell=True)
    parseVcf(open("output.vcf", 'r'), simpleVar, job['input']['compress_reference'], job['input']['compress_no_call'], job['input']['store_full_vcf'])


def buildCommand(job):
    
    command = "java org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.rg.bam -o output.vcf "
    command += " -out_mode " + (job['input']['output_mode'])
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

def runTrivialTest(contig_set, command):
    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    for i in range(len(names)):
        if sizes[i] > 0:
            chromosome = names[i]
            break
    command += ' -L ' + chromosome+':1-1'
    subprocess.call(command, shell=True)
    return extractHeader(open("output.vcf", 'r'))
    

def extractHeader(vcfFile):
    header = ''
    fileIter = vcfFile.__iter__()

    #Additional data will contain the extra format and info columns that are optional in VCF and may not be
    #   present in the VCF file. These are stored in an extended table 
    additionalColumns = []
    while 1:
        try:
            input = fileIter.next()
            if input[0] == "#":
                header += input
                #extract additional column header data
                if(input[1] != "#"):
                    tabSplit = input.split("\t")
                    return {'header': header, 'additionalColumns': tabSplit[7:]}
        except StopIteration:
            break
    
def reduceGatk():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    t.close(block=True)
    print "Closing Table"
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())
    
def parseVcf(vcfFile, simpleVar, compressNoCall, compressReference, storeFullVcf):

    #These prior variables are used for keeping track of contiguous reference/no-call
    #   in the event that compressReference or compressNoCall is True
    priorType = "None"
    priorPosition = -1

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
            
            if input[0] != "#":
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
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.extend(tabSplit[7:])
                            simpleVar.add_rows([entry])
                    else:
                        type = "Ref"
                        if compressReference == False:
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.extend(tabSplit[7:])
                            simpleVar.add_rows([entry])
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
                    entry = [chr, lo-overlap, lo+len(ref[overlap:]), type, ref[overlap:], alt, qual, coverage, int(genotypeQuality)]
                    entry.extend(tabSplit[7:])
                    simpleVar.add_rows([entry])
                if compressReference:
                    if priorType == "Ref" and type != priorType:
                        entry = [chr, priorPosition, hi, type, "", "", 0, 0, 0]
                        entry.extend(tabSplit[7:])
                        simpleVar.add_rows([entry])                        
                if compressNoCall:
                    if priorType == "No-call" and type != priorType:
                        entry = [chr, priorPosition, hi, type, "", "", 0, 0, 0]
                        entry.extend(tabSplit[7:])
                        simpleVar.add_rows([entry])
                if type != priorType:
                    priorType = type
                    priorPosition = lo  
        except StopIteration:
            break

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
    
def splitGenomeLength(contig_set, includeInterval, excludeInterval, chunkSize, splits):
    print contig_set
    print dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']
    
    print names
    print sizes
    print offsets
    
    commandList = []
    position = 0
    chromosome = 0
    currentChunk = 0
    
    for i in range(splits):
        commandList.append(" "+excludeInterval)
        
    includeDictionary = {}
    includeMatch = re.findall("(\w+):(\d+)-(\d+)", includeInterval)
    for x in includeMatch:
        if includeDictionary.get(x[0]) == None:
            includeDictionary[x[0]] = []
            includeDictionary[x[0]].append([int(x[1]), int(x[2])])
    
    while chromosome < len(names):
        if position + chunkSize >= sizes[chromosome]:
            print chromosome
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, sizes[chromosome]+1)
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, position+chunkSize+1)
            position += chunkSize
        currentChunk = (currentChunk+1)%splits
 
    return commandList
    
def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            print "List"
            print x
            min = lo
            max = hi
            if (lo >= x[0] and lo <= x[1]) or (hi <= x[1] and hi >= x[0]):
                if lo >= x[0] and lo <= x[1]:
                    min = lo
                elif lo <= x[0]:
                    min = x[0]
                if hi <= x[1] and hi >= x[0]:
                    max = hi
                elif hi >= x[1]:
                    max = x[1]
                command += " -L %s:%d-%d" % (chromosome, min, max)
    return command

    
    
    
    
    
    
    
    
    