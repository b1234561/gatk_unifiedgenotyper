#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator

def main():
    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()
    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")
    
    

    variants_schema = [
            {"name": "chr", "type": "string"}, 
            {"name": "lo", "type": "int32"},
            {"name": "hi", "type": "int32"},
            {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
            {"name": "ref", "type": "string"},
            {"name": "alt", "type": "string"},
            {"name": "qual", "type": "int32"},
            {"name": "coverage", "type": "int32"},
            {"name": "genotype_quality", "type": "int32"},    
        ]
    if job['input']['store_full_vcf']:
        variants_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])
    
    simpleVar = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details({'original_contigset':originalContigSet})
    
    
    reduceInput = {}
    commandList = splitGenomeLength(originalContigSet, job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
        
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:
            mapInput = {
                'mappings_table_id':mappingsTableId,
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'command': buildCommand(job),
                'compress_reference': job['input']['compress_reference'],
                'compress_no_call' : job['input']['compress_no_call'],
                'store_full_vcf' : job['input']['store_full_vcf'],
                'part_number' : i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapGatk").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}
            dxpy.DXJob(mapJobId).wait_on_done()

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()
    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})
    job['output'] = {'simplevar': {'job': reduceJobId, 'field': 'simplevar'}}
    
    
def mapGatk():
    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'
    
    print "Converting Contigset to Fasta"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['original_contig_set']), shell=True)
    print "Converting Table to SAM"
    subprocess.check_call("dx_mappingsTableToSam --table_id %s --output input.sam --region_index_offset -1 %s" % (job['input']['mappings_table_id'], job['input']['interval']), shell=True)
    print "Converting to BAM"
    subprocess.check_call("samtools view -bS input.sam > input.sorted.bam", shell=True)
    print "Adding Read Groups"
    subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.sorted.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1", shell=True)
    print "Indexing"
    subprocess.check_call("samtools index input.rg.bam", shell=True)
    
    

    writeGenomeDict(job['input']['original_contig_set'], "ref.dict")
    writeReferenceIndex(job['input']['original_contig_set'], "ref.fa.fai")
        
    command = job['input']['command'] + job['input']['interval']
    print command
    subprocess.call(command, shell=True)
    
    command = "dx_vcfToSimplevar --table_id %s --vcf_file output.vcf" % (job['input']['tableId'])
    if job['input']['compress_reference']:
        command += " --compress_reference"
    if job['input']['compress_no_call']:
        command += " --compress_no_call"
    if job['input']['store_full_vcf']:
        command += " --store_full_vcf"
    command += " --extract_header"
    print command
    print "In GATK"
    subprocess.call(command ,shell=True)

def buildCommand(job):
    
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.rg.bam -o output.vcf "
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


def writeGenomeDict(contig_set, dictFileName):
    print contig_set
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']

    dictFile = open(dictFileName, 'w')
    dictFile.write("@HD\tVN:1.0\tSN:unsorted\n")
    for i in range(len(sizes)):
        line = "@SQ\tSN:%s\tLN:%s\n" % (names[i], str(sizes[i]))
        dictFile.write(line)
    return

def writeReferenceIndex(contig_set, indexFileName):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    
    dictFile = open(indexFileName, 'w')
    priorOffset = 0
    priorLength = 0
    first = True
    for i in range(len(names)):
        offset = priorOffset + len(names[i])+ 2 + priorLength + math.ceil(float(priorLength)/50.0)
        if first:
            first = False
        else:
            offset += 1
        dictFile.write(names[i]+"\t"+str(sizes[i]) + "\t" + str(int(offset)) + "\t50\t51\n")
        priorOffset = offset
        priorLength = sizes[i]
    return

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
                    return header
        except StopIteration:
            break
    
def reduceGatk():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    t.close(block=True)
    print "Closing Table"
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())

def splitGenomeLength(contig_set, includeInterval, excludeInterval, chunkSize, splits):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

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
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, sizes[chromosome])
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, position+chunkSize)
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

    
    
    
    
    
    
    
    
    