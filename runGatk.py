#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of gatk_unifiedgenotyper (DNAnexus platform app).
#
# (The MIT Expat License)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

import dxpy
import subprocess, logging
import os, sys, re, math, operator

from multiprocessing import Pool, cpu_count

@dxpy.entry_point('main')
def main():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'

    if job['input']['output_mode'] == "EMIT_VARIANTS_ONLY":
        job['input']['infer_no_call'] = False

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings'][0]['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()

    #This controls the degree of parallelism in GATK
    reads = 0
    for x in job['input']['mappings']:
        reads += int(dxpy.DXGTable(x).describe()['length'])
    chunks = int(reads/job['input']['reads_per_job'])+1

    command = buildCommand(job)

    #callVariantsOnSample(mappingsTable, command)

    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    variants_schema = [
        {"name": "chr", "type": "string"},
        {"name": "lo", "type": "int32"},
        {"name": "hi", "type": "int32"},
        {"name": "ref", "type": "string"},
        {"name": "alt", "type": "string"},
        {"name": "qual", "type": "double"},
        {"name": "ids", "type": "string"}
         ]

    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    headerInfo = extractHeader("/tmp/header.txt", elevatedTags)
    description = {}
    samples = []

    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]

    formats = {}
    infos = {}
    filters = {}

    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}

    samples = []
    for i in range(len(job['input']['mappings'])):
        samples.append(dxpy.DXGTable(job['input']['mappings'][i]).describe()['name'].replace(" ", ""))
    numSamples = len(job['input']['mappings'])
    if job['input']['call_multiple_samples'] == False:
        numSamples = 1
        samples = ["Sample_0"]
    #For each sample, write the sample-specific columns
    for i in range(numSamples):
      variants_schema.extend([
        {"name": "genotype_"+str(i), "type": "string"},
        {"name": "phasing_"+str(i), "type": "string"},
        {"name": "type_"+str(i), "type": "string"},
        {"name": "variation_qual_"+str(i), "type": "double"},
        {"name": "genotype_qual_"+str(i), "type": "double"},
        {"name": "coverage_"+str(i), "type": "string"},
        {"name": "total_coverage_"+str(i), "type": "int32"}
      ])
      indices.append(dxpy.DXGTable.lexicographic_index([["type_"+str(i), "ASC"]], 'type_'+str(i)))
      for k, v in headerInfo['tags']['format'].iteritems():
        if "format_"+k not in elevatedTags:
          variants_schema.append({"name": "format_"+k+"_"+str(i), "type":translateTagTypeToColumnType(v)})

    variantsTable = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr", "lo", "hi", "gri")])
    tableId = variantsTable.get_id()
    variantsTable = dxpy.open_dxgtable(tableId)
    variantsTable.add_types(["Variants", "gri"])

    details = {'samples':samples, 'original_contigset':job['input']['reference'], 'original_mappings':job['input']['mappings'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info']}
    #if headerInfo.get('filters') != {}:
    #  details['filters'] = headerInfo['filters']
    variantsTable.set_details(details)

    if 'output_name' in job['input']:
        variantsTable.rename(job['input']['output_name'])
    elif (job['input']['genotype_likelihood_model'] == "SNP"):
        variantsTable.rename(mappingsTable.describe()['name'] + " SNP calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "INDEL"):
        variantsTable.rename(mappingsTable.describe()['name'] + " indel calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "BOTH"):
        variantsTable.rename(mappingsTable.describe()['name'] + " SNP and indel calls by GATK")
    else:
        variantsTable.rename(mappingsTable.describe()['name'] + " variant calls by GATK")

    reduceInput = {}
    #commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
    commandList = splitGenomeLengthLargePieces(originalContigSet, chunks)

    for i in range(len(commandList)):
        if len(commandList[i]) > 0:
            mapInput = {
                'mappings_tables': job['input']['mappings'],
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'command': buildCommand(job),
                'compress_reference': job['input']['compress_reference'],
                'infer_no_call': job['input']['infer_no_call'],
                'compress_no_call': job['input']['compress_no_call'],
                'intervals_to_include': job['input'].get('intervals_to_process'),
                'intervals_to_exclude': job['input'].get('intervals_to_exclude'),
                'intervals_merging': job['input']['intervals_merging'],
                'part_number': i,
                'samples': samples,
                'call_multiple_samples': job['input']['call_multiple_samples']
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapGatk").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()

    job['output'] = {'variants': {'job': reduceJobId, 'field': 'variants'}}

def runAndCatchGATKError(command, shell=True):
    # Added to capture any errors outputted by GATK
    try:
        subprocess.check_output(command, stderr=subprocess.STDOUT, shell=shell)
    except subprocess.CalledProcessError, e:
        print e 
        error = '\n'.join([l for l in e.output.splitlines() if l.startswith('##### ERROR MESSAGE:')])
        if error: 
            raise dxpy.AppError("App failed with GATK error. Please see logs for more information: {err}".format(err=error))
        else: 
            raise dxpy.AppInternalError("App failed with error. Please see logs for more information: {err}".format(err=e))         


def mapGatk():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:opt/jar/CreateSequenceDictionary.jar'
    print os.environ

    regionFile = open("regions.txt", 'w')
    regionFile.write(job['input']['interval'])

    regionFile.close()

    if job['input']['intervals_merging'] != "INTERSECTION" and job["input"].get("intervals_to_include") != None and job["input"].get("intervals_to_include") != "":
        job['input']['interval'] = splitUserInputRegions(job['input']['interval'], job['input']['intervals_to_include'], "-L")
        if job['input']['interval'] == '':
            job['output']['id'] = job['input']['tableId']
            return
        
    gatkIntervals = open("regions.interval_list", 'w')
    for x in re.findall("-L ([^:]*):(\d+)-(\d+)", job['input']['interval']):
        gatkIntervals.write(x[0] + ":" + x[1] + "-" + x[2] + "\n")
    gatkIntervals.close()

    print "Converting Contigset to Fasta"
    subprocess.check_call("dx-contigset-to-fasta %s ref.fa" % (job['input']['original_contig_set']), shell=True)

    for i in range(len(job['input']['mappings_tables'])):
        mappingsTableId = dxpy.DXGTable(job['input']['mappings_tables'][i]).get_id()

        if job['input']['call_multiple_samples'] == False:
            sample = "Sample_0"
        else:
            sample = job['input']['samples'][i]

        print "Converting Table to SAM"
        #subprocess.check_call("dx-mappings-to-sam %s --output input.sam --region_index_offset -1 --region_file regions.txt --sample %s" % (job['input']['mappings_table_id']), shell=True)
        subprocess.check_call("dx_mappings_to_sam.py %s --output input.%d.sam --region_index_offset -1 --region_file regions.txt --sample %s" % (mappingsTableId, i, sample), shell=True)
        if checkSamContainsRead("input.%d.sam" % i):
            print "Converting to BAM"
            subprocess.check_call("samtools view -bS input.%d.sam > input.%d.bam" % (i, i), shell=True)
        else:
            subprocess.check_call("samtools view -HbS input.%d.sam > input.%d.bam" % (i, i), shell=True)
        print "Sorting BAM"
        subprocess.check_call("samtools sort input.%d.bam input.%d.sorted" % (i, i), shell=True)
        print "Indexing"
        subprocess.check_call("samtools index input.%d.sorted.bam" % i, shell=True)
        job['input']['command'] += " -I input.%d.sorted.bam" % i
        
    print "Indexing Reference"
    subprocess.check_call("samtools faidx ref.fa", shell=True)
    runAndCatchGATKError("java -Xmx4g net.sf.picard.sam.CreateSequenceDictionary REFERENCE=ref.fa OUTPUT=ref.dict" ,shell=True)

    command = job['input']['command'] + job['input']['interval']

    if 'quality' not in dxpy.DXGTable(mappingsTableId).get_col_names():
        print "Quality scores not found in mappings table, adding default quality scores"
        command += " --defaultBaseQualities 20"

    print "In GATK"
    runAndCatchGATKError(command, shell=True)

    command = "dx_vcfToVariants2 --table_id %s --vcf_file output.vcf --region_file regions.txt" % (job['input']['tableId'])
    if job['input']['compress_reference']:
        command += " --compress_reference"
    if job['input']['infer_no_call']:
        command += " --infer_no_call"
    if job['input']['compress_no_call']:
        command += " --compress_no_call"

    print "Parsing Variants"
    subprocess.check_call(command, shell=True)

    job['output']['id'] = job['input']['tableId']
    
def buildCommand(job):

    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -o output.vcf -rf BadCigar"
    if job['input']['output_mode'] != "EMIT_VARIANTS_ONLY":
        command += " -out_mode " + (job['input']['output_mode'])
    if job['input']['call_confidence'] != 30.0:
        command += " -stand_call_conf " +str(job['input']['call_confidence'])
    if job['input']['emit_confidence'] != 30.0:
        command += " -stand_emit_conf " +str(job['input']['emit_confidence'])
    if job['input']['intervals_merging'] == "INTERSECTION":
        if job['input'].get('intervals_to_process') != None:
            command += " " + job['input']['intervals_to_process']
    if job['input'].get('intervals_to_exclude') != None:
        command += " " + job['input']['intervals_to_exclude']
    if job['input']['pcr_error_rate'] != 0.0001:
        command += " -pcr_error " +str(job['input']['pcr_error_rate'])
    if job['input']['heterozygosity'] != 0.001:
        command += " -hets " + str(job['input']['heterozygosity'])
    if job['input']['indel_heterozygosity'] != 0.000125:
        command += " -indelHeterozygosity " + str(job['input']['indel_heterozygosity'])
    if job['input']['genotype_likelihood_model'] != "SNP":
        command += " -glm " + job['input']['genotype_likelihood_model']
    if job['input']['minimum_base_quality'] != 17:
        command += " -mbq " + str(job['input']['minimum_base_quality'])
    if job['input']['max_alternate_alleles'] != 3:
        command += " -maxAlleles " + str(job['input']['max_alternate_alleles'])
    if job['input']['max_deletion_fraction'] != 0.05:
        command += " -deletions " + str(job['input']['max_deletion_fraction'])
    if job['input']['min_indel_count'] != 5:
        command += " -minIndelCnt " + str(job['input']['min_indel_count'])
    if job['input']['non_reference_probability_model'] != "EXACT":
        if job['input']['non_reference_probability_model'] != "GRID_SEARCH":
            raise dxpy.AppError("Option \"Probability Model\" must be either \"EXACT\" or \"GRID_SEARCH\". Found " + job['input']['non_reference_probability_model'] + " instead")
        command += " -pnrm " + str(job['input']['non_reference_probability_model'])

    threads = str(cpu_count())
    if job["input"].get("num_threads") != None:
        if job["input"]["num_threads"] < str(cpu_count) and job["input"]["num_threads"] > 0:
            threads = str(job["input"]["num_threads"])

    command += " --num_threads " + threads
    command += " -L regions.interval_list "

    if job['input']['downsample_to_coverage'] != 250:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    if job['input']['nondeterministic']:
        command += " -ndrs "

    if job['input']['calculate_BAQ'] != "OFF":
        if job['input']['calculate_BAQ'] != "CALCULATE_AS_NECESSARY" and job['input']['calculate_BAQ'] != "RECALCULATE":
            raise dxpy.AppError("Option \"Calculate BAQ\" must be either \"OFF\" or or \"CALCULATE_AS_NECESSARY\" \"RECALCULATE\". Found " + job['input']['calculate_BAQ'] + " instead")
        command += " -baq " + job['input']['calculate_BAQ']
        if job['input']['BAQ_gap_open_penalty'] != 40.0:
            command += " -baqGOP " + str(job['input']['BAQ_gap_open_penalty'])
    if job['input']['no_output_SLOD']:
        command += " -nosl "

    #print command
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

def reduceGatk():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    print "Closing Table"
    t.close()
    job['output']['variants'] = dxpy.dxlink(t.get_id())

def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
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

def splitGenomeLengthLargePieces(contig_set, chunks):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    for i in range(len(names)):
        print names[i]+":"+str(sizes[i])


    commandList = []
    for i in range(chunks):
        commandList.append('')
    position = 0
    chromosome = 0
    chunkSize = sum(sizes) / chunks
    currentChunk = 0
    currentLength = 0

    while chromosome < len(names):
        if position + (chunkSize - currentLength) >= sizes[chromosome]:
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, sizes[chromosome])
            currentLength += sizes[chromosome] - position
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, position+(chunkSize-currentLength)+1)
            position += (chunkSize-currentLength) + 1
            if currentChunk < chunks-1:
                currentChunk += 1
            currentLength = 0
    return commandList

def callVariantsOnSample(mappingsTable, command):
    subprocess.check_call("dx_mappingsTableToSam --table_id %s --output dummy.sam --end_row 100" % (job['input']['mappings_table_id']), shell=True)

def extractHeader(vcfFileName, elevatedTags):
    result = {'columns': '', 'tags' : {'format' : {}, 'info' : {} }, 'filters' : {}}
    for line in open(vcfFileName):
        tag = re.findall("ID=(\w+),", line)
        if len(tag) > 0:
          tagType = ''
          if line.count("FORMAT") > 0:
            tagType = 'format'
          elif line.count("INFO") > 0:
            tagType = 'info'
          elif line.count("FILTER") > 0:
            result['filters'][re.findall("ID=(\w+),")[0]] = re.findall('Description="(.*)"')[0]

          typ = re.findall("Type=(\w+),", line)
          if tagType != '':
            number = re.findall("Number=(\w+)", line)
            description = re.findall('Description="(.*)"', line)
            if len(number) == 0:
              number = ['.']
            if len(description) == 0:
              description = ['']
            if "format_"+tag[0] not in elevatedTags:
                result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
        if line[0] == "#" and line[1] != "#":
          result['columns'] = line.strip()
        if line == '' or line[0] != "#":
            break
    return result

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False

def translateTagTypeToColumnType(tag):
  if tag['type'] == "Flag":
    return "boolean"
  if tag['number'] != '1':
    return 'string'
  if tag['type'] == "Integer":
    return 'int32'
  if tag['type'] == "Float":
    return "double"
  return "string"

def splitUserInputRegions(jobRegions, inputRegions, prefix):
    
    jobList = re.findall("-L ([^:]*):(\d+)-(\d+)", jobRegions)    
    inputList = re.findall("-L ([^:]*):(\d+)-(\d+)", inputRegions)
    
    result = ""
    for x in inputList:
        for y in jobList:
            if(x[0] == y[0]):
                lo = max(int(x[1]), int(y[1]))
                hi = min(int(x[2]), int(y[2]))
                if hi > lo:
                    result += " %s %s:%d-%d" % (prefix, x[0], lo, hi)
                    
    return result
