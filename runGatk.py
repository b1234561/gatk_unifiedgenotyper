#!/usr/bin/env python

import dxpy
import subprocess, logging
import os, sys, re, math, operator

from multiprocessing import Pool, cpu_count

def main():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'

    if job['input']['output_mode'] == "EMIT_VARIANTS_ONLY":
        job['input']['infer_no_call'] = False

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()
    
    #This controls the degree of parallelism in GATK
    chunks = int(mappingsTable.describe()['length']/job['input']['reads_per_job'])+1
    
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
        {"name": "qual", "type": "float"},
        {"name": "filter", "type": "string"},
        {"name": "ids", "type": "string"}
         ]

    headerInfo = extractHeader("/tmp/header.txt")
    description = {}
    samples = []

    print headerInfo
    print headerInfo['tags']['format']

    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]
    
    formats = {}
    infos = {}
    filters = {}
    
    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}
    
    numSamples = 1
    #For each sample, write the sample-specific columns
    for i in range(numSamples):
      variants_schema.extend([
        {"name": "genotype_"+str(i), "type": "string"},
        {"name": "phasing_"+str(i), "type": "string"},
        {"name": "type_"+str(i), "type": "string"},
        {"name": "variation_qual_"+str(i), "type": "float"},
        {"name": "genotype_qual_"+str(i), "type": "float"},
        {"name": "coverage_"+str(i), "type": "string"},
        {"name": "total_coverage_"+str(i), "type": "int32"}
      ])
      indices.append(dxpy.DXGTable.lexicographic_index([["type_"+str(i), "ASC"]], 'type_'+str(i)))
      samples.append("Sample_0")
      for k, v in headerInfo['tags']['format'].iteritems():
        if "format_"+k not in elevatedTags:
          variants_schema.append({"name": "format_"+k+"_"+str(i), "type":translateTagTypeToColumnType(v)})

    variantsTable = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr", "lo", "hi", "gri")])
    tableId = variantsTable.get_id()
    variantsTable = dxpy.open_dxgtable(tableId)
    variantsTable.add_types(["Variants", "gri"])
    
    details = {'samples':samples, 'original_contigset':job['input']['reference'], 'original_mappings':job['input']['mappings'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info']}
    if headerInfo.get('filters') != {}:
      details['filters'] = headerInfo['filters']
    variantsTable.set_details(details)

    if 'output name' in job['input']:
        variantsTable.rename(job['input']['output name'])
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
                'mappings_table_id': mappingsTableId,
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'command': buildCommand(job),
                'compress_reference': job['input']['compress_reference'],
                'infer_no_call': job['input']['infer_no_call'],
                'compress_no_call': job['input']['compress_no_call'],
                'part_number': i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapGatk").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()

    job['output'] = {'variants': {'job': reduceJobId, 'field': 'variants'}}

def mapGatk():

    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:opt/jar/CreateSequenceDictionary.jar'
    print os.environ
    
    regionFile = open("regions.txt", 'w')
    regionFile.write(job['input']['interval'])
    regionFile.close()

    gatkIntervals = open("regions.interval_list", 'w')
    for x in re.findall("(\w+):(\d+)-(\d+)", job['input']['interval']):
        gatkIntervals.write(x[0] + ":" + x[1] + "-" + x[2] + "\n")
    gatkIntervals.close()

    print "Converting Contigset to Fasta"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['original_contig_set']), shell=True)

    print "Converting Table to SAM"
    subprocess.check_call("dx_mappingsToSam --table_id %s --output input.sam --region_index_offset -1 --region_file regions.txt" % (job['input']['mappings_table_id']), shell=True)

    if checkSamContainsRead("input.sam"):
        print "Converting to BAM"
        subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
        print "Indexing"
        subprocess.check_call("samtools index input.bam", shell=True)
        print "Indexing Dictionary"
        subprocess.check_call("samtools faidx ref.fa", shell=True)
        subprocess.call("java -Xmx4g net.sf.picard.sam.CreateSequenceDictionary REFERENCE=ref.fa OUTPUT=ref.dict" ,shell=True)

        command = job['input']['command'] + job['input']['interval']
        #print command
        print "In GATK"
        subprocess.call(command, shell=True)
        #command += " | "
    
        command = "dx_vcfToVariants2 --table_id %s --vcf_file output.vcf --region_file regions.txt" % (job['input']['tableId'])
        if job['input']['compress_reference']:
            command += " --compress_reference"
        if job['input']['infer_no_call']:
            command += " --infer_no_call"
        if job['input']['compress_no_call']:
            command += " --compress_no_call"
        
        print "Parsing Variants"
        subprocess.call(command, shell=True)

def buildCommand(job):
    
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.bam -o output.vcf "
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
    command += " --num_threads " + str(cpu_count())
    command += " -L regions.interval_list"

    if job['input']['downsample_to_coverage'] != 50000:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    if job['input']['nondeterministic']:
        command += " -ndrs "
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
    t.close(block=True)
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

def extractHeader(vcfFileName):
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
    return "float"
  return "string"
