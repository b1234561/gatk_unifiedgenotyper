#!/usr/bin/env python

import dxpy
import subprocess, logging
import os, sys, re, math, operator

def main():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()
    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    variants_schema = [{"name": "chr", "type": "string"},
                       {"name": "lo", "type": "int32"},
                       {"name": "hi", "type": "int32"},
                       {"name": "type", "type": "string"},    # TODO: change this type to uint once there is an abstraction method for enum?
                       {"name": "ref", "type": "string"},
                       {"name": "alt", "type": "string"},
                       {"name": "qual", "type": "int32"},
                       {"name": "coverage", "type": "int32"},
                       {"name": "genotype_quality", "type": "int32"}]

    if job['input']['store_full_vcf']:
        variants_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])

    simpleVar = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr", "lo", "hi", "gri")])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details({'original_contigset': originalContigSet})
    simpleVar.add_types(["SimpleVar", "gri"])

    if 'output name' in job['input']:
        simpleVar.rename(job['input']['output name'])
    elif (job['input']['genotype_likelihood_model'] == "SNP"):
        simpleVar.rename(mappingsTable.describe()['name'] + " SNP calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "INDEL"):
        simpleVar.rename(mappingsTable.describe()['name'] + " indel calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "BOTH"):
        simpleVar.rename(mappingsTable.describe()['name'] + " SNP and indel calls by GATK")
    else:
        simpleVar.rename(mappingsTable.describe()['name'] + " variant calls by GATK")

    reduceInput = {}
    #commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
    commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['maximum_chunks'])

    for i in range(len(commandList)):
        if len(commandList[i]) > 0:
            mapInput = {
                'mappings_table_id': mappingsTableId,
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'command': buildCommand(job),
                'compress_reference': job['input']['compress_reference'],
                'compress_no_call': job['input']['compress_no_call'],
                'store_full_vcf': job['input']['store_full_vcf'],
                'part_number': i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapGatk").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()

    job['output'] = {'simplevar': {'job': reduceJobId, 'field': 'simplevar'}}

def mapGatk():

    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:opt/jar/CreateSequenceDictionary.jar'
    

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
    subprocess.check_call("dx_mappingsTableToSam --table_id %s --output input.sam --region_index_offset -1 --region_file regions.txt" % (job['input']['mappings_table_id']), shell=True)

    if checkSamContainsRead("input.sam"):
        print "Converting to BAM"
        subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
        #print "Adding Read Groups"
        #subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.sorted.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1", shell=True)
        print "Indexing"
        subprocess.check_call("samtools index input.bam", shell=True)

    
        #subprocess.check_call("dx_writeReferenceIndex --contig_set %s --writeSamtoolsIndex ref.fa.fai --writePicardDictionary ref.dict" % (job['input']['original_contig_set']), shell=True)

        print "Indexing Dictionary"
        subprocess.check_call("samtools faidx ref.fa", shell=True)
        subprocess.call("java -Xmx4g net.sf.picard.sam.CreateSequenceDictionary REFERENCE=ref.fa OUTPUT=ref.dict" ,shell=True)


        command = job['input']['command'] + job['input']['interval']
        print command
        subprocess.call(command, shell=True)
        
        vcf = dxpy.dxlink(dxpy.upload_local_file("output.vcf")) 
        #print open("output.vcf", 'r').read()
        
        command = "dx_vcfToSimplevar --table_id %s --vcf_file output.vcf" % (job['input']['tableId'])
        if job['input']['compress_reference']:
            command += " --compress_reference"
        if job['input']['compress_no_call']:
            command += " --compress_no_call"
        if job['input']['store_full_vcf']:
            command += " --store_full_vcf"
        command += " --extract_header"
        print "In GATK"

        subprocess.call(command, shell=True)

    else:
        print "No reads in SAM"
        print job['input']['interval']

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
    command += " -L regions.interval_list"

    if job['input']['downsample_to_coverage'] != 50000:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    command += " -dt " + job['input']['downsampling_type']
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

def extractHeader(vcfFile):
    header = ''
    fileIter = vcfFile.__iter__()

    # Additional data will contain the extra format and info columns that
    # are optional in VCF and may not be present in the VCF file. These are
    # stored in an extended table.
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
    print "Closing Table"
    t.close()
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())

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

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False
