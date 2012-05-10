#!/usr/bin/env python
import os, sys, unittest, json, subprocess

import dxpy, dxpy.program_builder
from dxpy.exceptions import *

import subprocess

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")


def makeInputs():
    try:
        contigset_importer = dxpy.DXProgram(dxpy.find_data_objects(classname="program", properties={"name": "fasta_contigset_importer"}).next()['id'])
    except StopIteration:
        raise Exception("fasta_contigset_importer or LetterSpaceFileObjectToReadsTable not found, please upload them")

    sam = dxpy.upload_local_file(os.path.join(test_resources_dir, "reads.sam"), wait_on_close=True)
        
    genome_archive = dxpy.upload_local_file(os.path.join(test_resources_dir, "hg19_chr22.fa.xz"), wait_on_close=True)
    contigset_importer_input = {"name": "hg19_chr22", "sequence_file": dxpy.dxlink(genome_archive)}
    print "Running fasta_contigset_importer with", contigset_importer_input
    job = contigset_importer.run(contigset_importer_input)
    job.wait_on_done()
    contig_set = job.describe()["output"]["contig_set"]
    print contig_set

    return { 'mappings':{"$dnanexus_link": "gtable-9ybB2Y80000FKXKJB8yQ0009"}, "output_mode":"EMIT_VARIANTS_ONLY"}


class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_input = makeInputs()
        bundled_resources = dxpy.program_builder.upload_resources(src_dir)
        cls.program_id = dxpy.program_builder.upload_program(src_dir, bundled_resources, overwrite=True)
    
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_gatk(self):
        program_input = self.base_input
        
        job = dxpy.DXProgram(self.program_id).run(program_input)
        
        print "Waiting for job to complete"
        job.wait_on_done()
        print json.dumps(job.describe()["output"])


if __name__ == '__main__':
    unittest.main()
