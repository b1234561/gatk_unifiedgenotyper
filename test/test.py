#!/usr/bin/env python
import os, sys, unittest, json, subprocess

import dxpy, dxpy.program_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")

class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        bundled_resources = dxpy.program_builder.upload_resources(src_dir)
        cls.program_id = dxpy.program_builder.upload_program(src_dir, bundled_resources, overwrite=True)
    
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_gatk(self):
        sam = dxpy.upload_local_file(os.path.join(test_resources_dir, "reads.sam"), wait_on_close=True)
        reference_sequence = dxpy.upload_local_file(os.path.join(test_resources_dir, "ref.fa"), wait_on_close=True)
        reference_dict = dxpy.upload_local_file(os.path.join(test_resources_dir, "ref.dict"), wait_on_close=True)
        reference_index = dxpy.upload_local_file(os.path.join(test_resources_dir, "ref.fa.fai"), wait_on_close=True)
        
        program_input = {"sam": dxpy.dxlink(sam),
                         "reference_sequence": dxpy.dxlink(reference_sequence),
                         "reference_dict": dxpy.dxlink(reference_dict),
                         "reference_index": dxpy.dxlink(reference_index)
                         }
        
        print program_input
        
        job = dxpy.DXProgram(self.program_id).run(program_input)
        
        print "Waiting for job to complete"
        job.wait_on_done()
        print json.dumps(job.describe()["output"])


if __name__ == '__main__':
    unittest.main()
