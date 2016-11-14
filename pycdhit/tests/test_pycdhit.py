import unittest
#from pycdhit import cdhit_read
#from pycdhit import cdhit_cluster
from pycdhit import *




class TestCdhit(unittest.TestCase):
    
    def setUp(self):
        pass
        #add some files in future :)
    
    
    def test_read(self):
        c = cdhit_read()
        self.assertEqual(None, c.name)
        self.assertEqual(None, c.pb)
        self.assertEqual(None, c.length)
        c2 = cdhit_read('test_name', 90, '431aa')
        self.assertEqual('test_name', c2.name)
        self.assertEqual(90, c2.pb)
        self.assertEqual('431aa', c2.length)
        self.assertEqual({'len':'431aa', 'name':'test_name', 'pb':90}, c2.to_json())
        c.name = 'test_name'
        c.pb = 90
        c.length = '431aa'
        #self.assertEqual(c,c2)
        
    def test_cluster(self):
        c = cdhit_cluster('some name')
        self.assertEqual('some name', c.name)
        self.assertEqual(None, c.gene_name)
        self.assertEqual([], c.reads)
        c2 = cdhit_cluster('some name', [cdhit_read('g_name', '*', '435aa'), cdhit_read('test_name', 90, '431aa')], 'gname')
        self.assertEqual('some name', c2.name)
        self.assertEqual('gname', c2.gene_name)
        #self.assertEqual([cdhit_read('g_name', '*', '435aa'), cdhit_read('test_name', 90, '431aa')], c2.reads)
        e_json = {'gene':'gname', 'cluster':'some name', 'reads':[{'len':'435aa', 'name':'g_name', 'pb':'*'},{'len':'431aa', 'name':'test_name', 'pb':90}]}
        self.assertEqual(e_json, c2.to_json())
        self.assertEqual([],c.reads)
        c.append(cdhit_read('g_name', '*', '435aa'))
        self.assertEqual([cdhit_read('g_name', '*', '435aa')],c.reads)
        self.assertEqual(1,len(c.reads))
        self.assertEqual({'len':'435aa', 'name':'g_name', 'pb':'*'}, c.reads[0].to_json())
        
        
    def test_result(self):
        cr = cdhit_result()
        self.assertEqual(None, cr.name)
        self.assertEqual([], cr.data)
        self.assertEqual({'sample':None, 'data':[]})
        
    
    def test_set(self):
        pass
    
if __name__ == "__main__":
    unittest.main()