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
        self.assertEqual(c,c2)
        
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
        self.assertEqual('gname',c2.get_label())
        self.assertEqual(0, c.get_single_value())
        c.append(cdhit_read('g_name_x', 87, '1435aa'))
        self.assertEqual([cdhit_read('g_name', '*', '435aa'), cdhit_read('g_name_x', 87, '1435aa')],c.reads)
        self.assertEqual(87, c.get_single_value())
        
        
    def test_result(self):
        cr = cdhit_result()
        self.assertEqual(None, cr.name)
        self.assertEqual([], cr.data)
        self.assertEqual({'sample':None, 'data':[]},cr.to_json())
        cr.append(cdhit_cluster('Cluster 0', [cdhit_read('test3', '*', '900aa')], 'test3'))
        #self.assertEqual(cdhit_cluster('Cluster 0', ['test3', '*', '900aa'], 'test3'), cr.data)
        print cr.to_json()
        print cdhit_cluster('Cluster 0', [cdhit_read('test3', '*', '900aa')], 'test3').to_json()
        cr.append(cdhit_cluster('Cluster 1', [cdhit_read('test1', '*', '800aa'), cdhit_read('s1-test1', 90, '800aa')], 'test1'))
        '''
        >Cluster 0
0	900aa, >test3... *
>Cluster 1
0	800aa, >test1... *
1	800aa, >s1-test1... at 90.00%
>Cluster 2
0	700aa, >test5... *
1	700aa, >s1-test5... at 73.00%
>Cluster 3
0	500aa, >test2... *
1	500aa, >s1-test2... at 82.00%
>Cluster 4
0	100aa, >test4... *
1	100aa, >s1-test4... at 73.00%
'''
        
        
    
    def test_set(self):
        pass
    
if __name__ == "__main__":
    unittest.main()
