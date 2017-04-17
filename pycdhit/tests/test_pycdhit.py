import unittest
#from pycdhit import cdhit_read
#from pycdhit import cdhit_cluster
from pycdhit import *
import os




class TestCdhit(unittest.TestCase):
    
    def setUp(self):
        with open('test-1.clstr', 'w') as f:
            f.write(">Cluster 0\n"
                    "0	900aa, >test3... *\n"
                    ">Cluster 1\n"
                    "0	800aa, >test1... *\n"
                    "1	800aa, >s1-test1... at 90.00%\n"
                    ">Cluster 2\n"
                    "0	700aa, >test5... *\n"
                    "1	700aa, >s1-test5... at 73.00%\n"
                    ">Cluster 3\n"
                    "0	500aa, >test2... *\n"
                    "1	500aa, >s1-test2... at 82.00%\n"
                    ">Cluster 4\n"
                    "0	100aa, >test4... *\n"
                    "1	100aa, >s1-test4... at 73.00%")
    
    def tearDown(self):
        os.remove('test-1.clstr')
    
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

        cr.append(cdhit_cluster('Cluster 1', [cdhit_read('test1', '*', '800aa'), cdhit_read('s1-test1', '90.00', '800aa')], 'test1'))
        cr.append(cdhit_cluster('Cluster 2', [cdhit_read('test5', '*', '700aa'), cdhit_read('s1-test5', '73.00', '700aa')], 'test5'))
        cr.append(cdhit_cluster('Cluster 3', [cdhit_read('test2', '*', '500aa'), cdhit_read('s1-test2', '82.00', '500aa')], 'test2'))
        cr.append(cdhit_cluster('Cluster 4', [cdhit_read('test4', '*', '100aa'), cdhit_read('s1-test4', '73.00', '100aa')], 'test4'))
        tcr = cdhit_result()
        tcr.load_from_file('test-1.clstr')
        self.assertEqual(tcr.to_json(), cr.to_json())
        
    def test_result_to_df(self):
        cr = cdhit_result('test-df')
        cr.append(cdhit_cluster('Cluster 1', [cdhit_read('test1', '*', '800aa'), cdhit_read('s1-test1', 90.00, '800aa')], 'test1'))
        cr.append(cdhit_cluster('Cluster 2', [cdhit_read('test5', '*', '700aa'), cdhit_read('s1-test5', 73.00, '700aa')], 'test5'))
        cr.append(cdhit_cluster('Cluster 3', [cdhit_read('test2', '*', '500aa'), cdhit_read('s1-test2', 82.00, '500aa')], 'test2'))
        cr.append(cdhit_cluster('Cluster 4', [cdhit_read('test4', '*', '100aa'), cdhit_read('s1-test4', 73.00, '100aa')], 'test4'))

        df = pd.DataFrame(data={'test1': 90, 'test5': 73, 'test2': 82, 'test4': 73}, index=['test-df'])
        print df
        print '\n----\n'
        print cr.to_df()
        self.assertEqual(df, cr.to_df())
        
    def test_result_get_th_lables(self):
        cr = cdhit_result('test-th')
        cr.append(cdhit_cluster('Cluster 1', [cdhit_read('test1', '*', '800aa'), cdhit_read('s1-test1', 90.00, '800aa')], 'test1'))
        cr.append(cdhit_cluster('Cluster 2', [cdhit_read('test5', '*', '700aa'), cdhit_read('s1-test5', 73.00, '700aa')], 'test5'))
        cr.append(cdhit_cluster('Cluster 3', [cdhit_read('test2', '*', '500aa'), cdhit_read('s1-test2', 82.00, '500aa')], 'test2'))
        cr.append(cdhit_cluster('Cluster 4', [cdhit_read('test4', '*', '100aa'), cdhit_read('s1-test4', 73.00, '100aa')], 'test4'))
        self.assertEqual(['test1', 'test2'],cr.get_thold_labels(80,1))
        self.assertEqual(['test5', 'test4'],cr.get_thold_labels(80,-1))
        cr.append(cdhit_cluster('Cluster 5', [cdhit_read('test6', '*', '100aa')], 'test6'))
        self.assertEqual(['test5', 'test4', 'test6'],cr.get_thold_labels(80,-1))
        self.assertEqual(['test6'],cr.get_thold_labels(0,-1))
        tcr = cdhit_result()
        tcr.load_from_file('test-1.clstr')
        self.assertEqual(['test1', 'test2'],tcr.get_thold_labels(80,1))

        
        
    
    def test_set(self):
        pass
    
if __name__ == "__main__":
    unittest.main()
