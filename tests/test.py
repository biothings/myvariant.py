import unittest
import sys
import os
import types
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from pandas import DataFrame
    pandas_avail = True
except ImportError:
    pandas_avail = False
try:
    import requests_cache
    caching_avail = True
except ImportError:
    caching_avail = False
sys.path.insert(0, os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])
import myvariant
sys.stderr.write('"myvariant {0}" loaded from "{1}"\n'.format(myvariant.__version__, myvariant.__file__))


class TestMyVariantPy(unittest.TestCase):

    def setUp(self):
        self.mv = myvariant.MyVariantInfo()
        self.query_list1 = [
            'chr1:g.866422C>T',
            'chr1:g.876664G>A',
            'chr1:g.69635G>C',
            'chr1:g.69869T>A',
            'chr1:g.881918G>A',
            'chr1:g.865625G>A',
            'chr1:g.69892T>C',
            'chr1:g.879381C>T',
            'chr1:g.878330C>G'
        ]
        self.query_list2 = [
            'rs374802787',
            'rs1433078',
            'rs1433115',
            'rs377266517',
            'rs587640013',
            'rs137857980',
            'rs199710579',
            'rs186823979',
            'rs2276240',
            'rs372452565'
        ]

    def test_format_hgvs(self):
        self.assertEqual(myvariant.format_hgvs("1", 35366, "C", "T"),
                         'chr1:g.35366C>T')
        self.assertEqual(myvariant.format_hgvs("chr2", 17142, "G", "GA"),
                         'chr2:g.17142_17143insA')
        self.assertEqual(myvariant.format_hgvs("1", 10019, "TA", "T"),
                         'chr1:g.10020del')
        self.assertEqual(myvariant.format_hgvs("MT", 8270, "CACCCCCTCT", "C"),
                         'chrMT:g.8271_8279del')
        self.assertEqual(myvariant.format_hgvs("7", 15903, "G", "GC"),
                         'chr7:g.15903_15904insC')
        self.assertEqual(myvariant.format_hgvs("X", 107930849, "GGA", "C"),
                         'chrX:g.107930849_107930851delinsC')
        self.assertEqual(myvariant.format_hgvs("20", 1234567, "GTC", "GTCT"),
                         'chr20:g.1234569_1234570insT')

    def test_metadata(self):
        meta = self.mv.metadata()
        self.assertTrue("stats" in meta)
        self.assertTrue("total" in meta['stats'])

    def test_getvariant(self):
        v = self.mv.getvariant("chr9:g.107620835G>A")
        self.assertEqual(v['_id'], "chr9:g.107620835G>A")
        self.assertEqual(v['snpeff']['ann']['genename'], 'ABCA1')

        v = self.mv.getvariant("'chr1:g.1A>C'")   # something does not exist
        self.assertEqual(v, None)

    def test_getvariant_with_fields(self):
        v = self.mv.getvariant("chr9:g.107620835G>A", fields="dbnsfp,cadd,cosmic")
        self.assertTrue('_id' in v)
        self.assertTrue('dbnsfp' in v)
        self.assertTrue('cadd' in v)
        self.assertTrue('cosmic' in v)

    def test_getvariants(self):
        v_li = self.mv.getvariants(self.query_list1)
        self.assertEqual(len(v_li), 9)
        self.assertEqual(v_li[0]['_id'], self.query_list1[0])
        self.assertEqual(v_li[1]['_id'], self.query_list1[1])
        self.assertEqual(v_li[2]['_id'], self.query_list1[2])

        self.mv.step = 4
        # test input is a string of comma-separated ids
        v_li2 = self.mv.getvariants(','.join(self.query_list1))
        self.assertEqual(v_li, v_li2)
        # test input is a tuple
        v_li2 = self.mv.getvariants(tuple(self.query_list1))
        self.assertEqual(v_li, v_li2)

        # test input is a generator
        def _input(li):
            for x in li:
                yield x
        v_li2 = self.mv.getvariants(_input(self.query_list1))
        self.assertEqual(v_li, v_li2)
        self.mv.step = 1000

    def test_query(self):
        qres = self.mv.query('dbnsfp.genename:cdk2', size=5)
        self.assertTrue('hits' in qres)
        self.assertEqual(len(qres['hits']), 5)

    def test_query_hgvs(self):
        qres = self.mv.query('"NM_000048.3:c.566A>G"', size=5)
        self.assertTrue('hits' in qres)
        self.assertEqual(len(qres['hits']), 1)

    def test_query_rsid(self):
        qres = self.mv.query('dbsnp.rsid:rs58991260')
        self.assertTrue('hits' in qres)
        self.assertEqual(len(qres['hits']), 1)
        self.assertEqual(qres['hits'][0]['_id'], 'chr1:g.218631822G>A')
        qres2 = self.mv.query('rs58991260')
        # exclude _score field before comparison
        qres['hits'][0].pop('_score')
        qres2['hits'][0].pop('_score')
        self.assertEqual(qres['hits'], qres2['hits'])

    def test_query_symbol(self):
        qres = self.mv.query('snpeff.ann.genename:cdk2')
        self.assertTrue('hits' in qres)
        self.assertTrue(qres['total'] > 5000)
        self.assertEqual(qres['hits'][0]['snpeff']['ann'][0]['genename'], 'CDK2')

    def test_query_genomic_range(self):
        qres = self.mv.query('chr1:69000-70000')
        self.assertTrue('hits' in qres)
        self.assertTrue(qres['total'] >= 3)

    def test_query_fetch_all(self):
        qres = self.mv.query('chr1:69500-70000', fields="chrom")
        total = qres['total']

        qres = self.mv.query('chr1:69500-70000', fields="chrom", fetch_all=True)
        self.assertTrue(isinstance(qres, types.GeneratorType))
        self.assertEqual(total, len(list(qres)))

    def test_querymany(self):
        qres = self.mv.querymany(self.query_list1, verbose=False)
        self.assertEqual(len(qres), 9)

        self.mv.step = 4
        # test input as a string
        qres2 = self.mv.querymany(','.join(self.query_list1), verbose=False)
        self.assertEqual(qres, qres2)
        # test input as a tuple
        qres2 = self.mv.querymany(tuple(self.query_list1), verbose=False)
        self.assertEqual(qres, qres2)
        # test input as a iterator
        qres2 = self.mv.querymany(iter(self.query_list1), verbose=False)
        self.assertEqual(qres, qres2)
        self.mv.step = 1000

    def test_querymany_with_scopes(self):
        qres = self.mv.querymany(['rs58991260', 'rs2500'], scopes='dbsnp.rsid', verbose=False)
        self.assertEqual(len(qres), 2)

        qres = self.mv.querymany(['RCV000083620', 'RCV000083611', 'RCV000083584'], scopes='clinvar.rcv_accession', verbose=False)
        self.assertEqual(len(qres), 3)

        qres = self.mv.querymany(['rs2500', 'RCV000083611', 'COSM1392449'],
                                 scopes='clinvar.rcv_accession,dbsnp.rsid,cosmic.cosmic_id', verbose=False)
        self.assertEqual(len(qres), 3)

    def test_querymany_fields(self):
        ids = ['COSM1362966', 'COSM990046', 'COSM1392449']
        qres1 = self.mv.querymany(ids, scopes='cosmic.cosmic_id', fields=['cosmic.tumor_site', 'cosmic.cosmic_id'], verbose=False)
        self.assertEqual(len(qres1), 3)

        qres2 = self.mv.querymany(ids, scopes='cosmic.cosmic_id', fields='cosmic.tumor_site,cosmic.cosmic_id', verbose=False)
        self.assertEqual(len(qres2), 3)

        self.assertEqual(qres1, qres2)

    def test_querymany_notfound(self):
        qres = self.mv.querymany(['rs58991260', 'rs2500', 'NA_TEST'], scopes='dbsnp.rsid', verbose=False)
        self.assertEqual(len(qres), 3)
        self.assertEqual(qres[2], {"query": 'NA_TEST', "notfound": True})

    def test_querymany_dataframe(self):
        if not pandas_avail:
            from nose.plugins.skip import SkipTest
            raise SkipTest
        qres = self.mv.querymany(self.query_list2, scopes='dbsnp.rsid', fields='dbsnp', as_dataframe=True, verbose=False)
        self.assertTrue(isinstance(qres, DataFrame))
        self.assertTrue('dbsnp.vartype' in qres.columns)
        self.assertEqual(set(self.query_list2), set(qres.index))

    def test_querymany_step(self):
        qres1 = self.mv.querymany(self.query_list2, scopes='dbsnp.rsid', verbose=False)
        default_step = self.mv.step
        self.mv.step = 3
        qres2 = self.mv.querymany(self.query_list2, scopes='dbsnp.rsid', verbose=False)
        self.mv.step = default_step
        self.assertEqual(qres1, qres2)

    def test_get_fields(self):
        fields = self.mv.get_fields()
        self.assertTrue('dbsnp' in fields.keys())
        self.assertTrue('clinvar' in fields.keys())

    def test_caching(self):
        if not caching_avail:
            from nose.plugins.skip import SkipTest
            raise SkipTest

        def _getvariant():
            return self.mv.getvariant("chr9:g.107620835G>A")

        def _getvariants():
            return self.mv.getvariants(["chr9:g.107620835G>A", "chr1:g.876664G>A"])

        def _query():
            return self.mv.query("dbsnp.rsid:rs58991260")

        def _querymany():
            return self.mv.querymany(['rs58991260', 'rs2500'], scopes='dbsnp.rsid')

        def _cache_request(f):
            current_stdout = sys.stdout
            try:
                out = StringIO()
                sys.stdout = out
                r = f()
                output = out.getvalue().strip()              
            finally:
                sys.stdout = current_stdout

            return ('Returning cached result' in output, r)

        from_cache, pre_cache_r = _cache_request(_getvariant)
        self.assertFalse(from_cache)
        
        self.mv.set_caching('mvc')
        
        # populate cache
        from_cache, cache_fill_r = _cache_request(_getvariant)
        self.assertTrue(os.path.exists('mvc.sqlite'))
        self.assertFalse(from_cache)
        # is it from the cache?
        from_cache, cached_r = _cache_request(_getvariant)
        self.assertTrue(from_cache)
        
        self.mv.stop_caching()
        # same query should be live - not cached
        from_cache, post_cache_r = _cache_request(_getvariant)
        self.assertFalse(from_cache)

        self.mv.set_caching('mvc')
        # same query should still be sourced from cache
        from_cache, recached_r = _cache_request(_getvariant)
        self.assertTrue(from_cache)

        self.mv.clear_cache()
        # cache was cleared, same query should be live
        from_cache, clear_cached_r = _cache_request(_getvariant)
        self.assertFalse(from_cache)

        # all requests should be identical
        self.assertTrue(all([x == pre_cache_r for x in 
            [pre_cache_r, cache_fill_r, cached_r, post_cache_r, recached_r, clear_cached_r]]))

        # test getvariants POST caching
        from_cache, first_getvariants_r = _cache_request(_getvariants)
        self.assertFalse(from_cache)
        # should be from cache this time
        from_cache, second_getvariants_r = _cache_request(_getvariants)
        self.assertTrue(from_cache)

        # test query GET caching
        from_cache, first_query_r = _cache_request(_query)
        self.assertFalse(from_cache)
        # should be from cache this time
        from_cache, second_query_r = _cache_request(_query)
        self.assertTrue(from_cache)

        # test querymany POST caching
        from_cache, first_querymany_r = _cache_request(_querymany)
        self.assertFalse(from_cache)
        # should be from cache this time
        from_cache, second_querymany_r = _cache_request(_querymany)
        self.assertTrue(from_cache)

        self.mv.stop_caching()

        os.remove('mvc.sqlite')

if __name__ == '__main__':
    unittest.main()
