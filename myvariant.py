# -*- coding: utf-8 -*-
'''
Python Client for MyVariant.Info services
'''
from __future__ import print_function
import sys
import time
import requests
import json
try:
    from pandas import DataFrame
    from pandas.io.json import json_normalize
    df_avail = True
except:
    df_avail = False

__version__ = '2.2.1'

if sys.version_info[0] == 3:
    str_types = str
    from urllib.parse import urlencode
else:
    str_types = (str, unicode)
    from urllib import urlencode


def safe_str(s, encoding='utf-8'):
    '''if input is an unicode string, do proper encoding.'''
    try:
        _s = str(s)
    except UnicodeEncodeError:
        _s = s.encode(encoding)
    return _s


class MyVariantInfo():
    '''This is the client for MyVariant.info web services.
    Example:

        >>> mv = MyVariantInfo()

    '''
    def __init__(self, url='http://myvariant.info/v1'):
        self.url = url
        if self.url[-1] == '/':
            self.url = self.url[:-1]
        self.max_query = 1000
        # delay and step attributes are for batch queries.
        self.delay = 1
        self.step = 1000

    def _dataframe(self, var_obj, dataframe):
        """
        converts gene object to DataFrame (pandas)
        """
        if not df_avail:
            print("Error: pandas module must be installed for as_dataframe option.")
            return
        if dataframe not in [ "by_source", "normal"] :
            raise ValueError("return must be normal or by_source")
        if 'hits' in var_obj:
            if dataframe == "normal":
                df = json_normalize(var_obj['hits'])
            else:
                df = DataFrame.from_dict(var_obj['hits'])
        else:
            if dataframe == "normal":
                df = json_normalize(var_obj)
            else:
                df = DataFrame.from_dict(var_obj)
        return df

    def _get(self, url, params={}):
        debug = params.pop('debug', False)
        return_raw = params.pop('return_raw', False)
        headers = {'user-agent': "Python-requests_myvariant.py/%s (gzip)" % requests.__version__}
        res = requests.get(url, params=params, headers=headers)
        #if debug:
        #    return _url, res, con
        assert res.status_code == 200
        if return_raw:
            return res
        else:
            return res.json()

    def _post(self, url, params):
#        #if debug:
#        #    return url, res, con
        debug = params.pop('debug', False)
        return_raw = params.pop('return_raw', False)
        headers = {'content-type': 'application/x-www-form-urlencoded',
                   'user-agent': "Python-requests_myvariant.py/%s (gzip)" % requests.__version__}
        res = requests.post(url, data=params, headers=headers)
        assert res.status_code == 200
        if return_raw:
            return res
        else:
            return res.json()


    def _format_list(self, a_list, sep=','):
        if isinstance(a_list, (list, tuple)):
            _out = sep.join([safe_str(x) for x in a_list])
        else:
            _out = a_list     # a_list is already a comma separated string
        return _out

    def _repeated_query(self, query_fn, query_li, verbose=True, **fn_kwargs):
        step = min(self.step, self.max_query)
        if len(query_li) <= step:
            # No need to do series of batch queries, turn off verbose output
            verbose = False
        for i in range(0, len(query_li), step):
            is_last_loop = i+step >= len(query_li)
            if verbose:
                print("querying {0}-{1}...".format(i+1, min(i+step, len(query_li))), end="")
            query_result = query_fn(query_li[i:i+step], **fn_kwargs)

            yield query_result

            if verbose:
                print("done.")
            if not is_last_loop and self.delay:
                time.sleep(self.delay)

    @property
    def metadata(self):
        '''Return a dictionary of MyVariant.info metadata.

        Example:

        >>> metadata = mv.metadata

        '''
        _url = self.url+'/metadata'
        return self._get(_url)

    def getvariant(self, geneid, **kwargs):
        '''Return the gene object for the give geneid.
        This is a wrapper for GET query of "/gene/<geneid>" service.

        :param geneid: entrez/ensembl gene id, entrez gene id can be either
                       a string or integer
        :param fields: fields to return, a list or a comma-separated string.
                        If **fields="all"**, all available fields are returned
        :param species: optionally, you can pass comma-separated species names
                        or taxonomy ids
        :param email: optionally, pass your email to help us to track usage
        :param filter: alias for **fields** parameter

        :return: a gene object as a dictionary

        :ref: http://mygene.info/doc/annotation_service.html for available
             fields, extra *kwargs* and more.

        Example:

        >>> mv.getvariant(1017, email='abc@example.com')
        >>> mv.getvariant('1017', fields='symbol,name,entrezgene,refseq')
        >>> mv.getvariant('1017', fields='symbol,name,entrezgene,refseq.rna')
        >>> mv.getvariant('1017', fields=['symbol', 'name', 'pathway.kegg'])
        >>> mv.getvariant('ENSG00000123374', fields='all')

        .. Hint:: The supported field names passed to **fields** parameter can be found from
                  any full gene object (when **fields="all"**). Note that field name supports dot
                  notation for nested data structure as well, e.g. you can pass "refseq.rna" or
                  "pathway.kegg".
        '''
        #if fields:
        #    kwargs['fields'] = self._format_list(fields)
        if 'filter' in kwargs:
            kwargs['fields'] = self._format_list(kwargs['filter'])
        _url = self.url + '/variant/' + str(geneid)
        return self._get(_url, kwargs)

    def _getvariants_inner(self, geneids, **kwargs):
        _kwargs = {'ids': self._format_list(geneids)}
        _kwargs.update(kwargs)
        _url = self.url + '/variant/'
        return self._post(_url, _kwargs)

    def getvariants(self, ids, fields=None, **kwargs):
        '''Return the list of gene objects for the given list of geneids.
        This is a wrapper for POST query of "/gene" service.

        :param ids: a list or comm-sep HGVS ids
        :param fields: fields to return, a list or a comma-separated string.
                        If **fields="all"**, all available fields are returned
        :param email: optionally, pass your email to help us to track usage
        :param dataframe: return object as DataFrame (requires Pandas).
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.

        :return: a list of gene objects or a pandas DataFrame object (when **as_dataframe** is True)

        :ref: http://myvariant.info/doc/annotation_service.html for available
                fields, extra *kwargs* and more.

        Example:
        >>> vars = ['chr1:g.866422C>T',
                 'chr1:g.876664G>A',
                 'chr1:g.69635G>C',
                 'chr1:g.69869T>A',
                 'chr1:g.881918G>A',
                 'chr1:g.865625G>A',
                 'chr1:g.69892T>C',
                 'chr1:g.879381C>T',
                 'chr1:g.878330C>G']

        >>> mv.getvariants(vars, email='abc@example.com')
        >>> mv.getvariants(vars, fields="dbnsfp.cadd.phred")
        >>> mv.getvariants('chr1:g.876664G>A,chr1:g.881918G>A', fields="all")
        >>> mv.getvariants(['chr1:g.876664G>A', 'chr1:g.881918G>A'], as_dataframe=True)
        
        .. Hint:: A large list of more than 1000 input ids will be sent to the backend
                  web service in batches (1000 at a time), and then the results will be
                  concatenated together. So, from the user-end, it's exactly the same as
                  passing a shorter list. You don't need to worry about saturating our
                  backend servers.
        '''
        if isinstance(ids, str_types):
            ids = ids.split(',')
        if (not (isinstance(ids, (list, tuple)) and len(ids) > 0)):
            raise ValueError('input "variantids" must be non-empty list or tuple.')   
        if fields:
            kwargs['fields'] = self._format_list(fields)
        verbose = kwargs.pop('verbose', True)
        dataframe = kwargs.pop('dataframe', None)
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        query_fn = lambda ids: self._getvariants_inner(ids, **kwargs)
        out = []
        for hits in self._repeated_query(query_fn, ids, verbose=verbose):
            if return_raw:
                out.append(hits)   # hits is the raw response text
            else:
                out.extend(hits)
        if return_raw and len(out) == 1:
            out = out[0]
        if dataframe:
            out = self._dataframe(out, dataframe)
        return out

    def queryvariant(self, q, **kwargs):
        '''Return  the query result.
        This is a wrapper for GET query of "/query?q=<query>" service.

        :param q: a query string, detailed query syntax `here <http://mygene.info/doc/query_service.html#query-syntax>`_
        :param fields: fields to return, a list or a comma-separated string.
                        If **fields="all"**, all available fields are returned
        :param species: optionally, you can pass comma-separated species names
                        or taxonomy ids. Default: human,mouse,rat.
        :param size:   the maximum number of results to return (with a cap
                       of 1000 at the moment). Default: 10.
        :param skip:   the number of results to skip. Default: 0.
        :param sort:   Prefix with "-" for descending order, otherwise in ascending order.
                       Default: sort by matching scores in decending order.
        :param entrezonly: if True, return only matching entrez genes, otherwise, including matching
                           Ensemble-only genes (those have no matching entrez genes).
        :param email: optionally, pass your email to help us to track usage
        :param as_dataframe: if True, return object as DataFrame (requires Pandas).
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.

        :return: a dictionary with returned gene hits or a pandas DataFrame object (when **as_dataframe** is True)

        :ref: http://mygene.info/doc/query_service.html for available
              fields, extra *kwargs* and more.

        Example:

        >>> mv.queryvariant('cdk2')
        >>> mv.queryvariant('reporter:1000_at')
        >>> mv.queryvariant('symbol:cdk2', species='human')
        >>> mv.queryvariant('symbol:cdk*', species=10090, size=5, as_dataframe=True)
        >>> mv.queryvariant('q=chrX:151073054-151383976', species=9606)

        '''
        as_dataframe = kwargs.pop('as_dataframe', False)
        kwargs.update({'q': q})
        _url = self.url + '/query'
        out = self._get(_url, kwargs)
        if as_dataframe:
            out = self._as_dataframe(out, False)
        return out

    def _queryvariants_inner(self, qterms, **kwargs):
        _kwargs = {'q': self._format_list(qterms)}
        _kwargs.update(kwargs)
        _url = self.url + '/query'
        return self._post(_url, _kwargs)

    def queryvariants(self, q, scopes=None, **kwargs):
        '''Return the batch query result.
        This is a wrapper for POST query of "/query" service.

        :param qterms: a list of query terms, or a string of comma-separated query terms.
        :param scopes:  type of types of identifiers, either a list or a comma-separated fields to specify type of
                       input qterms, e.g. "entrezgene", "entrezgene,symbol", ["ensemblgene", "symbol"]
                       refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
                       of fields.
        :param fields: fields to return, a list or a comma-separated string.
                        If **fields="all"**, all available fields are returned
        :param species: optionally, you can pass comma-separated species names
                          or taxonomy ids. Default: human,mouse,rat.
        :param entrezonly:  if True, return only matching entrez genes, otherwise, including matching
                             Ensemble-only genes (those have no matching entrez genes).

        :param returnall:   if True, return a dict of all related data, including dup. and missing qterms
        :param verbose:     if True (default), print out infomation about dup and missing qterms
        :param email: optionally, pass your email to help us to track usage
        :param as_dataframe: if True, return object as DataFrame (requires Pandas).
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.

        :return: a list of gene objects or a pandas DataFrame object (when **as_dataframe** is True)

        :ref: http://mygene.info/doc/query_service.html for available
              fields, extra *kwargs* and more.

        Example:

        >>> mv.queryvariants(['DDX26B', 'CCDC83'], scopes='symbol', species=9606)
        >>> mv.queryvariants(['1255_g_at', '1294_at', '1316_at', '1320_at'], scopes='reporter')
        >>> mv.queryvariants(['NM_003466', 'CDK2', 695, '1320_at', 'Q08345'],
        ...              scopes='refseq,symbol,entrezgene,reporter,uniprot', species='human')
        >>> mv.queryvariants(['1255_g_at', '1294_at', '1316_at', '1320_at'], scopes='reporter',
        ...              fields='ensembl.gene,symbol', as_dataframe=True)

        .. Hint:: :py:meth:`queryvariants` is perfect for doing id mappings.

        .. Hint:: Just like :py:meth:`getvariants`, passing a large list of ids (>1000) to :py:meth:`queryvariants` is perfectly fine.

        '''
        if isinstance(q, str_types):
            qterms = q.split(',')
        if (not (isinstance(q, (list, tuple)) and len(q) > 0)):
            raise ValueError('input "qterms" must be non-empty list or tuple.')

        if scopes:
            kwargs['scopes'] = self._format_list(scopes)
        if 'fields' in kwargs:
            kwargs['fields'] = self._format_list(kwargs['fields'])
        returnall = kwargs.pop('returnall', False)
        verbose = kwargs.pop('verbose', True)
        dataframe = kwargs.pop('dataframe', None)
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        out = []
        li_missing = []
        li_dup = []
        li_query = []
        query_fn = lambda q: self._queryvariants_inner(q, **kwargs)
        for hits in self._repeated_query(query_fn, q, verbose=verbose):
            if return_raw:
                out.append(hits)   # hits is the raw response text
            else:
                out.extend(hits)
                for hit in hits:
                    if hit.get('notfound', False):
                        li_missing.append(hit['query'])
                    else:
                        li_query.append(hit['query'])

        if verbose:
            print("Finished.")
        if return_raw:
            if len(out) == 1:
                out = out[0]
            return out
        if dataframe:
            out = self._dataframe(out, dataframe)

        # check dup hits
       # if li_query:
       #     li_dup = [(query, cnt) for query, cnt in list_itemcnt(li_query) if cnt > 1]
       # del li_query

        if verbose:
            if li_dup:
                print("{0} input query terms found dup hits:".format(len(li_dup)))
                print("\t"+str(li_dup)[:100])
            if li_missing:
                print("{0} input query terms found no hit:".format(len(li_missing)))
                print("\t"+str(li_missing)[:100])
        if returnall:
            return {'out': out, 'dup': li_dup, 'missing': li_missing}
        else:
            if verbose and (li_dup or li_missing):
                print('Pass "returnall=True" to return complete lists of duplicate or missing query terms.')
            return out

mv=MyVariantInfo()
