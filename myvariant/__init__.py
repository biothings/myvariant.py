# -*- coding: utf-8 -*-
'''
Python Client for MyVariant.Info services
'''
from __future__ import print_function
import sys
import time
import requests
from itertools import islice
from collections import Iterable
try:
    from pandas import DataFrame
    from pandas.io.json import json_normalize
    df_avail = True
except:
    df_avail = False

__version__ = '0.1.1'

if sys.version_info[0] == 3:
    str_types = str
else:
    str_types = (str, unicode)


def alwayslist(value):
    '''If input value if not a list/tuple type, return it as a single value list.

    Example:

    >>> x = 'abc'
    >>> for xx in alwayslist(x):
    ...     print xx
    >>> x = ['abc', 'def']
    >>> for xx in alwayslist(x):
    ...     print xx

    '''
    if isinstance(value, (list, tuple)):
        return value
    else:
        return [value]


def safe_str(s, encoding='utf-8'):
    '''if input is an unicode string, do proper encoding.'''
    try:
        _s = str(s)
    except UnicodeEncodeError:
        _s = s.encode(encoding)
    return _s


def list_itemcnt(list):
    '''Return number of occurrence for each type of item in the list.'''
    x = {}
    for item in list:
        if item in x:
            x[item] += 1
        else:
            x[item] = 1
    return [(i, x[i]) for i in x]


def iter_n(iterable, n, with_cnt=False):
    '''
    Iterate an iterator by chunks (of n)
    if with_cnt is True, return (chunk, cnt) each time
    '''
    it = iter(iterable)
    if with_cnt:
        cnt = 0
    while True:
        chunk = tuple(islice(it, n))
        if not chunk:
            return
        if with_cnt:
            cnt += len(chunk)
            yield (chunk, cnt)
        else:
            yield chunk


def get_hgvs(input_vcf):
    '''From the input vcf file (filename or file handle), return a generator
       of genomic based HGVS ids.
    '''
    if isinstance(input_vcf, str_types):
        # if input_vcf is a string, open it as a file
        in_f = open(input_vcf)
    else:
        # otherwise it should be a file handle already
        in_f = input_vcf
    for row in in_f:
        if row[0] == '#':
            continue
        row = row.strip()
        if row:
            row = row.split('\t')
            if row[0].lower().startswith('chr'):
                row[0] = row[0][3:]
            yield get_hgvs_from_vcf(row[0], row[1], row[3], row[4])


def get_hgvs_from_vcf(chr, pos, ref, alt):
    '''get a valid hgvs name from VCF-style "chr, pos, ref, alt" data.'''
    if len(ref) == len(alt) == 1:
        # this is a SNP
        hgvs = 'chr{0}:g.{1}{2}>{3}'.format(chr, pos, ref, alt)
    elif len(ref) > 1 and len(alt) == 1:
        # this is a deletion:
        if ref[0] == alt:
            start = int(pos) + 1
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}del'.format(chr, start, end)
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
    elif len(ref) == 1 and len(alt) > 1:
        # this is a insertion
        if alt[0] == ref:
            hgvs = 'chr{0}:g.{1}_{2}ins'.format(chr, pos, int(pos) + 1)
            ins_seq = alt[1:]
            hgvs += ins_seq
        else:
            hgvs = 'chr{0}:g.{1}delins{2}'.format(chr, pos, alt)
    elif len(ref) > 1 and len(alt) > 1:
        end = int(pos) + len(alt) - 1
        hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
    else:
        raise ValueError("Cannot convert {} into HGVS id.".format((chr, pos, ref, alt)))
    return hgvs


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

    def _dataframe(self, var_obj, dataframe, df_index=True):
        """
        converts gene object to DataFrame (pandas)
        """
        if not df_avail:
            print("Error: pandas module must be installed for as_dataframe option.")
            return
        # if dataframe not in ["by_source", "normal"]:
        if dataframe not in [1, 2]:
            raise ValueError("dataframe must be either 1 (using json_normalize) or 2 (using DataFrame.from_dict")
        if 'hits' in var_obj:
            if dataframe == 1:
                df = json_normalize(var_obj['hits'])
            else:
                df = DataFrame.from_dict(var_obj['hits'])
        else:
            if dataframe == 1:
                df = json_normalize(var_obj)
            else:
                df = DataFrame.from_dict(var_obj)
        if df_index:
            df = df.set_index('query')
        return df

    def _get(self, url, params={}):
        debug = params.pop('debug', False)
        return_raw = params.pop('return_raw', False)
        headers = {'user-agent': "Python-requests_myvariant.py/%s (gzip)" % requests.__version__}
        res = requests.get(url, params=params, headers=headers)
        if debug:
            return res
        assert res.status_code == 200
        if return_raw:
            return res.text
        else:
            return res.json()

    def _post(self, url, params):
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

    def _repeated_query_old(self, query_fn, query_li, verbose=True, **fn_kwargs):
        '''This is deprecated, query_li can only be a list'''
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

    def _repeated_query(self, query_fn, query_li, verbose=True, **fn_kwargs):
        '''run query_fn for input query_li in a batch (self.step).
           return a generator of query_result in each batch.
           input query_li can be a list/tuple/iterable
        '''
        step = min(self.step, self.max_query)
        i = 0
        for batch, cnt in iter_n(query_li, step, with_cnt=True):
            if verbose:
                print("querying {0}-{1}...".format(i+1, cnt), end="")
            i = cnt
            query_result = query_fn(batch, **fn_kwargs)
            yield query_result
            if verbose:
                print("done.")
            if self.delay:
                time.sleep(self.delay)

    @property
    def metadata(self):
        '''Return a dictionary of MyVariant.info metadata.

        Example:

        >>> metadata = mv.metadata

        '''
        _url = self.url+'/metadata'
        return self._get(_url)

    def getvariant(self, vid, fields=None, **kwargs):
        '''Return the variant object for the give HGVS-based variant id.
        This is a wrapper for GET query of "/variant/<hgvsid>" service.

        :param geneid: entrez/ensembl gene id, entrez gene id can be either
                       a string or integer
        :param fields: fields to return, a list or a comma-separated string.
                        If not provided or **fields="all"**, all available fields are returned

        Example:

        >>> mv.getvariant('chr9:g.107620835G>A')
        >>> mv.getvariant('chr9:g.107620835G>A', fields='dbnsfp.genename')
        >>> mv.getvariant('chr9:g.107620835G>A', fields=['dbnsfp.genename', 'cadd.phred'])
        >>> mv.getvariant('chr9:g.107620835G>A', fields='all')

        .. Hint:: The supported field names passed to **fields** parameter can be found from
                  any full variant object (without **fields**, or **fields="all"**). Note that field name supports dot
                  notation for nested data structure as well, e.g. you can pass "dbnsfp.genename" or
                  "cadd.phred".
        '''
        if fields:
            kwargs['fields'] = self._format_list(fields)
        _url = self.url + '/variant/' + str(vid)
        return self._get(_url, kwargs)

    def _getvariants_inner(self, geneids, **kwargs):
        _kwargs = {'ids': self._format_list(geneids)}
        _kwargs.update(kwargs)
        _url = self.url + '/variant/'
        return self._post(_url, _kwargs)

    def getvariants(self, vids, fields=None, **kwargs):
        '''Return the list of variant annotation objects for the given list of hgvs-base varaint ids.
        This is a wrapper for POST query of "/variant" service.

        :param ids: a list/tuple/iterable or a string of comm-sep HGVS ids.
        :param fields: fields to return, a list or a comma-separated string.
                        If **fields="all"**, all available fields are returned
        :param as_dataframe: if True or 1 or 2, return object as DataFrame (requires Pandas).
                                  True or 1: using json_normalize
                                  2        : using DataFrame.from_dict
                                  otherwise: return original json
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.

        :return: a list of variant objects or a pandas DataFrame object (when **as_dataframe** is True)

        :ref: http://docs.myvariant.info/en/latest/doc/variant_annotation_service.html.

        Example:

        >>> vars = ['chr1:g.866422C>T',
        ...         'chr1:g.876664G>A',
        ...         'chr1:g.69635G>C',
        ...         'chr1:g.69869T>A',
        ...         'chr1:g.881918G>A',
        ...         'chr1:g.865625G>A',
        ...         'chr1:g.69892T>C',
        ...         'chr1:g.879381C>T',
        ...         'chr1:g.878330C>G']
        >>> mv.getvariants(vars, fields="cadd.phred")
        >>> mv.getvariants('chr1:g.876664G>A,chr1:g.881918G>A', fields="all")
        >>> mv.getvariants(['chr1:g.876664G>A', 'chr1:g.881918G>A'], as_dataframe=True)

        .. Hint:: A large list of more than 1000 input ids will be sent to the backend
                  web service in batches (1000 at a time), and then the results will be
                  concatenated together. So, from the user-end, it's exactly the same as
                  passing a shorter list. You don't need to worry about saturating our
                  backend servers.

        .. Hint:: If you need to pass a very large list of input ids, you can pass a generator
                  instead of a full list, which is more memory efficient.
        '''
        if isinstance(vids, str_types):
            vids = vids.split(',') if vids else []
        if (not (isinstance(vids, (list, tuple, Iterable)))):
            raise ValueError('input "vids" must be a list, tuple or iterable.')
        if fields:
            kwargs['fields'] = self._format_list(fields)
        verbose = kwargs.pop('verbose', True)
        dataframe = kwargs.pop('as_dataframe', None)
        df_index = kwargs.pop('df_index', True)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        query_fn = lambda vids: self._getvariants_inner(vids, **kwargs)
        out = []
        for hits in self._repeated_query(query_fn, vids, verbose=verbose):
            if return_raw:
                out.append(hits)   # hits is the raw response text
            else:
                out.extend(hits)
        if return_raw and len(out) == 1:
            out = out[0]
        if dataframe:
            out = self._dataframe(out, dataframe, df_index=df_index)
        return out

    def query(self, q, **kwargs):
        '''Return  the query result.
        This is a wrapper for GET query of "/query?q=<query>" service.

        :param q: a query string, detailed query syntax `here <http://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax>`_
        :param fields: fields to return, a list or a comma-separated string.
                        If not provided or **fields="all"**, all available fields are returned
        :param size:   the maximum number of results to return (with a cap
                       of 1000 at the moment). Default: 10.
        :param skip:   the number of results to skip. Default: 0.
        :param sort:   Prefix with "-" for descending order, otherwise in ascending order.
                       Default: sort by matching scores in decending order.
        :param as_dataframe: if True or 1 or 2, return object as DataFrame (requires Pandas).
                                  True or 1: using json_normalize
                                  2        : using DataFrame.from_dict
                                  otherwise: return original json

        :return: a dictionary with returned gene hits or a pandas DataFrame object (when **as_dataframe** is True)

        :ref: http://docs.myvariant.info/en/latest/doc/variant_query_service.html.

        Example:

        >>> mv.query('_exists_:dbsnp AND _exists_:cosmic')
        >>> mv.query('dbnsfp.polyphen2.hdiv.score:>0.99 AND chrom:1')
        >>> mv.query('cadd.phred:>50')
        >>> mv.query('dbnsfp.genename:CDK2', size=5)
        >>> mv.query('chrX:151073054-151383976')

        '''
        dataframe = kwargs.pop('as_dataframe', None)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        kwargs.update({'q': q})
        _url = self.url + '/query'
        out = self._get(_url, kwargs)
        if dataframe:
            out = self._dataframe(out, dataframe, df_index=False)
        return out

    def _queryvariants_inner(self, qterms, **kwargs):
        _kwargs = {'q': self._format_list(qterms)}
        _kwargs.update(kwargs)
        _url = self.url + '/query'
        return self._post(_url, _kwargs)

    def querymany(self, qterms, scopes=None, **kwargs):
        '''Return the batch query result.
        This is a wrapper for POST query of "/query" service.

        :param qterms: a list/tuple/iterable of query terms, or a string of comma-separated query terms.
        :param scopes:  type of types of identifiers, either a list or a comma-separated fields to specify type of
                       input qterms, e.g. "dbsnp.rsid", "clinvar.rcv_accession", ["dbsnp.rsid", "cosmic.cosmic_id"]
                       refer to "http://docs.myvariant.info/en/latest/doc/data.html#available-fields" for full list
                       of fields.
        :param fields: fields to return, a list or a comma-separated string.
                        If not provided or **fields="all"**, all available fields are returned
        :param returnall:   if True, return a dict of all related data, including dup. and missing qterms
        :param verbose:     if True (default), print out information about dup and missing qterms
        :param as_dataframe: if True or 1 or 2, return object as DataFrame (requires Pandas).
                                  True or 1: using json_normalize
                                  2        : using DataFrame.from_dict
                                  otherwise: return original json
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.
        :return: a list of matching variant objects or a pandas DataFrame object.
        :ref: http://docs.myvariant.info/en/latest/doc/variant_query_service.html for available
              fields, extra *kwargs* and more.

        Example:

        >>> mv.querymany(['rs58991260', 'rs2500'], scopes='dbsnp.rsid')
        >>> mv.querymany(['RCV000083620', 'RCV000083611', 'RCV000083584'], scopes='clinvar.rcv_accession')
        >>> mv.querymany(['COSM1362966', 'COSM990046', 'COSM1392449'], scopes='cosmic.cosmic_id', fields='cosmic')
        >>> mv.querymany(['COSM1362966', 'COSM990046', 'COSM1392449'], scopes='cosmic.cosmic_id',
        ...              fields='cosmic.tumor_site', as_dataframe=True)

        .. Hint:: :py:meth:`queryvariants` is perfect for query variants based different ids, e.g. rsid, clinvar ids, etc.

        .. Hint:: Just like :py:meth:`getvariants`, passing a large list of ids (>1000) to :py:meth:`queryvariants` is perfectly fine.

        .. Hint:: If you need to pass a very large list of input qterms, you can pass a generator
                  instead of a full list, which is more memory efficient.

        '''
        if isinstance(qterms, str_types):
            qterms = qterms.split(',') if qterms else []
        if (not (isinstance(qterms, (list, tuple, Iterable)))):
            raise ValueError('input "qterms" must be a list, tuple or iterable.')

        if scopes:
            kwargs['scopes'] = self._format_list(scopes)
        if 'fields' in kwargs:
            kwargs['fields'] = self._format_list(kwargs['fields'])
        returnall = kwargs.pop('returnall', False)
        verbose = kwargs.pop('verbose', True)
        dataframe = kwargs.pop('as_dataframe', None)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        df_index = kwargs.pop('df_index', True)
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        out = []
        li_missing = []
        li_dup = []
        li_query = []
        query_fn = lambda qterms: self._queryvariants_inner(qterms, **kwargs)
        for hits in self._repeated_query(query_fn, qterms, verbose=verbose):
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
            out = self._dataframe(out, dataframe, df_index=df_index)

        # check dup hits
        if li_query:
            li_dup = [(query, cnt) for query, cnt in list_itemcnt(li_query) if cnt > 1]
        del li_query

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
