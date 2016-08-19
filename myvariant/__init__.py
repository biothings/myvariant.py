'''
Python Client for MyVariant.Info services
'''
from __future__ import print_function
import sys
import os
import time
from itertools import islice
from collections import Iterable

import requests
try:
    from pandas import DataFrame
    from pandas.io.json import json_normalize
    df_avail = True
except:
    df_avail = False

try:
    import requests_cache
    caching_avail = True
except:
    caching_avail = False

__version__ = '0.3.1'

if sys.version_info[0] == 3:
    str_types = str
else:
    str_types = (str, unicode)


class ScanError(Exception):
    # for errors in scan search type
    pass


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


def get_hgvs_from_vcf(input_vcf):
    '''From the input VCF file (filename or file handle), return a generator
       of genomic based HGVS ids.

       :param input_vcf: input VCF file, can be a filename or a file handle

       :returns: a generator of genomic based HGVS ids. To get back a list
                 instead, using *list(get_hgvs_from_vcf("your_vcf_file"))*

       .. NOTE:: This is a lightweight VCF parser to return valid genomic-based
                 HGVS ids from the *input_vcf* file. For more sophisticated VCF
                 parser, consider using `PyVCF <https://pypi.python.org/pypi/PyVCF>`_
                 module.

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
            for alt in row[4].split(','):
                yield format_hgvs(row[0], row[1], row[3], alt)


def _normalized_vcf(chrom, pos, ref, alt):
    """If both ref/alt are > 1 base, and there are overlapping from the left,
       we need to trim off the overlapping bases.
       In the case that ref/alt is like this:
           CTTTT/CT    # with >1 overlapping bases from the left
       ref/alt should be normalized as TTTT/T, more examples:
            TC/TG --> C/G
       and pos should be fixed as well.
    """
    for i in range(max(len(ref), len(alt))):
        _ref = ref[i] if i < len(ref) else None
        _alt = alt[i] if i < len(alt) else None
        if _ref is None or _alt is None or _ref != _alt:
            break

    # _ref/_alt cannot be both None, if so,
    # ref and alt are exactly the same,
    # something is wrong with this VCF record
    assert not (_ref is None and _alt is None)

    _pos = int(pos)
    if _ref is None or _alt is None:
        # if either is None, del or ins types
        _pos = _pos + i - 1
        _ref = ref[i-1:]
        _alt = alt[i-1:]
    else:
        # both _ref/_alt are not None
        _pos = _pos + i
        _ref = ref[i:]
        _alt = alt[i:]

    return (chrom, _pos, _ref, _alt)


def format_hgvs(chrom, pos, ref, alt):
    '''get a valid hgvs name from VCF-style "chrom, pos, ref, alt" data.

    Example:

        >>> myvariant.format_hgvs("1", 35366, "C", "T")
        >>> myvariant.format_hgvs("2", 17142, "G", "GA")
        >>> myvariant.format_hgvs("MT", 8270, "CACCCCCTCT", "C")
        >>> myvariant.format_hgvs("X", 107930849, "GGA", "C")

    '''
    chrom = str(chrom)
    if chrom.lower().startswith('chr'):
        # trim off leading "chr" if any
        chrom = chrom[3:]
    if len(ref) == len(alt) == 1:
        # this is a SNP
        hgvs = 'chr{0}:g.{1}{2}>{3}'.format(chrom, pos, ref, alt)
    elif len(ref) > 1 and len(alt) == 1:
        # this is a deletion:
        if ref[0] == alt:
            start = int(pos) + 1
            end = int(pos) + len(ref) - 1
            if start == end:
                hgvs = 'chr{0}:g.{1}del'.format(chrom, start)
            else:
                hgvs = 'chr{0}:g.{1}_{2}del'.format(chrom, start, end)
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chrom, pos, end, alt)
    elif len(ref) == 1 and len(alt) > 1:
        # this is an insertion
        if alt[0] == ref:
            hgvs = 'chr{0}:g.{1}_{2}ins'.format(chrom, pos, int(pos) + 1)
            ins_seq = alt[1:]
            hgvs += ins_seq
        else:
            hgvs = 'chr{0}:g.{1}delins{2}'.format(chrom, pos, alt)
    elif len(ref) > 1 and len(alt) > 1:
        if ref[0] == alt[0]:
            # if ref and alt overlap from the left, trim them first
            _chrom, _pos, _ref, _alt = _normalized_vcf(chrom, pos, ref, alt)
            return format_hgvs(_chrom, _pos, _ref, _alt)
        else:
            end = int(pos) + len(alt) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chrom, pos, end, alt)
    else:
        raise ValueError("Cannot convert {} into HGVS id.".format((chrom, pos, ref, alt)))
    return hgvs


class MyVariantInfo:
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
        self.scroll_size = 1000
        # raise requests.exceptions.HTTPError for status_code > 400
        #   but not for 404 on getvariant
        #   set to False to surpress the exceptions.
        self.raise_for_status = True
        self.default_user_agent = "myvariant.py/%s python-requests/%s" % (__version__, requests.__version__)
        self._cached = False


    def _dataframe(self, var_obj, dataframe, df_index=True):
        """
        converts variant object to DataFrame (pandas)
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

    def _get(self, url, params={}, none_on_404=False, verbose=True):
        debug = params.pop('debug', False)
        return_raw = params.pop('return_raw', False)
        headers = {'user-agent': self.default_user_agent}
        res = requests.get(url, params=params, headers=headers)
        from_cache = getattr(res, 'from_cache', False)
        if debug:
            return from_cache, res
        if none_on_404 and res.status_code == 404:
            return from_cache, None
        if self.raise_for_status:
            # raise requests.exceptions.HTTPError if not 200
            res.raise_for_status()
        if return_raw:
            return from_cache, res.text
        ret = res.json()
        return from_cache, ret

    def _post(self, url, params, verbose=True):
        return_raw = params.pop('return_raw', False)
        headers = {'content-type': 'application/x-www-form-urlencoded',
                   'user-agent': self.default_user_agent}
        res = requests.post(url, data=params, headers=headers)
        from_cache = getattr(res, 'from_cache', False)
        if self.raise_for_status:
            # raise requests.exceptions.HTTPError if not 200
            res.raise_for_status()
        if return_raw:
            return from_cache, res
        ret = res.json()
        return from_cache, ret

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
            from_cache, query_result = query_fn(batch, **fn_kwargs)
            yield query_result
            if verbose:
                cache_str = " {0}".format(self._from_cache_notification) if from_cache else ""
                print("done.{0}".format(cache_str))
            if self.delay:
                time.sleep(self.delay)

    @property
    def _from_cache_notification(self):
        ''' Notification to alert user that a cached result is being returned.'''
        return "[ from cache ]"

    def metadata(self, verbose=True, **kwargs):
        '''Return a dictionary of MyVariant.info metadata.

        Example:

        >>> metadata = mv.metadata()

        '''
        _url = self.url+'/metadata'
        from_cache, ret = self._get(_url, params=kwargs, verbose=verbose)
        if verbose and from_cache:
            print(self._from_cache_notification)
        return ret

    def set_caching(self, cache_db='myvariant_cache', verbose=True, **kwargs):
        ''' Installs a local cache for all requests.

            **cache_db** is the path to the local sqlite cache database.'''
        if caching_avail:
            requests_cache.install_cache(cache_name=cache_db, allowable_methods=('GET', 'POST'), **kwargs)
            self._cached = True
            if verbose:
                print('[ Future queries will be cached in "{0}" ]'.format(os.path.abspath(cache_db + '.sqlite')))
        else:
            print("Error: The requests_cache python module is required to use request caching.")
            print("See - https://requests-cache.readthedocs.io/en/latest/user_guide.html#installation")
        return

    def stop_caching(self):
        ''' Stop caching.'''
        if self._cached and caching_avail:
            requests_cache.uninstall_cache()
            self._cached = False
        return

    def clear_cache(self):
        ''' Clear the globally installed cache. '''
        try:
            requests_cache.clear()
        except:
            pass

    def get_fields(self, search_term=None, verbose=True):
        ''' Wrapper for http://myvariant.info/v1/metadata/fields

            **search_term** is a case insensitive string to search for in available field names.
            If not provided, all available fields will be returned.


        Example:

        >>> mv.get_fields()
        >>> mv.get_fields("rsid")
        >>> mv.get_fields("sift")

        .. Hint:: This is useful to find out the field names you need to pass to **fields** parameter of other methods.

        '''
        _url = self.url + '/metadata/fields'
        if search_term:
            params = {'search': search_term}
        else:
            params = {}
        from_cache, ret = self._get(_url, params=params, verbose=verbose)
        for (k, v) in ret.items():
            # Get rid of the notes column information
            if "notes" in v:
                del v['notes']
        if verbose and from_cache:
            print(self._from_cache_notification)
        return ret

    def getvariant(self, vid, fields=None, **kwargs):
        '''Return the variant object for the give HGVS-based variant id.
        This is a wrapper for GET query of "/variant/<hgvsid>" service.

        :param vid: an HGVS-based variant id. `More about HGVS id <http://docs.myvariant.info/en/latest/doc/data.html#id-field>`_.
        :param fields: fields to return, a list or a comma-separated string.
                       If not provided or **fields="all"**, all available fields
                       are returned. See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_
                       for all available fields.

        :return: a variant object as a dictionary, or None if vid is not found.

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
        verbose = kwargs.pop('verbose', True)
        if fields:
            kwargs['fields'] = self._format_list(fields)
        _url = self.url + '/variant/' + str(vid)
        from_cache, ret = self._get(_url, kwargs, none_on_404=True, verbose=verbose)
        if verbose and from_cache:
            print(self._from_cache_notification)
        return ret

    def _getvariants_inner(self, vids, verbose=True, **kwargs):
        _kwargs = {'ids': self._format_list(vids)}
        _kwargs.update(kwargs)
        _url = self.url + '/variant/'
        return self._post(_url, _kwargs, verbose=verbose)

    def _variants_generator(self, query_fn, vids, verbose=True, **kwargs):
        ''' Convenience function to yield a batch of hits one at a yime. '''
        for hits in self._repeated_query(query_fn, vids, verbose=verbose):
            for hit in hits:
                yield hit

    def getvariants(self, vids, fields=None, **kwargs):
        '''Return the list of variant annotation objects for the given list of hgvs-base varaint ids.
        This is a wrapper for POST query of "/variant" service.

        :param ids: a list/tuple/iterable or a string of comma-separated HGVS ids.
                    `More about hgvs id <http://docs.myvariant.info/en/latest/doc/data.html#id-field>`_.
        :param fields: fields to return, a list or a comma-separated string.
                       If not provided or **fields="all"**, all available fields
                       are returned. See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_
                       for all available fields.
        :param as_generator:  if True, will yield the results in a generator.
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
        generator = kwargs.pop('as_generator', False)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        query_fn = lambda vids: self._getvariants_inner(vids, verbose=verbose, **kwargs)
        if generator:
            return self._variants_generator(query_fn, vids, verbose=verbose, **kwargs)
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

        :param q: a query string, detailed query syntax `here <http://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax>`_.
        :param fields: fields to return, a list or a comma-separated string.
                       If not provided or **fields="all"**, all available fields
                       are returned. See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_
                       for all available fields.
        :param size:   the maximum number of results to return (with a cap
                       of 1000 at the moment). Default: 10.
        :param skip:   the number of results to skip. Default: 0.
        :param sort:   Prefix with "-" for descending order, otherwise in ascending order.
                       Default: sort by matching scores in decending order.
        :param as_dataframe: if True or 1 or 2, return object as DataFrame (requires Pandas).
                                  True or 1: using json_normalize
                                  2        : using DataFrame.from_dict
                                  otherwise: return original json
        :param fetch_all: if True, return a generator to all query results (unsorted).  This can provide a very fast
                          return of all hits from a large query.
                          Server requests are done in blocks of 1000 and yielded individually.  Each 1000 block of
                          results must be yielded within 1 minute, otherwise the request will expire at server side.

        :return: a dictionary with returned variant hits or a pandas DataFrame object (when **as_dataframe** is True)
                 or a generator of all hits (when **fetch_all** is True)
        :ref: http://docs.myvariant.info/en/latest/doc/variant_query_service.html.

        Example:

        >>> mv.query('_exists_:dbsnp AND _exists_:cosmic')
        >>> mv.query('dbnsfp.polyphen2.hdiv.score:>0.99 AND chrom:1')
        >>> mv.query('cadd.phred:>50')
        >>> mv.query('dbnsfp.genename:CDK2', size=5)
        >>> mv.query('dbnsfp.genename:CDK2', fetch_all=True)
        >>> mv.query('chrX:151073054-151383976')

        .. Hint:: By default, **query** method returns the first 10 hits if the matched hits are >10. If the total number
                  of hits are less than 1000, you can increase the value for **size** parameter. For a query returns
                  more than 1000 hits, you can pass "fetch_all=True" to return a `generator <http://www.learnpython.org/en/Generators>`_
                  of all matching hits (internally, those hits are requested from the server-side in blocks of 1000).
        '''
        verbose = kwargs.pop('verbose', True)
        kwargs.update({'q': q})
        fetch_all = kwargs.get('fetch_all')
        if fetch_all in [True, 1]:
            return self._fetch_all(verbose=verbose, **kwargs)
        dataframe = kwargs.pop('as_dataframe', None)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        _url = self.url + '/query'
        from_cache, out = self._get(_url, kwargs, verbose=verbose)
        if verbose and from_cache:
            print(self._from_cache_notification)
        if dataframe:
            out = self._dataframe(out, dataframe, df_index=False)
        return out

    def _fetch_all(self, verbose=True, **kwargs):
        ''' Function that returns a generator to results.  Assumes that 'q' is in kwargs.'''
        # get the total number of hits and start the scroll_id
        _url = self.url + '/query'
        # function to get the next batch of results, automatically disables cache if we are caching
        def _batch():
            if caching_avail and self._cached:
                self._cached = False
                with requests_cache.disabled():
                    from_cache, ret = self._get(_url, params=kwargs, verbose=verbose)
                self._cached = True
            else:
                from_cache, ret = self._get(_url, params=kwargs, verbose=verbose)
            return ret
        batch = _batch()
        if verbose:
            print("Fetching {0} variant(s) . . .".format(batch['total']))
        for key in ['q', 'fetch_all']:
            kwargs.pop(key)
        while not batch.get('error', '').startswith('No results to return.'):
            if 'error' in batch:
                print(batch['error'])
                break
            if '_warning' in batch and verbose:
                print(batch['_warning'])
            for hit in batch['hits']:
                yield hit
            kwargs.update({'scroll_id': batch['_scroll_id']})
            batch = _batch()

    def _querymany_inner(self, qterms, verbose=True, **kwargs):
        _kwargs = {'q': self._format_list(qterms)}
        _kwargs.update(kwargs)
        _url = self.url + '/query'
        return self._post(_url, params=_kwargs, verbose=verbose)

    def querymany(self, qterms, scopes=None, **kwargs):
        '''Return the batch query result.
        This is a wrapper for POST query of "/query" service.

        :param qterms: a list/tuple/iterable of query terms, or a string of comma-separated query terms.
        :param scopes: specify the type (or types) of identifiers passed to **qterms**, either a list or a comma-separated fields to specify type of
                       input qterms, e.g. "dbsnp.rsid", "clinvar.rcv_accession", ["dbsnp.rsid", "cosmic.cosmic_id"].
                       See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_ for full list
                       of supported fields.
        :param fields: fields to return, a list or a comma-separated string.
                       If not provided or **fields="all"**, all available fields
                       are returned. See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_
                       for all available fields.
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

        .. Hint:: :py:meth:`querymany` is perfect for query variants based different ids, e.g. rsid, clinvar ids, etc.

        .. Hint:: Just like :py:meth:`getvariants`, passing a large list of ids (>1000) to :py:meth:`querymany` is perfectly fine.

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
        query_fn = lambda qterms: self._querymany_inner(qterms, verbose=verbose, **kwargs)
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
