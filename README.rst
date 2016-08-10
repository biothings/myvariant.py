.. image:: https://badge.fury.io/py/myvariant.svg
    :target: https://pypi.python.org/pypi/myvariant

.. image:: https://img.shields.io/pypi/pyversions/myvariant.svg
    :target: https://pypi.python.org/pypi/myvariant

.. image:: https://img.shields.io/pypi/format/myvariant.svg
    :target: https://pypi.python.org/pypi/myvariant

Intro
=====

MyVariant.Info_ provides simple-to-use REST web services to query/retrieve genetic variant annotation data. It's designed with simplicity and performance emphasized. *myvariant*, is an easy-to-use Python wrapper to access MyVariant.Info_ services.

.. _MyVariant.Info: http://myvariant.info
.. _requests: https://pypi.python.org/pypi/requests

Requirements
============
    python >=2.6 (including python3)

    requests_ (install using "pip install requests")

Optional dependencies
======================
    `pandas <http://pandas.pydata.org>`_ (install using "pip install pandas") is required for returning a list of variant objects as `DataFrame <http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe>`_.

Installation
=============

    Option 1
          ::

           pip install myvariant

    Option 2
          download/extract the source code and run::

           python setup.py install

    Option 3
          install the latest code directly from the repository::

            pip install -e git+https://github.com/sulab/myvariant.py

Version history
===============

    `CHANGES.txt <https://raw.githubusercontent.com/sulab/myvariant.py/master/CHANGES.txt>`_

Tutorial
=========

    `Access ClinVar Data from MyVariant.info Services <https://cdn.rawgit.com/SuLab/myvariant.info/master/docs/ipynb/myvariant_clinvar_demo.html>`_ (the raw ipynb file is `here <https://raw.githubusercontent.com/SuLab/myvariant.info/master/docs/ipynb/myvariant_clinvar_demo.ipynb>`_)

TODO

Documentation
=============

    http://myvariant-py.readthedocs.org/

Usage
=====

.. code-block:: python

    In [1]: import myvariant

    In [2]: mv = myvariant.MyVariantInfo()

    In [3]: mv.getvariant("chr7:g.140453134T>C")
    Out[3]:  #output below is collapsed
    {"_id": "chr7:g.140453134T>C",
     "_version": 1,
     "cadd": {...},
     "cosmic": {...},
     "dbnsfp": {...},
     "dbsnp": {...},
     "docm": {...},
     "mutdb": {...},
     "snpeff": {...},
     "vcf": {
        "alt": "C",
        "position": "140453134",
        "ref": "T"
     }}

    In [4]: mv.getvariant("chr7:g.140453134T>C", fields='cosmic,snpeff')
    Out[4]:
    {'_id': 'chr7:g.140453134T>C',
     '_version': 1,
     'snpeff': {'ann': {'transcript_biotype': 'Coding',
       'gene_id': 'BRAF',
       'effect': 'missense_variant',
       'putative_impact': 'MODERATE',
       'cds': {'length': '2301', 'position': '1801'},
       'feature_type': 'transcript',
       'gene_name': 'BRAF',
       'feature_id': 'NM_004333.4',
       'hgvs_p': 'p.Lys601Glu',
       'hgvs_c': 'c.1801A>G',
       'rank': '15',
       'total': '18',
       'protein': {'length': '766', 'position': '601'},
       'cdna': {'length': '2946', 'position': '1862'}}},
     'cosmic': {'mut_freq': 0.07,
      'alt': 'G',
      'mut_nt': 'A>G',
      'tumor_site': 'upper_aerodigestive_tract',
      'ref': 'A',
      'chrom': '7',
      'hg19': {'start': 140453134, 'end': 140453134},
      'cosmic_id': 'COSM478'}
     }

    In [5]: mv.getvariant("chr7:g.140453134T>C", fields=['cosmic.tumor_site', 'snpeff.ann.gene_name'])
    Out[5]:
    {'_id': 'chr7:g.140453134T>C',
     '_version': 1,
     'snpeff': {'ann': {'gene_name': 'BRAF'}},
     'cosmic': {'tumor_site': 'upper_aerodigestive_tract'}
    }

    In [6]: mv.getvariants(['chr1:g.866422C>T', 'chr1:g.876664G>A','chr1:g.69635G>C'])
    Out[6]:
    [{'_id': 'chr1:g.866422C>T',
       ...
     },
     {'_id': 'chr1:g.876664G>A',
      ...
     },
     {'_id': 'chr1:g.69635G>C',
      ...
     }]

    In [7]: mv.getvariants(['chr1:g.866422C>T', 'chr1:g.876664G>A','chr1:g.69635G>C'],
    fields='cadd.phred,dbsnp.rsid')
    Out[7]:
    [{'query': 'chr1:g.866422C>T',
      '_id': 'chr1:g.866422C>T',
      'dbsnp': {'rsid': 'rs139210662'},
      'cadd': {'phred': 14.31}},
     {'query': 'chr1:g.876664G>A',
      '_id': 'chr1:g.876664G>A',
      'dbsnp': {'rsid': 'rs571654307'},
      'cadd': {'phred': 9.971}},
     {'query': 'chr1:g.69635G>C',
      '_id': 'chr1:g.69635G>C',
      'dbsnp': {'rsid': 'rs541766448'},
      'cadd': {'phred': 6.123}}]

    In [8]: mv.getvariants(['chr1:g.866422C>T', 'chr1:g.876664G>A','chr1:g.69635G>C'],
    fields='cadd.phred,dbsnp.rsid', as_dataframe=True)
    Out[8]:
                                   _id  cadd.phred   dbsnp.rsid
    query
    chr1:g.866422C>T  chr1:g.866422C>T      14.310  rs139210662
    chr1:g.876664G>A  chr1:g.876664G>A       9.971  rs571654307
    chr1:g.69635G>C    chr1:g.69635G>C       6.123  rs541766448

    In [9]: mv.query('dbsnp.rsid:rs58991260', fields='dbsnp')
    Out[9]:
    {'total': 1,
     'hits': [{'_score': 17.48471,
       '_id': 'chr1:g.218631822G>A',
       'dbsnp': {'class': 'SNV',
        'gmaf': 0.02157,
        'vartype': 'snp',
        'flags': ['ASP', 'G5', 'G5A', 'GNO', 'KGPhase1', 'KGPhase3', 'SLO'],
        'var_subtype': 'ts',
        'alleles': [{'freq': 0.9784, 'allele': 'G'},
         {'freq': 0.02157, 'allele': 'A'}],
        'allele_origin': 'unspecified',
        'chrom': '1',
        'hg19': {'start': 218631822, 'end': 218631823},
        'validated': True,
        'dbsnp_build': 129,
        'alt': 'A',
        'rsid': 'rs58991260',
        'ref': 'G'}}],
     'took': 24,
     'max_score': 17.48471}


    In [10]: mv.query('snpeff.ann.gene_name:cdk2 AND dbnsfp.polyphen2.hdiv.pred:D',
    fields='dbnsfp.polyphen2.hdiv')
    Out[10]:
    {'total': 1188,
     'hits': [{'dbnsfp': {'polyphen2': {'hdiv': {'rankscore': 0.89865,
          'pred': 'D',
          'score': 1.0}}},
       '_score': 8.343648,
       '_id': 'chr12:g.56359720C>T'},
      {'dbnsfp': {'polyphen2': {'hdiv': {'rankscore': 0.89865,
          'pred': 'D',
          'score': [1.0, 0.957, 0.998]}}},
       '_score': 8.343648,
       '_id': 'chr12:g.56360819G>C'},

       ...

      {'dbnsfp': {'polyphen2': {'hdiv': {'rankscore': 0.89865,
          'pred': 'D',
          'score': 1.0}}},
       '_score': 8.343648,
       '_id': 'chr12:g.56360853G>A'}],
       'took': 3521,
       'max_score': 8.343648}


    In [11]: mv.query('chr1:69000-70000', fields='cadd.phred')
    Out[11]:
    {'total': 3,
     'hits': [
      {'_score': 14.155852, '_id': 'chr1:g.69428T>G', 'cadd': {'phred': 12.14}},
      {'_score': 14.148425, '_id': 'chr1:g.69511A>G', 'cadd': {'phred': 8.98}},
      {'_score': 3.5420983, '_id': 'chr1:g.69538G>A', 'cadd': {'phred': 7.339}}],
     'took': 725,
     'max_score': 14.155852}

    In [12]: mv.querymany(['rs58991260', 'rs2500'], scopes='dbsnp.rsid', fields='dbsnp')
    Finished.
    Out[12]:
    [{'query': 'rs58991260',
      '_id': 'chr1:g.218631822G>A',
      'dbsnp': {'class': 'SNV',
       'gmaf': 0.02157,
       'vartype': 'snp',
       'flags': ['ASP', 'G5', 'G5A', 'GNO', 'KGPhase1', 'KGPhase3', 'SLO'],
       'var_subtype': 'ts',
       'alleles': [{'freq': 0.9784, 'allele': 'G'},
        {'freq': 0.02157, 'allele': 'A'}],
       'allele_origin': 'unspecified',
       'chrom': '1',
       'hg19': {'start': 218631822, 'end': 218631823},
       'validated': True,
       'dbsnp_build': 129,
       'alt': 'A',
       'rsid': 'rs58991260',
       'ref': 'G'}},
     {'query': 'rs2500',
      '_id': 'chr11:g.66397320A>G',
      'dbsnp': {'class': 'SNV',
       'vartype': 'snp',
       'flags': ['ASP', 'INT', 'RV', 'U3'],
       'var_subtype': 'ts',
       'alleles': [{'allele': 'A'}, {'allele': 'G'}],
       'allele_origin': 'unspecified',
       'chrom': '11',
       'hg19': {'start': 66397320, 'end': 66397321},
       'dbsnp_build': 36,
       'alt': 'G',
       'ref': 'A',
       'rsid': 'rs2500',
       'validated': False}}]

    In [13]: mv.querymany(['RCV000083620', 'RCV000083584'],
    scopes='clinvar.rcv_accession', fields='clinvar')
    Finished.
    Out[13]:
    [{'query': 'RCV000083620',
      'clinvar': {'type': 'single nucleotide variant',
       'gene': {'id': 5009, 'symbol': 'OTC'},
       'origin': 'unknown',
       'last_evaluated': 'None',
       'other_ids': 'dbSNP:72558473;',
       'clinvar_id': 97371,
       'hgvs': {'genomic': ['NG_008471.1:g.64470C>T',
         'NC_000023.11:g.38411952C>T',
         'NC_000023.10:g.38271205C>T'],
        'coding': 'NM_000531.5:c.958C>T'},
       'chrom': 'X',
       'cytogenic': 'Xp11.4',
       'name': 'NM_000531.5(OTC):c.958C>T (p.Arg320Ter)',
       'number_submitters': 1,
       'alt': 'T',
       'hg19': {'start': 38271205, 'end': 38271205},
       'allele_id': 103263,
       'rcv_accession': 'RCV000083620',
       'review_status': 'classified by single submitter',
       'clinical_significance': 'Pathogenic',
       'rsid': 'rs72558473',
       'ref': 'C'},
      '_id': 'chrX:g.38271205C>T'},
     {'query': 'RCV000083584',
      'clinvar': {'type': 'Deletion',
       'gene': {'id': 5009, 'symbol': 'OTC'},
       'origin': 'unknown',
       'last_evaluated': 'None',
       'other_ids': 'dbSNP:72558452;',
       'clinvar_id': 97337,
       'hgvs': {'genomic': ['NG_008471.1:g.61493_61495delGAG',
         'NC_000023.11:g.38408975_38408977delGAG',
         'NC_000023.10:g.38268228_38268230delGAG'],
        'coding': 'NM_000531.5:c.817_819delGAG'},
       'chrom': 'X',
       'cytogenic': 'Xp11.4',
       'name': 'NM_000531.5(OTC):c.817_819delGAG (p.Glu273del)',
       'number_submitters': 1,
       'alt': '-',
       'hg19': {'start': 38268228, 'end': 38268230},
       'allele_id': 103229,
       'rcv_accession': 'RCV000083584',
       'review_status': 'classified by single submitter',
       'clinical_significance': 'Pathogenic',
       'rsid': 'rs72558452',
       'ref': 'GAG'},
      '_id': 'chrX:g.38268228_38268230del'}]

    In [14]: mv.querymany(['rs2500', 'RCV000083611', 'COSM1392449'],
    scopes='clinvar.rcv_accession,dbsnp.rsid,cosmic.cosmic_id', fields='vcf', as_dataframe=1)
    Finished.
    Out[14]:
                                  _id vcf.alt vcf.position vcf.ref
    query
    rs2500        chr11:g.66397320A>G       G     66397320       A
    RCV000083611   chrX:g.38271176A>G       G     38271176       A
    COSM1392449   chr19:g.30935013C>T       T     30935013       C


    In [15]: mv.querymany(['rs58991260', 'rs2500', 'NA_TEST'], scopes='dbsnp.rsid', fields='dbsnp')
    Finished.
    1 input query terms found no hit:
            ['NA_TEST']
    Pass "returnall=True" to return complete lists of duplicate or missing query terms.
    Out[15]:
    [{'query': 'rs58991260',
      '_id': 'chr1:g.218631822G>A',
      'dbsnp': {'class': 'SNV',
       'gmaf': 0.02157,
       'vartype': 'snp',
       'flags': ['ASP', 'G5', 'G5A', 'GNO', 'KGPhase1', 'KGPhase3', 'SLO'],
       'var_subtype': 'ts',
       'alleles': [{'freq': 0.9784, 'allele': 'G'},
        {'freq': 0.02157, 'allele': 'A'}],
       'allele_origin': 'unspecified',
       'chrom': '1',
       'hg19': {'start': 218631822, 'end': 218631823},
       'validated': True,
       'dbsnp_build': 129,
       'alt': 'A',
       'rsid': 'rs58991260',
       'ref': 'G'}},
     {'query': 'rs2500',
      '_id': 'chr11:g.66397320A>G',
      'dbsnp': {'class': 'SNV',
       'vartype': 'snp',
       'flags': ['ASP', 'INT', 'RV', 'U3'],
       'var_subtype': 'ts',
       'alleles': [{'allele': 'A'}, {'allele': 'G'}],
       'allele_origin': 'unspecified',
       'chrom': '11',
       'hg19': {'start': 66397320, 'end': 66397321},
       'dbsnp_build': 36,
       'alt': 'G',
       'ref': 'A',
       'rsid': 'rs2500',
       'validated': False}},
     {'query': 'NA_TEST', 'notfound': True}]


Contact
========
Drop us any feedback at: help@myvariant.info or on twitter `@myvariantinfo <https://twitter.com/myvariantinfo>`_.
