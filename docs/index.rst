.. MyVariant.py documentation master file, created by
   sphinx-quickstart on Thu Jul 30 17:55:51 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. _MyVariant.Info: http://myvariant.info
.. _biothings_client: https://pypi.org/project/biothings-client/
.. _myvariant: https://pypi.org/project/myvariant/

Welcome to MyVariant.py's documentation!
========================================

MyVariant.Info_ provides simple-to-use REST web services to query/retrieve variant annotation data.
It's designed with simplicity and performance emphasized. *myvariant*, is an easy-to-use Python wrapper
to access MyVariant.Info_ services.

.. Note::
    As of v1.0.0, myvariant_ Python package is now a thin wrapper of underlying biothings_client_ package,
    a universal Python client for all `BioThings APIs <http://biothings.io>`_, including MyVariant.info_.
    The installation of myvariant_ will install biothings_client_ automatically. The following code snippets
    are essentially equivalent:

    * Continue using myvariant_ package

        .. code-block:: python

            In [1]: import myvariant
            In [2]: mv = myvariant.MyVariantInfo()

    * Use biothings_client_ package directly

        .. code-block:: python

            In [1]: from biothings_client import get_client
            In [2]: mv = get_client('variant')

    After that, the use of ``mv`` instance is exactly the same.


.. toctree::
   :maxdepth: 2
   index

Requirements
============

    python >=2.7 (including python3)

    (Python 2.6 might still work, not it's not supported any more since v1.0.0)

    biothings_client_ (>=0.2.0, install using "pip install biothings_client")

Optional dependencies
======================

    `pandas <http://pandas.pydata.org>`_ (install using "pip install pandas") is required for returning a list
    of gene objects as `DataFrame <http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe>`_.

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

            pip install -e git+https://github.com/biothings/myvariant.py

Version history
===============

    `CHANGES.txt <https://raw.githubusercontent.com/biothings/myvariant.py/master/CHANGES.txt>`_

Tutorial
=========

    `Access ClinVar Data from MyVariant.info Services <https://cdn.rawgit.com/biothings/myvariant.info/master/docs/ipynb/myvariant_clinvar_demo.html>`_ (the raw ipynb file is `here <https://raw.githubusercontent.com/biothings/myvariant.info/master/docs/ipynb/myvariant_clinvar_demo.ipynb>`_)

API
======

.. py:module:: myvariant
.. autofunction:: alwayslist
.. autofunction:: get_hgvs_from_vcf
.. autofunction:: format_hgvs
.. autoclass:: MyVariantInfo
    :members:
    :inherited-members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
