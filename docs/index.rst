.. MyVariant.py documentation master file, created by
   sphinx-quickstart on Thu Jul 30 17:55:51 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. _MyVariant.Info: http://myvariant.info
.. _requests: https://pypi.python.org/pypi/requests

Welcome to MyVariant.py's documentation!
========================================

MyVariant.Info_ provides simple-to-use REST web services to query/retrieve variant annotation data. It's designed with simplicity and performance emphasized. *myvariant*, is an easy-to-use Python wrapper to access MyVariant.Info_ services.

.. toctree::
   :maxdepth: 2
   index

Requirements
============
    python >=2.6 (including python3)

    requests_ (install using "pip install requests")

Optional dependencies
======================
    `pandas <http://pandas.pydata.org>`_ (install using "pip install pandas") is required for returning a list of gene objects as `DataFrame <http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe>`_.

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

.. * `ID mapping using mygene module in Python <http://nbviewer.ipython.org/6771106>`_

TODO


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
