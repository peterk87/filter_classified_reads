=======================
filter_classified_reads
=======================


.. image:: https://img.shields.io/pypi/v/filter_classified_reads.svg
        :target: https://pypi.python.org/pypi/filter_classified_reads

.. image:: https://travis-ci.com/peterk87/filter_classified_reads.svg?branch=master
    :target: https://travis-ci.com/peterk87/filter_classified_reads

.. image:: https://readthedocs.org/projects/filter-classified-reads/badge/?version=latest
        :target: https://filter-classified-reads.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Filter for reads from taxa of interest using Kraken2/Centrifuge classification results.


* Free software: Apache Software License 2.0
* Documentation: https://filter-classified-reads.readthedocs.io.


Features
--------

* Filter for union of reads classified to taxa of interest Kraken2_ and Centrifuge_ (by default filter for Viral reads (taxid=10239))
* Output unclassified reads along with reads from taxa of interest *or* exlude them with `--exclude-unclassified`
* screed_ for quickly filtering reads

Usage
-----

Paired-end reads with classification results by both Kraken2_ and Centrifuge_

.. code-block::

    filter_classified_reads -i /path/to/reads/R1.fq \
                            -I /path/to/reads/R2.fq \
                            -o  /path/to/reads/R1.filtered.fq \
                            -O  /path/to/reads/R2.filtered.fq \
                            -k  /path/to/kraken2/results.tsv \
                            -K  /path/to/kraken2/kreport.tsv \
                            -c  /path/to/centrifuge/results.tsv \
                            -C  /path/to/centrifuge/kreport.tsv \


Using test data in `tests/data/`:

.. code-block::

    $ filter_classified_reads -i tests/data/SRR8207674_1.viral_unclassified.seqtk_seed42_n10000.fastq.gz \
                              -I tests/data/SRR8207674_2.viral_unclassified.seqtk_seed42_n10000.fastq.gz \
                              -o r1.fq \
                              -O r2.fq \
                              -k tests/data/SRR8207674-kraken2_results.tsv \
                              -K tests/data/SRR8207674-kraken2_report.tsv \
                              -c tests/data/SRR8207674-centrifuge_results.tsv \
                              -C tests/data/SRR8207674-centrifuge_kreport.tsv

You should see the following log information:

.. code-block::

    2019-04-16 13:40:34,114 INFO: Parsing centrifuge results into DataFrame [in target_classified_reads.py:49]
    2019-04-16 13:40:34,168 INFO: Parsed n=12281 centrifuge result records into DataFrame from "tests/data/SRR8207674-centrifuge_results.tsv" [in target_classified_reads.py:57]
    2019-04-16 13:40:34,172 INFO: Parsed n=298 centrifuge Kraken-style report records into DataFrame from "tests/data/SRR8207674-centrifuge_kreport.tsv" [in target_classified_reads.py:60]
    2019-04-16 13:40:34,177 INFO: Found 7129 unclassified reads from Centrifuge results [in target_classified_reads.py:65]
    2019-04-16 13:40:34,242 INFO: Found 231 unique viral Taxonomy IDs [in target_classified_reads.py:98]
    2019-04-16 13:40:34,245 INFO: Found 2181 target reads from centrifuge results [in target_classified_reads.py:101]
    2019-04-16 13:40:34,245 INFO: Parsing kraken2 results into DataFrame [in target_classified_reads.py:49]
    2019-04-16 13:40:34,289 INFO: Parsed n=20000 kraken2 result records into DataFrame from "tests/data/SRR8207674-kraken2_results.tsv" [in target_classified_reads.py:57]
    2019-04-16 13:40:34,293 INFO: Parsed n=139 kraken2 Kraken-style report records into DataFrame from "tests/data/SRR8207674-kraken2_report.tsv" [in target_classified_reads.py:60]
    2019-04-16 13:40:34,295 INFO: Found 1737 unclassified reads from Centrifuge results [in target_classified_reads.py:65]
    2019-04-16 13:40:34,325 INFO: Found 26 unique viral Taxonomy IDs [in target_classified_reads.py:98]
    2019-04-16 13:40:34,331 INFO: Found 8345 target reads from kraken2 results [in target_classified_reads.py:101]
    2019-04-16 13:40:34,332 INFO: Found N=1701 common unclassified reads by all classification methods. [in cli.py:110]
    2019-04-16 13:40:34,333 INFO: Total viral reads=8357 [in util.py:37]
    2019-04-16 13:40:34,333 INFO: Centrifuge found n=12 target reads not found with Kraken2 [in util.py:38]
    2019-04-16 13:40:34,333 INFO: Kraken2 found n=6176 target reads not found with Centrifuge [in util.py:40]
    2019-04-16 13:40:34,338 INFO: N=1701 reads unclassified by both Centrifuge and Kraken2. [in util.py:62]
    2019-04-16 13:40:34,345 INFO: Writing n=9999 filtered reads from "tests/data/SRR8207674_1.viral_unclassified.seqtk_seed42_n10000.fastq.gz" to "r1.fq" [in cli.py:129]
    2019-04-16 13:40:34,957 INFO: Writing n=9999 filtered reads from "tests/data/SRR8207674_2.viral_unclassified.seqtk_seed42_n10000.fastq.gz" to "r2.fq" [in cli.py:134]
    2019-04-16 13:40:35,459 INFO: Done! [in cli.py:137]



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _Kraken2: https://ccb.jhu.edu/software/kraken2/
.. _Centrifuge: https://ccb.jhu.edu/software/centrifuge/manual.shtml
.. _screed: https://screed.readthedocs.io/en/latest/screed.html
