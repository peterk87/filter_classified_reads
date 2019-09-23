#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['click==7.0',
                'pandas==0.25.1',
                'numpy==1.17.2',
                'attrs>=19.1.0']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Peter Kruczkiewicz",
    author_email='peter.kruczkiewicz@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    description="Filter for reads from taxa of interest using Kraken2/Centrifuge classification results",
    entry_points={
        'console_scripts': [
            'filter_classified_reads=filter_classified_reads.cli:main',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='filter_classified_reads',
    name='filter_classified_reads',
    packages=find_packages(include=['filter_classified_reads']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/peterk87/filter_classified_reads',
    version='0.2.0',
    zip_safe=False,
)
