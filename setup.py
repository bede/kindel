import re
import sys

from setuptools import setup


__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('kindel/__init__.py').read()).group(1)


if sys.version_info[0] < 3:
      sys.exit('Kindel requires Python >= 3.5')


CLASSIFIERS = ['Environment :: Console',
               'Environment :: MacOS X',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Natural Language :: English',
               'Operating System :: POSIX :: Linux',
               'Operating System :: MacOS :: MacOS X',
               'Programming Language :: Python :: 3.5',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics']


setup(name='kindel',
      version=__version__,
      description='Indel-aware consensus calling for DNA alignments in BAM format',
      url='https://github.com/bede/kindel',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['kindel'],
      zip_safe=True,
      install_requires=['argh>=0.26.2',
                        'tqdm>=4.11.2',
                        'numpy>=1.12.0',
                        'scipy>=0.19.0',
                        'pandas>=0.19.2',
                        'biopython>=1.68',
                        'simplesam>=0.0.4'],
      extras_require={'plots': ['plotly>=2.0.0']},
      entry_points = {'console_scripts':['kindel=kindel.cli:main']})
