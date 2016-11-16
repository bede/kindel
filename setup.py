import re
import sys

from setuptools import setup

__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('kindel/__init__.py').read()).group(1)

if sys.version_info[0] < 3:
      sys.exit('Kindel requires Python 3.0 or greater. Please upgrade')

setup(name='kindel',
      version=__version__,
      description='Indel-aware offline consensus calling for DNA alignments in SAM/BAM format',
      url='https://github.com/bede/kindel',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['kindel'],
      zip_safe=True,
      install_requires=['argh>=0.26.2','biopython>=1.67','simplesam>=0.0.4', 'tqdm>=4.7.4'],
      entry_points = {'console_scripts':['kindel=kindel.cli:main']})
