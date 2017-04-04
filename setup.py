import re
import sys

from setuptools import setup

__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('kindel/__init__.py').read()).group(1)

if sys.version_info[0] < 3:
      sys.exit('Kindel requires Python 3.0 or greater')

setup(name='kindel',
      version=__version__,
      description='Indel-aware offline consensus calling for DNA alignments in SAM/BAM format',
      url='https://github.com/bede/kindel',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['kindel'],
      zip_safe=True,
      install_requires=['argh==0.26.2',
                        'tqdm==4.11.2',
                        'numpy==1.12.0',
                        'scipy==0.19.0',
                        'pandas==0.19.2',
                        'plotly==2.0.5',
                        'biopython==1.68',
                        'simplesam==0.0.4'],
      entry_points = {'console_scripts':['kindel=kindel.cli:main']})
