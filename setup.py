from setuptools import setup

setup(name='kindel',
      version='0.1.3',
      description='Indel-aware offline consensus calling for DNA alignments in SAM/BAM format',
      url='https://github.com/bede/kindel',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['kindel'],
      zip_safe=True,
      install_requires=['argh>=0.26.2','biopython>=1.67','simplesam>=0.0.4', 'tqdm>=4.7.4'],
      entry_points = {'console_scripts':['kindel=kindel.cli:main']}
      )
