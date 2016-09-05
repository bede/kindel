from setuptools import setup

setup(name='kindel',
      version='0.1.2',
      description='Indel-aware offline consensus calling for DNA alignments in SAM/BAM format',
      url='https://github.com/bede/kindel',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['kindel'],
      zip_safe=True,
      install_requires=['argh','biopython','simplesam'],
      entry_points = {'console_scripts':['kindel=kindel.cli:main']}
      )