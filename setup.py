from setuptools import setup, find_packages

setup(
    name='GFHunter',
    version='1.0.1',
    description=(
    	'GFHunter is a long read sequence transcriptome alignment-based fusion genes detection tool.'
    	 ),
    long_description=open('README.md').read(),
    author='Zhenhao Lu',
    author_email='luzhenhao.hit@gmail.com',
    maintainer='Zhenhao Lu',
    maintainer_email='luzhenhao.hit@gmail.com',
    license='MIT License',
    url='https://github.com/luzhenhao-HIT/GFHunter',
    packages=find_packages(),
    entry_points={
    'console_scripts': [
        'GFHunter = GFHunter.main:Parser_set',
        ]
    },
    classifiers=[
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ],
    python_requires='>=3',
    install_requires=[
		'numpy',
        'pysam',
		'scipy',
		'intervaltree',
		'pyabpoa'
    ]
)
