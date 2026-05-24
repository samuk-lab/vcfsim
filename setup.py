from setuptools import setup

requirements = [
    'numpy',
    'msprime',
    'tskit',
]

setup(
    name='vcfsim',
    version='1.0.29.alpha',
    packages=['vcfsim'],
    entry_points={
        'console_scripts': [
            'vcfsim=vcfsim.__main__:main'
        ]
    },
    url='https://github.com/samuk-lab/vcfsim',
    license='MIT',
    author='Paimon Goulart',
    author_email='paimongoulart@gmail.com',
    description='Coalescent all-sites VCF simulator with configurable missingness',
    install_requires=requirements,
    keywords='pixy vcf simulation coalescent population-genetics',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ]
)
