from setuptools import setup
from fastmlst import __version__

setup(name='fastmlst',
    version=__version__,
    description='A multi-core tool for multilocus sequence typing of draft genome assemblies using PubMLST typing schemes',
    url='https://github.com/EnzoAndree/FastMLST',
    author='Enzo Guerrero-Araya',
    author_email='biologoenzo@gmail.com',
    license='GPLv3',
    packages=['fastmlst'],
    package_data={'fastmlst': ['bundle/scheme_catalog.json', 'bundle/scheme_catalog_meta.json']},
    install_requires=['tqdm',
                      'pandas',
                      'biopython',
                      'requests',
                      'requests-oauthlib'],
    entry_points={
        'console_scripts': [
            'fastmlst = fastmlst.cli:main'
        ]
    },
    zip_safe=False)
