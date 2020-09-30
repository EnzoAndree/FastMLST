from setuptools import setup

setup(name='fastmlst',
    version='0.0.8b',
    description='A Fast Multilocus Sequence Typing scan against PubMLST typing schemes',
    url='https://github.com/EnzoAndree/FastMLST',
    author='Enzo Guerrero-Araya',
    author_email='biologoenzo@gmail.com',
    license='GPLv3',
    packages=['fastmlst'],
    install_requires=['tqdm', 'pandas', 'biopython'],
    scripts=['bin/fastmlst'],
    zip_safe=False)
