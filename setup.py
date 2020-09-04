from setuptools import setup

setup(name='fastmlst',
    version='0.0.3',
    description='A Fast Multilocus Sequence Typing scan against PubMLST typing schemes',
    url='https://github.com/EnzoAndree/FastMLST',
    author='Enzo Guerrero-Araya',
    author_email='biologoenzo@gmail.com',
    license='GPLv3',
    packages=['fastmlst'],
    scripts=['bin/fastmlst'],
    zip_safe=False)
