from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name = 'PrediProt',
    version = '1.0.16',
    description = 'Protein complex structure predictor',
    url='https://pypi.org/project/PrediProt/',
    long_description = readme(),
    long_description_content_type = 'text/markdown',
    author = 'The Dream Team',
    license = 'MIT',
    packages = ['prediprot','prediprot.imports'],
    install_requires = ['argparse'],
    python_requires='>=3.6',
    scripts=['prediprot/PrediProt', 'prediprot/PrediProt_getfasta', 'prediprot/PrediProt_GUI'],
    include_package_data=True,
    zip_safe=False)
