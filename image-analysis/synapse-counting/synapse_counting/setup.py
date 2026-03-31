from setuptools import setup, find_packages

setup(
    name="synapse_counting",
    version="0.1.4",
    author="Cydric Geyskens",
    author_email="cydric.geyskens@gmail.com",
    description="A custom package for counting synapses",
    long_description=open('README.md').read(),
    packages=find_packages(),
    long_description_content_type="text/markdown",
    python_requires='>=3.11',
    license="MIT",
    install_requires=[
        "scikit-image>=0.25.2", 
        "numpy>=2.2.5", 
        "matplotlib>=3.10.0", 
        "pandas>=2.2.3", 
        "czifile>=2019.7.2.1", 
        "scipy>=1.15.3", 
        "pillow>=11.1.0", 
        "seaborn>=0.13.2", 
        "lxml>=5.4.0"
    ], 
    classifiers=[
    "Operating System :: Unix",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: Microsoft :: Windows",
    ]
)

