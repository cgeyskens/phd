from setuptools import setup, find_packages

setup(
    name="synapse_counting",
    version="0.0.2",
    author="Cydric Geyskens",
    author_email="cydric.geyskens@gmail.com",
    description="A custom package for counting synapses",
    long_description=open('README.md').read(),
    packages=find_packages(),
    long_description_content_type="text/markdown",
    python_requires='>=3.11',
    license="MIT",
    install_requires=["scikit-image", "numpy", "matplotlib", "pandas", "czifile", "scipy", "pillow", "seaborn", "lxml"],
    classifiers=[
    "Operating System :: Unix",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: Microsoft :: Windows",
    ]
)

