import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="openplaning",
    version="0.1.2",
    author="Esteban L. Castro-Feliciano",
    author_email="ecastro@crown-hydro.com",
    description="Hydrodynamic evaluation of planing hulls based on the Savitsky empirical methods.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elcf/python-openplaning",
    project_urls={
        "Bug Tracker": "https://github.com/elcf/python-openplaning/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages = ['openplaning'],
    python_requires=">=3.6",
    install_requires=[
        'numpy',
        'scipy',
        'ndmath',
    ],
    package_data={'openplaning':['tables/*.csv']},
)
