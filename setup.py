import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
        'laspy',
        'fiona ~= 1.8.3',
        'rasterio',
        'scipy',
        'h5py',
        'pandas',
        'numpy',
        'simplekml',
	'pyproj',
        'matplotlib']

setuptools.setup(
        name="phoreal",
        version = "3.31.1",
        author = "University of Texas CSR",
        author_email = "phoreal@csr.utexas.edu",
        license = 'UT AUSTIN RESEARCH LICENCE',
        description = "A python package for reading and analysis of ATL03/ATL08 ICESat-2 data",
        long_description = long_description,
        long_description_content_type = "text/markdown",
        url = "https://github.com/icesat-2UT/PhoREAL",
        project_urls ={
            "Bug Tracker": "https://github.com/icesat-2UT/PhoREAL/issues",
            },
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: OS Idependent",
        ],
        packages=['phoreal'],
	python_requires=">=3.6",
        include_package_data=True,
        install_requires=requirements,
)
