import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="edit-distance-clustering",
    version="1.0.0",
    author="Your Name",
    author_email="your_email@example.com",
    description="A package for clustering sequences based on their edit distance",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your_username/edit-distance-clustering",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
