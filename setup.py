from setuptools import setup, find_packages
 
setup(
    name="scPANDA",
    version="1.0.0",
    author="Chang-Xiao Li, Can Huang, Dong-Sheng Chen*",
    author_email="xtlichangxiao@student.pumc.edu.cn",
    description=" a three-layer PAN-blood single-cell transcriptomic data annotator.",
    long_description="scPANDA is a tool designed to infer cell types in blood single-cell RNA sequencing (scRNA-seq) datasets using a three-layer annotation approach.",
    license="MIT",
    url="https://github.com/Cavalpionn/scPANDA/",
 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Science/Research",
    ],
    packages=setuptools.find_packages(),
    python_requires='>=3.10',
    install_requires=['numpy','scanpy','matplotlib'],
)
