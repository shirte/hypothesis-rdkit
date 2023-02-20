from setuptools import setup, find_packages


setup(
    name="hypothesis-rdkit",
    version="0.1",
    maintainer="Steffen Hirte",
    maintainer_email="shirte@users.noreply.github.com",
    packages=find_packages(),
    url="https://github.com/shirte/hypothesis-rdkit",
    long_description=open("README.md").read(),
    install_requires=["hypothesis", "rdkit"],
    entry_points={"hypothesis": {"_ = hypothesis_rdkit.hook:_hypothesis_setup_hook"}},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Software Development :: Testing",
        "Framework :: Hypothesis",
    ],
)
