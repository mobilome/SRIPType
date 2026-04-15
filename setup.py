from setuptools import setup, find_packages

setup(
    name="sriptype",
    version="0.1.0",
    packages=find_packages(include=["sriptype_modules", "sriptype_modules.*"]),
    package_data={
        "sriptype_modules": ["data/*"],
    },
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "sriptype=sriptype_modules.cli:main",
        ],
    },
)
