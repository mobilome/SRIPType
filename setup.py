from setuptools import setup, find_packages

setup(
    name="scgt",
    version="1.0.0",
    packages=find_packages(include=["scgt_modules", "scgt_modules.*"]),
    package_data={
        "scgt_modules": ["data/*"],
    },
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "scgt=scgt_modules.cli:main",
        ],
    },
)
