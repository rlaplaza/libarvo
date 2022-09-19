from skbuild import setup

URL = "https://github.com/rlaplaza/libarvo"
DESCRIPTION = "Python library for arvo"
LONG_DESCRIPTION = f"""\
{DESCRIPTION}. For more information, see the [project repository]({URL}).
"""

setup(
    name="libarvo",
    version="0.1.0",
    author="R. Laplaza",
    author_email="rlaplaza@duck.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    packages=["libarvo"],
    python_requires=">=3.8",
    install_requires=["numpy>=1.20"],
    include_package_data=True,
    cmake_args=["-DSKBUILD=ON"],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
