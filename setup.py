using_setuptools = True
try:
    from setuptools import setup
except ImportError:
    using_setuptools = False
    from distutils.core import setup
import configuration


if __name__ == "__main__":
    if using_setuptools:
        setup(
            name=configuration.name,
            version=configuration.version,
            packages=configuration.packages,
            url=configuration.url,
            license=configuration.license,
            author=configuration.author,
            author_email=configuration.author_email,
            description=configuration.description,
            keywords = configuration.keywords,
            install_requires = configuration.install_requires,
            download_url=configuration.download_url
        )

    else:
        setup(
            name=configuration.name,
            version=configuration.version,
            packages=configuration.packages,
            url=configuration.url,
            license=configuration.license,
            author=configuration.author,
            author_email=configuration.author_email,
            description=configuration.description,
            keywords = configuration.keywords,
            download_url=configuration.download_url
        )
