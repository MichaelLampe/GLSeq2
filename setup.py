using_setuptools = True

from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
import configuration
import os


if __name__ == "__main__":
    for scheme in INSTALL_SCHEMES.values():
        scheme['data'] = scheme['purelib'] + "/glseq"

    scripts = "GLSeqScripts/Main"
    interface = "User Interface"

    setup(
        name=configuration.name,
        version=configuration.version,
        packages=configuration.packages,
        url=configuration.url,
        license=configuration.license,
        author=configuration.author,
        data_files=[
                    (scripts, [scripts + "/" + f for f in os.listdir(scripts)]),
                    (interface, [interface + "/GLSeq2_UI.jar", interface + "/README.txt"]),
                    ],
        author_email=configuration.author_email,
        description=configuration.description,
        keywords = configuration.keywords,
        download_url=configuration.download_url
    )
