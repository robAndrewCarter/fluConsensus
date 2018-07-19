from setuptools import setup

setup(
        name='fluConsensus',
        version='0.5',
        long_description=__doc__,
        packages=['fluConsensus'],
        include_package_data=True,
        package_data={
            'fluConsensus':['external/*']
        },
        zip_safe=False,
        scripts=['scripts/fluConsensus'],
        install_requires=['Flask']
)

