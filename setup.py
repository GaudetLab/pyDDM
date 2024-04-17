from setuptools import setup, find_packages
setup(name='pyDDM',
      version='1.0',
      author='Sam Berry',
      author_email="sberry@g.harvard.edu",
      description = "Calculate distance difference matrices to compare protein conformational changes",
      py_modules=["pyddm"],
      packages=find_packages(),
      install_requires=[
        'numpy',
        'pandas',
        'scikit-learn',
        'scipy'
      ],
        dependency_links=[
        'git+https://github.com/samberry19/bioviper.git@v1.2.3#egg=bioviper'
      ]
      )
