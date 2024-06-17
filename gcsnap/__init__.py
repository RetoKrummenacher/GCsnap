"""
Gcsnap.py
"""


# import commands for all modules
# __all__ = ["file1", "file2", "file3"]




# The project.toml file you're seeing is a configuration file for a Python project. It's a newer standard that aims to replace setup.py, requirements.txt, MANIFEST.in, and other build-related files. It's part of PEP 621 and PEP 518.

# When you run python setup.py install, Python's setuptools will look for a setup.py file by default. If it doesn't find one, it will look for a pyproject.toml file (not project.toml). If it finds a pyproject.toml file, it will use the build system specified in that file to build and install the project.

# In your case, the build system is specified as setuptools.build_meta, which is a backend that supports PEP 621 and PEP 518. This means that setuptools will use the information in the project.toml file to build and install the project.

# The [project.scripts] section in the project.toml file specifies that the GCsnap command should run the main function in the gcsnap.GCsnap module. This is equivalent to the entry_points argument in a setup.py file.

# The [tool.setuptools] section specifies that the package includes the gcsnap module. This is equivalent to the packages argument in a setup.py file.

# The [project.dependencies] section specifies the dependencies of the project. This is equivalent to the install_requires argument in a setup.py file.

# So, when you run python setup.py install, setuptools will use the information in the project.toml file to install the gcsnap module, its dependencies, and the GCsnap command.

# Please note that the project.toml file should be named pyproject.toml according to PEP 518. If it's not, you might need to rename it.