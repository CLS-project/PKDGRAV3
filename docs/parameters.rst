==========
Parameters
==========

This section describes all of the parameters that can be set in the input parameter file.
This file is actually parsed by a Python interpreter, so any valid Python expression can
be used to construct the parameter values.

During normal "simulation mode", the code will check all variables that were set in the
parameter files, and if it encounters one that it doesn't recognise, then it will issue
and error and abort the program. This is to catch misspelled parameters.

If you need to use intermediate variables to calculate parameter values, then prefix them
with an underscore; the code ignores any variables that start with an underscore.

.. make_parameters:: ../parameters.toml

