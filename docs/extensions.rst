==========
Extensions
==========

.. index:: accuracy
  :name: accuracy-module

Accuracy
--------

The accuracy module is designed to alter the accuracy parameters in a dynmaic way.
For cosmological runs it is **strongly** recommended to set `dTheta` and `nReplicas`
in the following way::

    from accuracy import classic_theta_switch,classic_replicas_switch

    dTheta          = classic_theta_switch()
    nReplicas       = classic_replicas_switch()

.. automodule:: accuracy
   :members:


PKDGRAV Interface
-----------------

Instead of simulation mode, you can enter analysis mode by importing the ``PKDGRAV``
module. Once imported you are then responsible for calling the correct methods
to achieve your analysis tasks. For example, the following script can be used
to measure the power spectrum of a periodic simulation output.

.. literalinclude:: ../tools/measurepk.py

The complete list of methods is as follows.

.. automodule:: PKDGRAV
   :members:

Cosmology Interface
-------------------

The classes SimpleCosmology and ClassCosmology expose the interface for determining
cosmological values used in the code.

The complete list of classes and methods is as follows.

.. automodule:: cosmology
   :members:

