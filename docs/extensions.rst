==========
Extensions
==========

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
