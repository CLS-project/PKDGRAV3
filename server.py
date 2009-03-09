#
# pkdgrav2 server script
#
# Run this python script from within pkdgrav2 (-script),
# or "achScript" to start the python server.
#
# The client can then connect with:
#
#     import Pyro.core
#     msr = Pyro.core.getProxyForURI("PYRO://ip:port/guid")
#
# Replace the connection URI with what is printed in the script below.
# MSR functions can now be called remotely, as in:
#
#     msr.Load("test.std")
#     msr.DeepestPotential()
#

import msr
import Pyro.core

class remoteMSR(msr.MSR,Pyro.core.ObjBase):
    def __init__(self):
        Pyro.core.ObjBase.__init__(self)
        msr.MSR.__init__(self)

Pyro.core.initServer()
daemon=Pyro.core.Daemon()
uri=daemon.connect(remoteMSR(),"msr")

print "The daemon runs on port:",daemon.port
print "The object's uri is:",uri

daemon.requestLoop()

