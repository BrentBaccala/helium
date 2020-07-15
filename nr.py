#!/usr/bin/env python
#
# Simple python CORBA server (taken from omniorb examples)

import sys
from omniORB import CORBA, PortableServer
import _GlobalIDL, _GlobalIDL__POA

import requests

class Echo_i (_GlobalIDL__POA.Echo):
    def echoString(self, mesg):
        print("echoString() called with message:", mesg)
        return mesg
    def vecfunc(self, x):
        val = requests.post("http://localhost:4443/vecfunc", json=x).json()
        return val

orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
poa = orb.resolve_initial_references("RootPOA")

ei = Echo_i()
eo = ei._this()

print(orb.object_to_string(eo))

poaManager = poa._get_the_POAManager()
poaManager.activate()

orb.run()
