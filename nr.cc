/*
 * Code taken from both Numerical Recipes introduction and omniorb examples.
 *
 * Ubuntu packages needed:    apt install omniidl libomniorb4-dev
 * Build stubs:               omniidl -bpython nr.idl
 * Compile:                   g++ -o nr nr.cc nrSK.cc -lomniORB4 -lomnithread
 */

#include <iostream>
#include <iomanip>
#include "NR3/nr3.h"
#include "NR3/calendar.h"
#include "NR3/moment.h"
#include "NR3/ludcmp.h"
#include "NR3/qrdcmp.h"
#include "NR3/roots_multidim.h"

#include "nr.hh"

Echo_var echoref;

VecDoub iv(long seed)
{
  Echo::VecDoub *retval;

  retval = echoref->iv(seed);

  VecDoub retdata(retval->length());

  for (int i=0; i < retval->length(); i++) {
    retdata[i] = (*retval)[i];
  }

  return retdata;
}

VecDoub vecfunc(VecDoub_I &x)
{
  Echo::VecDoub echodata;
  echodata.length(x.size());
  for (int i=0; i < x.size(); i++) {
    echodata[i] = x[i];
  }

  Echo::VecDoub *retval;

  retval = echoref->vecfunc(echodata);

  VecDoub retdata(retval->length());

  for (int i=0; i < retval->length(); i++) {
    retdata[i] = (*retval)[i];
  }

  return retdata;
}

int main(int argc, char** argv)
{
  try {
    CORBA::ORB_var orb = CORBA::ORB_init(argc, argv);

    if (argc != 2) {
      cerr << "usage:  eg2_clt <object reference>" << endl;
      return 1;
    }

    CORBA::Object_var obj = orb->string_to_object(argv[1]);

    echoref = Echo::_narrow(obj);

    if (CORBA::is_nil(echoref)) {
      cerr << "Can't narrow reference to type Echo (or it was nil)." << endl;
      return 1;
    }

    Bool check;

    VecDoub data = iv(0);

    newt(data, check, vecfunc);

    cout << "newt = ";
    for (int i=0; i<data.size();i++) {
      std::cout << data[i] << ' ';
    }
    cout << endl;

    orb->destroy();
  }
  catch (CORBA::TRANSIENT&) {
    cerr << "Caught system exception TRANSIENT -- unable to contact the "
         << "server." << endl;
  }
  catch (CORBA::SystemException& ex) {
    cerr << "Caught a CORBA::" << ex._name() << endl;
  }
  catch (CORBA::Exception& ex) {
    cerr << "Caught CORBA::Exception: " << ex._name() << endl;
  }
  return 0;
}
