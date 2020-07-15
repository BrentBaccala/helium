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

VecDoub vecfunc(VecDoub_I &x)
{
  return x;
}

Int mymain(void)
{
  const Int NTOT=20;
  Int i,jd,nph=2;
  Doub frac,ave,vrnce;
  VecDoub data(NTOT);
  Bool check;
  for (i=0;i<NTOT;i++) {
    flmoon(i,nph,jd,frac);
    data[i]=jd;
  }

  newt(data, check, vecfunc);

  cout << "newt = ";
  for (i=0;i<NTOT;i++) {
    std::cout << data[i] << ' ';
  }
  cout << endl;

  return 0;
}

#include "nr.hh"

// This is the object implementation.

class Echo_i : public POA_Echo
{
public:
  inline Echo_i() {}
  virtual ~Echo_i() {}
  virtual char* echoString(const char* mesg);
};


char* Echo_i::echoString(const char* mesg)
{
  // Memory management rules say we must return a newly allocated
  // string.
  return CORBA::string_dup(mesg);
}


//////////////////////////////////////////////////////////////////////

// This function acts as a client to the object.

static void hello(Echo_ptr e)
{
  if( CORBA::is_nil(e) ) {
    cerr << "hello: The object reference is nil!" << endl;
    return;
  }

  CORBA::String_var src = (const char*) "Hello!";
  // String literals are (char*) rather than (const char*) on some
  // old compilers.  Thus it is essential to cast to (const char*)
  // here to ensure that the string is copied, so that the
  // CORBA::String_var does not attempt to 'delete' the string
  // literal.

  CORBA::String_var dest = e->echoString(src);

  cout << "I said, \"" << (char*)src << "\"." << endl
       << "The Echo object replied, \"" << (char*)dest <<"\"." << endl;
}

int
colocated_main(int argc, char **argv)
{
  CORBA::ORB_ptr orb = CORBA::ORB_init(argc, argv, "omniORB4");

  CORBA::Object_var       obj = orb->resolve_initial_references("RootPOA");
  PortableServer::POA_var poa = PortableServer::POA::_narrow(obj);

  PortableServer::Servant_var<Echo_i> myecho = new Echo_i();
  PortableServer::ObjectId_var myechoid = poa->activate_object(myecho);

  Echo_var myechoref = myecho->_this();

  PortableServer::POAManager_var pman = poa->the_POAManager();
  pman->activate();

  hello(myechoref);

  orb->destroy();
  return 0;
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

    Echo_var echoref = Echo::_narrow(obj);

    if (CORBA::is_nil(echoref)) {
      cerr << "Can't narrow reference to type Echo (or it was nil)." << endl;
      return 1;
    }

    for (CORBA::ULong count=0; count<10; count++)
      hello(echoref);

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
