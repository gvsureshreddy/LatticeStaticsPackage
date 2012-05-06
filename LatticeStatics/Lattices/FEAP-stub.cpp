#include <iostream>
#include <cstdlib>

using std::cerr;

void error_and_exit()
{
   cerr << "Error FEAP object is not designed for stand-alone use.  Link with Feap8_3.a\n";
   exit(-1);
}

extern "C" void bfbfeap_main_(char const* const ffin, int* ffinlen,
                              int* const numDOFperNode, int* const numSpcDIM,
                              int* const numNodeInMesh, int* const numEqns)
{error_and_exit();}
extern "C" void bfbfeap_get_eqn_id_(int* bfb_id) {error_and_exit();}
extern "C" void bfbfeap_get_eqn_bc_(int* bfb_bc) {error_and_exit();}
extern "C" void bfbfeap_get_nodal_solution_(double* bfb_u) {error_and_exit();}
extern "C" void bfbfeap_set_nodal_solution_(double* bfb_u) {error_and_exit();}
extern "C" void bfbfeap_get_nodal_coords_(double* bfb_x) {error_and_exit();}
extern "C" void bfbfeap_get_potential_energy_(double* bfb_epl) {error_and_exit();}
extern "C" void bfbfeap_get_reduced_residual_(double* bfb_rd) {error_and_exit();}
extern "C" void bfbfeap_get_reduced_tang_(double* bfb_tang) {error_and_exit();}
extern "C" void bfbfeap_call_ener_() {error_and_exit();}
extern "C" void bfbfeap_call_form_() {error_and_exit();}
extern "C" void bfbfeap_call_tang_() {error_and_exit();}
extern "C" void plstop_() {error_and_exit();}
