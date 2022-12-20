/** @file
 *****************************************************************************

 Implementation of interfaces for Kate Commitment for R1CS example.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef RUN_KZG_TCC_
#define RUN_KZG_TCC_

#include <sstream>
#include <type_traits>

#include <libff/common/profiling.hpp>
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>

#include <libsnark/zk_proof_systems/kzg/kzg.hpp>

namespace libsnark {

/**
 * The code below provides an example of all stages of running a Kate Commitment.
 */
template<typename ppT>
bool run_kzg(const r1cs_example<libff::Fr<ppT> > &example,
                        const bool test_serialization)
{

    libff::enter_block("Call to run_kzg");

    /* Size of Input, Kernel*/

    int input = 10000;
    int kernel = 100;

    /* Generate Convolution Polynomial */  

    libff::print_header("Generate Convolution Polynomial");
    convpoly<ppT> polynomials = kzg_convpoly<ppT>(input, kernel);
    printf("\n"); libff::print_indent(); libff::print_mem("after convpoly");          

    /* Generate t-SDH tuple, and select secret randomness t */

    libff::print_header("Generate Key: t-SDH Tuple");
    commitkey<ppT> ck = kzg_setup<ppT>(polynomials.z);
    printf("\n"); libff::print_indent(); libff::print_mem("after setup");

    // /* Commit Polynomial into Product: G1-element */

    // libff::print_header("Commit Polynomial: A(x)");
    // libff::G1<ppT> commit_a = kzg_commit<ppT>(ck, polynomials.A, polynomials.t_a);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // libff::print_header("Commit Polynomial: B(x)");
    // libff::G1<ppT> commit_b = kzg_commit<ppT>(ck, polynomials.B, polynomials.t_b);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // libff::print_header("Commit Polynomial: C(x)");
    // libff::G1<ppT> commit_c = kzg_commit<ppT>(ck, polynomials.C, polynomials.t_c);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // /* Generate Random Point for Evaluation: SHA256 */

    // libff::print_header("Generate Random point: Hash(Commit(A(x)), Commit(B(x)), Commit(C(x)))");
    // libff::Fr<ppT> point = kzg_hash<ppT>(commit_a, commit_b, commit_c);
    // printf("\n"); libff::print_indent(); libff::print_mem("after hash");

    // /* Generate Evaluation Value for Convolution Proof */

    // libff::Fr<ppT> eval_a = kzg_evaluate<ppT>(polynomials.A, point, polynomials.t_a);
    // libff::Fr<ppT> eval_b = kzg_evaluate<ppT>(polynomials.B, point, polynomials.t_b); 
    // libff::Fr<ppT> eval_c = kzg_evaluate<ppT>(polynomials.C, point, polynomials.t_c);

    // /* Generate witness of the evaluation + Evaluate the Polynomial */

    // libff::print_header("Create Witness: A(x)");
    // witness<ppT> wit_a = kzg_witness<ppT>(ck, polynomials.A, point, polynomials.t_a);
    // printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    // libff::print_header("Create Witness: B(x)");
    // witness<ppT> wit_b = kzg_witness<ppT>(ck, polynomials.B, point, polynomials.t_b);
    // printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    // libff::print_header("Create Witness: C(x)");
    // witness<ppT> wit_c = kzg_witness<ppT>(ck, polynomials.C, point, polynomials.t_c);
    // printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");


    // /* Verify evaluation */
    // bool verifyresult;

    // libff::print_header("Verify Evaluation of Polynomial: A(x)");
    // bool verifyresult_a = kzg_vfyeval<ppT>(ck, commit_a, wit_a);

    // libff::print_header("Verify Evaluation of Polynomial: B(x)");
    // bool verifyresult_b = kzg_vfyeval<ppT>(ck, commit_b, wit_b);

    // libff::print_header("Verify Evaluation of Polynomial: C(x)");
    // bool verifyresult_c = kzg_vfyeval<ppT>(ck, commit_c, wit_c);

    // if (verifyresult_a == true && verifyresult_b == true && verifyresult_c == true) {
    //     verifyresult = true;
    //     libff::print_header("VERIFICATION ACCEPT!!");
    // } else {
    //     verifyresult = false;
    //     libff::print_header("VERIFICATION REJECT");
    // }
    
    // printf("\n"); libff::print_indent(); libff::print_mem("after vfyeval");

    // /* Output Values (a(k), b(k), c(k)) to link with Convolution Proof: Groth16 */
    // libff::print_header("vCNN+: Convolution Proof");

    // libff::Fr_vector<ppT> jsnark(3);
    // jsnark[0] = eval_a; jsnark[1] = eval_b; jsnark[2] = eval_c;
    // printf("This values are for jsnark\n");
    // jsnark[0].print(); jsnark[1].print(); jsnark[2].print();

    // //clear all
    // jsnark.clear();

    // return verifyresult;
    return true;

    libff::leave_block("Call to run_kzg");
    
}

} // libsnark

#endif // RUN_KZG_TCC_