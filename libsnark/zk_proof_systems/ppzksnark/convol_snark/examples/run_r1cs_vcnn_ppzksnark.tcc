/** @file
 *****************************************************************************

 Implementation of functionality that runs the R1CS GG-ppzkSNARK for
 a given R1CS example.

 See run_r1cs_gg_ppzksnark.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef RUN_R1CS_VCNN_PPZKSNARK_TCC_
#define RUN_R1CS_VCNN_PPZKSNARK_TCC_

#include <sstream>
#include <type_traits>

#include <libff/common/profiling.hpp>

#include <libsnark/common/default_types/r1cs_vcnn_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/r1cs_conv_ppzksnark.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/r1cs_legosnark.hpp>

#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/kzg/kzg.hpp>


namespace libsnark {
template<typename ppT>
bool run_r1cs_vcnn_ppzksnark(const r1cs_example<libff::Fr<ppT> > &example, 
                             const bool test_serialization)
{
    // libff::enter_block("Call to evaluation proof");

    // //Eval - SNARK

    // libff::print_header("Evaluation Proof: Kate-Commitment (KZG10)");
   
    // /* Generate Evaluation Value for Convolution Proof */
    // size_t a_point = example.constraint_system.convol_input_height[0];
    // size_t b_point = example.constraint_system.convol_kernel_height[0];
    // size_t c_point = example.auxiliary_input.size();
    // size_t depth = 3;
    // size_t d_point = b_point / depth;

    // printf("A size: %d, B size: %d, C size: %d, D: %d\n", a_point, b_point, c_point, d_point);

    // libff::Fr_vector<ppT> poly_a(example.primary_input.begin(), 
    //                              example.primary_input.begin()+(a_point));
    // libff::Fr_vector<ppT> poly_b(example.primary_input.begin()+(a_point), 
    //                              example.primary_input.end());
    // libff::Fr_vector<ppT> poly_c(example.auxiliary_input.begin(), 
    //                              example.auxiliary_input.end());

    // /* Generate t-SDH tuple, and select secret randomness t */

    // libff::print_header("Generate Key: t-SDH Tuple");
    // commitkey<ppT> ck = kzg_setup<ppT>(c_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after setup");

    // /* Commit Polynomial into Product: G1-element */

    // libff::print_header("Commit Polynomial: A(x)");
    // libff::G1<ppT> commit_a = kzg_commit<ppT>(ck, poly_a, a_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // libff::print_header("Commit Polynomial: B(x)");
    // libff::G1<ppT> commit_b = kzg_commit<ppT>(ck, poly_b, b_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // libff::print_header("Commit Polynomial: C(x)");
    // libff::G1<ppT> commit_c = kzg_commit<ppT>(ck, poly_c, c_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    // /* Generate Random Point for Evaluation: SHA256 */

    // libff::print_header("Generate Random point: Hash(Commit(C(x)))");
    // libff::Fr<ppT> point = kzg_hash<ppT>(commit_c);
    // printf("\n"); libff::print_indent(); libff::print_mem("after hash");

    // /* Generate Evaluation Value for Convolution Proof */

    // libff::print_header("Evaluate Polynomial: A(point), B(point), C(point)");
    // libff::Fr<ppT> eval_a = kzg_evaluate<ppT>(poly_a, point, a_point);
    // libff::Fr<ppT> eval_b = kzg_evaluate<ppT>(poly_b, point, b_point); 
    // libff::Fr<ppT> eval_c = kzg_evaluate<ppT>(poly_c, point, c_point);

    // /* Generate witness of the evaluation + Evaluate the Polynomial */

    // libff::print_header("Create Witness: A(x)");
    // witness<ppT> wit_a = kzg_witness<ppT>(ck, poly_a, point, a_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    // libff::print_header("Create Witness: B(x)");
    // witness<ppT> wit_b = kzg_witness<ppT>(ck, poly_b, point, b_point);
    // printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    // libff::print_header("Create Witness: C(x)");
    // witness<ppT> wit_c = kzg_witness<ppT>(ck, poly_c, point, c_point);
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
    //     libff::print_header("VERIFICATION ACCEPT");
    // } else {
    //     verifyresult = false;
    //     libff::print_header("VERIFICATION REJECT");
    // }

    // //Conv. SNARK
    // /* Output Values (a(k), b(k), c(k)) to link with Convolution Proof(jsnark): Groth16 */
    // libff::Fr_vector<ppT> jsnark(3);
    // jsnark[0] = eval_a; jsnark[1] = eval_b; jsnark[2] = eval_c;
    // printf("\nThis values are for jsnark\n");
    // jsnark[0].print(); jsnark[1].print(); jsnark[2].print();

    // libff::leave_block("Call to evaluation proof");

    // return verifyresult;

    libff::enter_block("Call to evaluation proof");

    //Eval - SNARK: Version for multi_depth.

    libff::print_header("Evaluation Proof: Kate-Commitment (KZG10)");
   
    /* Generate Evaluation Value for Convolution Proof */
    size_t a_point = example.constraint_system.convol_input_height[0];
    size_t b_point = example.constraint_system.convol_kernel_height[0];
    size_t c_point = example.auxiliary_input.size();
    size_t depth = example.constraint_system.num_convol - 1;
    size_t d_point = b_point / depth; //size of kernel
    size_t e_point = c_point / depth;

    printf("A size: %d, B size: %d, C size: %d, D: %d, E: %d\n", a_point, b_point, c_point, d_point, e_point);

    libff::Fr_vector<ppT> poly_a(example.primary_input.begin(), 
                                 example.primary_input.begin()+(a_point));
    libff::Fr_vector<ppT> poly_b(example.primary_input.begin()+(a_point), 
                                 example.primary_input.end());
    libff::Fr_vector<ppT> poly_c(example.auxiliary_input.begin(), 
                                 example.auxiliary_input.end());
    // libff::Fr_vector<ppT> poly_d(example.auxiliary_input.begin()+(c_point), 
    //                              example.auxiliary_input.end());                             

    /* Generate t-SDH tuple, and select secret randomness t */

    libff::print_header("Generate Key: t-SDH Tuple");
    commitkey<ppT> ck = kzg_setup<ppT>(c_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after setup");

    /* Commit Polynomial into Product: G1-element */

    libff::print_header("Commit Polynomial: A(x)");
    libff::G1<ppT> commit_a = kzg_commit<ppT>(ck, poly_a, a_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    libff::print_header("Commit Polynomial: B(x)");
    libff::G1<ppT> commit_b = kzg_commit<ppT>(ck, poly_b, b_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    libff::print_header("Commit Polynomial: C(x)");

    //before committing poly C, we need to rearrange the C_polys for Convolution Groth-16 proof.
    size_t final_c_point = e_point + d_point * (depth - 1);
    libff::Fr_vector<ppT> final_poly_c; final_poly_c.reserve(final_c_point);
    
    for(size_t i = 0; i < depth; i++) {
        for(size_t j = 0; j < e_point; j++) {
            // printf("poly_c: %d\n", j + (d_point * i));
            // printf("poly_e: %d\n", j + (e_point * i));
            final_poly_c[j + (d_point * i)] = final_poly_c[j + d_point * i] + poly_c[j + (e_point * i)];
        }
    }

    // //TEST.
    // int temp = 3;
    // int temp2 = 4;
    // int temp3 = 2;

    // int tempoly[3]; tempoly[0] = 1; tempoly[1] = 2; tempoly[2] = 3;
    // int tempoly2[6]; tempoly2[0] = 1; tempoly2[1] = 2; tempoly2[2] = 3; tempoly2[3] = 4; tempoly2[4] = 5; tempoly2[5] = 6;

    // int tempoly3[12]; tempoly3[0] = 1; tempoly3[1] = 4; tempoly3[2] = 7; tempoly3[3] = 6; tempoly3[4] = 3; tempoly3[5] = 10;
    // tempoly3[6] = 17; tempoly3[7] = 12; tempoly3[8] = 5; tempoly3[9] = 16; tempoly3[10] = 27; tempoly3[11] = 18;
    
    // int result[8] = {0, };

    // for(size_t i = 0; i < 3; i++) {
    //     for(size_t j = 0; j < temp2; j++) {
    //         printf("A: %d, B: %d\n", j + (temp3 * i), j + (temp2 * i));
    //         printf("A: %d, B: %d\n", result[j + (temp3 * i)],tempoly3[j + (temp2 * i)]);
    //         result[j + (temp3 * i)] += tempoly3[j + (temp2 * i)];
    //     }
    // }
    // for(size_t i = 0; i < 8; i++) {
    //     printf("%d ", result[i]);
    //     printf("\n");
    // }


    // for(size_t i = 0; i < a_point; i++) {
    //     poly_a[i].print();
    // }
    // printf("%d\n");
    // libff::G1_vector<ppT> commits_c; commits_c.reserve(depth);
    // for(size_t i = 0; i < depth; i++) {
    //     libff::Fr_vector<ppT> poly_c_depth(poly_c.begin()+(depth * a_point),
    //                                        poly_c.begin()+((depth+1) * a_point));

    //     commits_c[i] = kzg_commit<ppT>(ck, poly_c_depth, d_point);
    //     poly_c_depth.clear();
    // }
    libff::G1<ppT> commit_c = kzg_commit<ppT>(ck, final_poly_c, final_c_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    /* Generate Random Point for Evaluation: SHA256 */

    libff::print_header("Generate Random point: Hash(Commit(C(x)))");
    libff::Fr<ppT> point = kzg_hash<ppT>(commit_c);
    printf("\n"); libff::print_indent(); libff::print_mem("after hash");

    /* Generate Evaluation Value for Convolution Proof */

    libff::print_header("Evaluate Polynomial: A(point), B(point), C(point)");
    libff::Fr<ppT> eval_a = kzg_evaluate<ppT>(poly_a, point, a_point);
    libff::Fr<ppT> eval_b = kzg_evaluate<ppT>(poly_b, point, b_point); 
    libff::Fr<ppT> eval_c = kzg_evaluate<ppT>(final_poly_c, point, final_c_point);

    /* Generate witness of the evaluation + Evaluate the Polynomial */

    libff::print_header("Create Witness: A(x)");
    witness<ppT> wit_a = kzg_witness<ppT>(ck, poly_a, point, a_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    libff::print_header("Create Witness: B(x)");
    witness<ppT> wit_b = kzg_witness<ppT>(ck, poly_b, point, b_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    libff::print_header("Create Witness: C(x)");
    witness<ppT> wit_c = kzg_witness<ppT>(ck, final_poly_c, point, final_c_point);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");


    /* Verify evaluation */
    bool verifyresult;

    libff::print_header("Verify Evaluation of Polynomial: A(x)");
    bool verifyresult_a = kzg_vfyeval<ppT>(ck, commit_a, wit_a);

    libff::print_header("Verify Evaluation of Polynomial: B(x)");
    bool verifyresult_b = kzg_vfyeval<ppT>(ck, commit_b, wit_b);

    libff::print_header("Verify Evaluation of Polynomial: C(x)");
    bool verifyresult_c = kzg_vfyeval<ppT>(ck, commit_c, wit_c);

    if (verifyresult_a == true && verifyresult_b == true && verifyresult_c == true) {
        verifyresult = true;
        libff::print_header("VERIFICATION ACCEPT");
    } else {
        verifyresult = false;
        libff::print_header("VERIFICATION REJECT");
    }

    //Conv. SNARK
    /* Output Values (a(k), b(k), c(k)) to link with Convolution Proof(jsnark): Groth16 */
    libff::Fr_vector<ppT> jsnark(3);
    jsnark[0] = eval_a; jsnark[1] = eval_b; jsnark[2] = eval_c;
    printf("\nThis values are for jsnark\n");
    jsnark[0].print(); jsnark[1].print(); jsnark[2].print();

    libff::leave_block("Call to evaluation proof");

    return verifyresult;

}

} // libsnark

#endif // RUN_R1CS_VCNN_PPZKSNARK_TCC_