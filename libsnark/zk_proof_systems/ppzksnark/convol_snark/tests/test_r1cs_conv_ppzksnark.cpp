/** @file
 *****************************************************************************
 Test program that exercises the ppzkSNARK (first generator, then
 prover, then verifier) on a synthetic R1CS instance.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <cassert>
#include <cstdio>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/examples/run_r1cs_gg_ppzksnark.hpp>

using namespace libsnark;

template<typename ppT>
void test_r1cs_gg_ppzksnark(size_t kernel_h, size_t input_h)//size_t kernel_w, size_t input_h, size_t input_w)//size_t num_constraints,
                         //size_t input_size)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    std::vector<libff::Fr<ppT>> a;
    for(size_t i=0;i<kernel_h;i++){
        // a.push_back(libff::Fr<ppT>(i+1));
        a.push_back(libff::Fr<ppT>::random_element());
    }
    std::vector<libff::Fr<ppT>> x;
    for(size_t i=0;i<input_h;i++){
        // x.push_back(libff::Fr<ppT>(i+1));
        x.push_back(libff::Fr<ppT>::random_element());
    }

    const bool test_serialization = true;
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_origin_convol_example<libff::Fr<ppT> >(input_h, x, kernel_h, a);
    //r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_image_convol<libff::Fr<ppT> >(kernel_h, kernel_w, input_h,input_w);
    //r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_field_input<libff::Fr<ppT> >(num_constraints, input_size);
    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

int main(int argc, char* argv[])
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    libff::start_profiling();
    if(argc < 2){
		std::cout<<"wrong input"<<std::endl;
	}

    test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>(std::stoi(argv[1]), std::stoi(argv[2]));//, std::stoi(argv[3]), std::stoi(argv[4]));
    //test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>(4, 2);
}

