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

#include <libsnark/common/default_types/r1cs_vcnn_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/examples/run_r1cs_vcnn_ppzksnark.hpp>

using namespace libsnark;

template<typename ppT>
void test_r1cs_vcnn_ppzksnark(size_t input_h, size_t kernel_h)//size_t kernel_w, size_t input_h, size_t input_w)//size_t num_constraints,
                         //size_t input_size)
{
    libff::print_header("(enter) Test R1CS VCNN-ppzkSNARK");

    libff::print_header("Generate Convolution Example");

    //Quantization Depth: 값을 2^bit 이하로 제한시킨다.
    int quanti = 8;
    int quantint = 1;

    for(size_t i = 0; i < quanti; i++) {
        quantint = quantint * 2;
    }
    printf("%d\n", quantint);
    
    std::vector<libff::Fr<ppT>> a;
    libff::Fr<ppT> temp = libff::Fr<ppT>::one();
    for(size_t i=0;i<kernel_h;i++){
        int random = rand();
        random = random % quantint;
        for(size_t j=0;j<random;j++) {
            temp = temp + libff::Fr<ppT>::one();
        }
        // a.push_back(libff::Fr<ppT>(i+1));
        // a.push_back(libff::Fr<ppT>::random_element());
        // temp.print();
        a.push_back(temp);
        temp = libff::Fr<ppT>::one();
    }

    std::vector<libff::Fr<ppT>> x;
    libff::Fr<ppT> temp2 = libff::Fr<ppT>::one();
    for(size_t i=0;i<input_h;i++){
        int random = rand();
        random = random % quantint;
        for(size_t j=0;j<random;j++) {
            temp2 = temp2 + libff::Fr<ppT>::one();
        }
        // temp2.print();
        x.push_back(temp2);
        temp2 = libff::Fr<ppT>::one();
        // x.push_back(libff::Fr<ppT>(i+1));
        // x.push_back(libff::Fr<ppT>::random_element());
    }

    const bool test_serialization = true;
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_origin_convol_example<libff::Fr<ppT> >(input_h, x, kernel_h, a);
    //r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_image_convol<libff::Fr<ppT> >(kernel_h, kernel_w, input_h,input_w);
    //r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_field_input<libff::Fr<ppT> >(num_constraints, input_size);
    const bool bit = run_r1cs_vcnn_ppzksnark<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test R1CS VCNN-ppzkSNARK");
}

int main(int argc, char* argv[])
{
    default_r1cs_vcnn_ppzksnark_pp::init_public_params();
    libff::start_profiling();
    if(argc < 2){
		std::cout<<"wrong input"<<std::endl;
	}

    test_r1cs_vcnn_ppzksnark<default_r1cs_vcnn_ppzksnark_pp>(std::stoi(argv[1]), std::stoi(argv[2]));//, std::stoi(argv[3]), std::stoi(argv[4]));
    //test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>(4, 2);
}
