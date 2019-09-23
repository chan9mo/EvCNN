/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libsnark/gadgetlib2/adapters.hpp>
#include <libsnark/gadgetlib2/integration.hpp>

namespace libsnark {

linear_combination<libff::Fr<libff::default_ec_pp> > convert_gadgetlib2_linear_combination(const gadgetlib2::GadgetLibAdapter::linear_combination_t &lc)
{
    typedef libff::Fr<libff::default_ec_pp> FieldT;
    typedef gadgetlib2::GadgetLibAdapter GLA;

    linear_combination<FieldT> result = lc.second * variable<FieldT>(0);
    for (const GLA::linear_term_t &lt : lc.first)
    {
        //std::cout<<"fir : "<<lt.first <<" sec : "<<lt.second<<std::endl;
        result = result + lt.second * variable<FieldT>(lt.first+1);
    }

    return result;
}

/*
linear_combination<libff::Fr<libff::default_ec_pp> > convert_gadgetlib2_linear_combination(const gadgetlib2::GadgetLibAdapter::linear_combination_t &lc, bool convol)
{
    typedef libff::Fr<libff::default_ec_pp> FieldT;
    typedef gadgetlib2::GadgetLibAdapter GLA;

    //linear_combination<FieldT> result(0); //= lc.second * variable<FieldT>(0);
    bool first = true;
    linear_combination<FieldT> result();
    for (const GLA::linear_term_t &lt : lc.first)
    {
        if(first){
            result = lt.second * variable<FieldT>(lt.first+1);
        }
        else{
        //std::cout<<"fir : "<<lt.first <<" sec : "<<lt.second.as_ulong()<<std::endl;
        result = result + lt.second * variable<FieldT>(lt.first+1);
        }
    }

    return result;
}
*/

r1cs_constraint_system<libff::Fr<libff::default_ec_pp> > get_constraint_system_from_gadgetlib2(const gadgetlib2::Protoboard &pb)
{
    typedef libff::Fr<libff::default_ec_pp> FieldT;
    typedef gadgetlib2::GadgetLibAdapter GLA;

    r1cs_constraint_system<FieldT> result;
    const GLA adapter;

    GLA::protoboard_t converted_pb = adapter.convert(pb);
    for (const GLA::constraint_t &constr : converted_pb.first)
    {
        result.constraints.emplace_back(r1cs_constraint<FieldT>(convert_gadgetlib2_linear_combination(std::get<0>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<1>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<2>(constr))));
    }
    //The number of variables is the highest index created.
    //TODO: If there are multiple protoboards, or variables not assigned to a protoboard, then getNextFreeIndex() is *not* the number of variables! See also in get_variable_assignment_from_gadgetlib2.
    const size_t num_variables = GLA::getNextFreeIndex();
    result.primary_input_size = pb.numInputs();
    result.auxiliary_input_size = num_variables - pb.numInputs();
    return result;
}

r1cs_constraint_system<libff::Fr<libff::default_ec_pp> > get_constraint_convol_system_from_gadgetlib2(const gadgetlib2::Protoboard &pb)
{
    typedef libff::Fr<libff::default_ec_pp> FieldT;
    typedef gadgetlib2::GadgetLibAdapter GLA;

    r1cs_constraint_system<FieldT> result;
    const GLA adapter;

    GLA::protoboard_t converted_pb = adapter.convert(pb);
    for (const GLA::constraint_t &constr : converted_pb.first)
    {
        /*
        std::cout<<"lin1=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<0>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"lin2=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<1>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"lin3=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<2>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"lin4=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<3>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"lin5=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<4>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"lin6=>";
        for(GLA::linear_term_t lt : std::get<0>(std::get<5>(constr))){
            std::cout<<"idx : "<<std::get<0>(lt)<<", value : "<<(std::get<1>(lt)).as_ulong()<<std::endl;
        }
        std::cout<<std::endl;
        */
        result.constraints.emplace_back(r1cs_constraint<FieldT>(convert_gadgetlib2_linear_combination(std::get<0>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<1>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<2>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<3>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<4>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<5>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<6>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<7>(constr)),
                                                                convert_gadgetlib2_linear_combination(std::get<8>(constr))));
    }
    //The number of variables is the highest index created.
    //TODO: If there are multiple protoboards, or variables not assigned to a protoboard, then getNextFreeIndex() is *not* the number of variables! See also in get_variable_assignment_from_gadgetlib2.
    const size_t num_variables = GLA::getNextFreeIndex();
    std::cout<<"num val : "<<num_variables<<std::endl;
    std::cout<<"num input : "<<pb.numInputs()<<std::endl;
    result.primary_input_size = pb.numInputs();
    result.auxiliary_input_size = num_variables - pb.numInputs();
    for(size_t i=0;i<pb.convol_size;i++){
        result.convol_outputs_size.insert(result.convol_outputs_size.end(), pb.convol_outputs_size[i]);
        result.convol_outputs_size2.insert(result.convol_outputs_size2.end(), pb.convol_outputs_size2[i]);
        result.convol_dimensions.insert(result.convol_dimensions.end(), pb.convol_dimensions[i]);
        result.convol_input_height.insert(result.convol_input_height.end(), pb.convol_input_height[i]);
        result.convol_input_width.insert(result.convol_input_width.end(), pb.convol_input_width[i]);
        result.convol_kernel_height.insert(result.convol_kernel_height.end(), pb.convol_kernel_height[i]);
        result.convol_kernel_width.insert(result.convol_kernel_width.end(), pb.convol_kernel_width[i]);
    }
    result.num_convol = pb.convol_size;

    return result;
}

r1cs_variable_assignment<libff::Fr<libff::default_ec_pp> > get_variable_assignment_from_gadgetlib2(const gadgetlib2::Protoboard &pb)
{
    typedef libff::Fr<libff::default_ec_pp> FieldT;
    typedef gadgetlib2::GadgetLibAdapter GLA;

    //The number of variables is the highest index created. This is also the required size for the assignment vector.
    //TODO: If there are multiple protoboards, or variables not assigned to a protoboard, then getNextFreeIndex() is *not* the number of variables! See also in get_constraint_system_from_gadgetlib2.
    const size_t num_vars = GLA::getNextFreeIndex();
    const GLA adapter;
    std::cout<<"num var : "<<num_vars<<std::endl;
    r1cs_variable_assignment<FieldT> result(num_vars, FieldT::zero());
    VariableAssignment assignment = pb.assignment();

    //Go over all assigned values of the protoboard, from every variable-value pair, put the value in the variable.index place of the new assignment.
    for(VariableAssignment::iterator iter = assignment.begin(); iter != assignment.end(); ++iter){
    	result[GLA::getVariableIndex(iter->first)] = adapter.convert(iter->second);
    }
    
//    std::cout<<"var gadget : ";
//    for(size_t i=0;i<result.size();i++){
//        //result[i].print3();
//        std::cout<<result[i].as_ulong()<<"\t";
//    }
//    std::cout<<std::endl;
    

    return result;
}

}
