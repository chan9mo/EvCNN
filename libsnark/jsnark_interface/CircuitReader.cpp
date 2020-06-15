/*
 * CircuitReader.cpp
 *
 *      Author: Ahmed Kosba
 */

#include "CircuitReader.hpp"

CircuitReader::CircuitReader(char* arithFilepath, char* inputsFilepath,
		ProtoboardPtr pb) {

	this->pb = pb;

	parseAndEval(arithFilepath, inputsFilepath);
	constructCircuit(arithFilepath);
	mapValuesToProtoboard();

	wireLinearCombinations.clear();
	wireValues.clear();
	variables.clear();
	variableMap.clear();
	zeropMap.clear();
	zeroPwires.clear();
}

void CircuitReader::parseAndEval(char* arithFilepath, char* inputsFilepath) {

	 libff::enter_block("Parsing and Evaluating the circuit");

	ifstream arithfs(arithFilepath, ifstream::in);
	ifstream inputfs(inputsFilepath, ifstream::in);
	string line;

	if (!arithfs.good()) {
		printf("Unable to open circuit file %s \n", arithFilepath);
		exit(-1);
	}


	if (!inputfs.good()) {
		printf("Unable to open input file %s \n", inputsFilepath);
		exit(-1);
	} else {
		char* inputStr;
		while (getline(inputfs, line)) {
			if (line.length() == 0) {
				continue;
			}
			char wireId[80];
			inputStr = new char[line.size()];
			if (2 == sscanf(line.c_str(), "%s %s", wireId, inputStr)) {
				//size_t value;
				//inputStr >> value;
				size_t value;
				istringstream iss_i(inputStr, istringstream::in);
				iss_i >> value;

				wireValues[wireId] = FieldT(value);//readFieldElementFromHex(inputStr);
				// std::cout<<"wireValues["<<wireId<<"] = "<<wireValues[wireId].as_ulong()<<std::endl;
			} else {
				printf("Error in Input\n");
				exit(-1);
			}
			delete[] inputStr;
		}
		inputfs.close();
	}

	char type[200];
	char* inputStr;
	char* outputStr;
	char* stateStr;
	unsigned int numGateInputs, numGateOutputs, numStates;

	char wireId[80];

	FieldT oneElement = FieldT::one();
	FieldT zeroElement = FieldT::zero();
	FieldT negOneElement = FieldT(-1);
	FieldT inverse4 = FieldT(4).inverse();
	FieldT inverse2 = FieldT(2).inverse();
	std::cout<<"inv4 : ";
	inverse4.print3();
	std::cout<<"\ninv2 : ";
	inverse2.print();	


	std::cout<<"-1 : "<<negOneElement.as_ulong()<<std::endl;
	negOneElement.print();
	//std::cout<<"4^-1 : "<<inverse4.as_ulong()<<std::endl;
	//inverse4.print();
	//std::cout<<"200("<<FieldT(200).as_ulong()<<") * "<<"4^-1("<<inverse4.as_ulong()<<") = "<<(inverse4*FieldT(200)).as_ulong()<<std::endl;


	long long evalTime;
	long long begin, end;
	evalTime = 0;

	// Parse the circuit: few lines were imported from Pinocchio's code.

	//getline(arithfs, line);
	while (getline(arithfs, line)) {
		if (line.length() == 0) {
			continue;
		}
		if (line[0] == '#') 
			continue;

		inputStr = new char[line.size()];
		outputStr = new char[line.size()];
		stateStr = new char[line.size()];

		if (1 == sscanf(line.c_str(), "input %s", wireId)) {
			inputWireIds.push_back(wireId);
		} else if (1 == sscanf(line.c_str(), "nizkinput %s", wireId)) {
			nizkWireIds.push_back(wireId);
		} else if (1 == sscanf(line.c_str(), "output %s", wireId)) {
			outputWireIds.push_back(wireId);
			wireUseCounters[wireId]++;
		} else if (1 == sscanf(line.c_str(), "cminput %s", wireId)){
			cminputWireIds.push_back(wireId);
		}
		else if(7 == sscanf(line.c_str(), "%s in %d <%[^>]> out %d <%[^>]> state %d <%[^>]>", type,
						&numGateInputs, inputStr, &numGateOutputs, outputStr, &numStates, stateStr) ){

			
			//std::cout<<"convol 2d 3~~~~\n";
			istringstream iss_i(inputStr, istringstream::in);
			std::vector<FieldT> inValues;
			std::vector<Wire> outWires;
			Wire inWireId;
			//"3 4 5 6 7"
			while (iss_i >> inWireId) {
				//inWireId = 3
				wireUseCounters[inWireId]++;
				inValues.push_back(wireValues[inWireId]);
				//inValue <= wireValue[3] = (w3 =5)
			}
			readIds(outputStr, outWires);
			

			short opcode;
			FieldT constant;
			if(strcmp(type,"convol") == 0 && numStates == 4){
				//std::cout<<"convol 2d 2\n";
				opcode = CONVOL2D_OPCODE;

				istringstream iss_st(stateStr, istringstream::in);
				size_t inputHeight, inputWidth, kernelHeight, kernelWidth;
				iss_st >> inputHeight;
				iss_st >> inputWidth;
				iss_st >> kernelHeight;
				iss_st >> kernelWidth;
				for(size_t j=0; j<(inputWidth+kernelWidth-1);j++){
					for(size_t k=0; k<(inputHeight+kernelHeight-1);k++){
						FieldT y = FieldT::zero();
						for(size_t w=0; w<inputWidth;w++){
							for(size_t h=0;h<inputHeight;h++){
								size_t k1 = j - w;
								size_t k2 = k - h;
								if(( k1 >=0 && k1 < kernelWidth) && ( k2 >=0 && k2< kernelHeight))
								{
											//std::cout<<"w,h,k1,k2 = "<<w<<h<<k1<<k2<<std::endl;
											//std::cout<<"k["<<k1<<"*"<<kernelHeight<<"+"<<k2<<"]("<<(inValues[k1*(kernelHeight)+k2]).as_ulong()<<")"
											//"*x["<<w<<"*"<<inputHeight<<"+"<<h<<"]("<<(inValues[(w*inputHeight+h)]).as_ulong()<<") = "<<
											//(inValues[inputHeight*inputWidth + k1*(kernelHeight)+k2]*inValues[(w*inputHeight+h)]).as_ulong()<<"\t";
									y += inValues[inputHeight*inputWidth + k1*(kernelHeight)+k2]*inValues[(w*inputHeight+h)];
								}
							}
						}
						//std::cout<<"outwire["<<j*(inputHeight+kernelHeight-1) + k<<"] ="<<y.as_ulong()<<"\n";
						//outWire[] = output wire index
						wireValues[outWires[j*(inputHeight+kernelHeight-1) + k]] = y;
					}
				}
				/*
				for(size_t i=0;i<(inputWidth+kernelWidth-1);i++){
					for(size_t j=0;j<(inputHeight+kernelHeight-1);j++)
						std::cout<<outWires[i*(inputHeight+kernelHeight-1) + j]<<"\t";
				}
				std::cout<<std::endl;	
				*/
			}
		}
		else if (5
				== sscanf(line.c_str(), "%s in %u <%[^>]> out %u <%[^>]>", type,
						&numGateInputs, inputStr, &numGateOutputs, outputStr)) {

			istringstream iss_i(inputStr, istringstream::in);
			std::vector<FieldT> inValues;
			std::vector<Wire> outWires;
			Wire inWireId;
			while (iss_i >> inWireId) {
				wireUseCounters[inWireId]++;
				inValues.push_back(wireValues[inWireId]);
			}
			readIds(outputStr, outWires);

			short opcode;
			FieldT constant;
			if (strcmp(type, "add") == 0) {
				opcode = ADD_OPCODE;
			} else if (strcmp(type, "mul") == 0) {
				opcode = MUL_OPCODE;
			} else if (strcmp(type, "xor") == 0) {
				opcode = XOR_OPCODE;
			} else if (strcmp(type, "or") == 0) {
				opcode = OR_OPCODE;
			} else if (strcmp(type, "assert") == 0) {
				wireUseCounters[outWires[0]]++;
				opcode = CONSTRAINT_OPCODE;
			} else if (strcmp(type, "pack") == 0) {
				opcode = PACK_OPCODE;
			} else if (strcmp(type, "zerop") == 0) {
				opcode = NONZEROCHECK_OPCODE;
			} else if (strcmp(type, "split") == 0) {
				opcode = SPLIT_OPCODE;
			} else if (strstr(type, "const-mul-neg-")) {
				opcode = MULCONST_OPCODE;
				char* constStr = type + sizeof("const-mul-neg-") - 1;
				constant = readFieldElementFromHex(constStr) * negOneElement;
			} else if (strstr(type, "const-mul-")) {
				opcode = MULCONST_OPCODE;
				char* constStr = type + sizeof("const-mul-") - 1;
				constant = readFieldElementFromHex(constStr);
			} else if(strstr(type,"convol1d")){
				opcode = CONVOL_OPCODE;
			} else {
				printf("Error: unrecognized line: %s\n", line.c_str());
				assert(0);
			}

			// TODO: separate evaluation from parsing completely to get accurate evaluation cost
			//	 Calling  libff::get_nsec_time(); repetitively as in the old version adds much overhead 
			// TODO 2: change circuit format to enable skipping some lines during evaluation
			//       Not all intermediate wire values need to be computed in this phase
			// TODO 3: change circuit format to make common constants defined once			
	
			//begin = libff::get_nsec_time();
			if (opcode == ADD_OPCODE) {
				FieldT sum;
				for (auto &v : inValues)
					sum += v;
				wireValues[outWires[0]] = sum;
			} else if (opcode == MUL_OPCODE) {
				wireValues[outWires[0]] = inValues[0] * inValues[1];
			} else if (opcode == XOR_OPCODE) {
				wireValues[outWires[0]] =
						(inValues[0] == inValues[1]) ? zeroElement : oneElement;
			} else if (opcode == OR_OPCODE) {
				wireValues[outWires[0]] =
						(inValues[0] == zeroElement
								&& inValues[1] == zeroElement) ?
								zeroElement : oneElement;
			} else if (opcode == NONZEROCHECK_OPCODE) {
				wireValues[outWires[1]] =
						(inValues[0] == zeroElement) ? zeroElement : oneElement;
			} else if (opcode == PACK_OPCODE) {
				FieldT sum, coeff;
				FieldT two = oneElement;
				for (auto &v : inValues) {
					sum += two * v;
					two += two;
				}
				wireValues[outWires[0]] = sum;
			} else if (opcode == SPLIT_OPCODE) {
				int size = outWires.size();
				FElem inVal = inValues[0];
				for (int i = 0; i < size; i++) {
					wireValues[outWires[i]] = inVal.getBit(i, R1P);
				}
			} else if (opcode == MULCONST_OPCODE) {
				wireValues[outWires[0]] = constant * inValues[0];
				//std::cout<<FieldT("16416182153879456416684804308942956316411273300312025757773653139931856371713").as_ulong()<<"*"<<FieldT(4).as_ulong()<<"="<<(FieldT(4)*FieldT("16416182153879456416684804308942956316411273300312025757773653139931856371713")).as_ulong()<<std::endl;
				//std::cout<<"Hex "<<readFieldElementFromHex("244B3AD628E5381F4A3C3448E1210245DE26EE365B4B146CF2E9782EF4000001").as_ulong()<<"*"<<FieldT(4).as_ulong()<<"="<<(FieldT(4)*readFieldElementFromHex("244B3AD628E5381F4A3C3448E1210245DE26EE365B4B146CF2E9782EF4000001")).as_ulong()<<std::endl;
			} else if(opcode == CONVOL_OPCODE){
				size_t num_in = inValues[0].as_ulong();
				size_t num_ker = inValues[num_in+1].as_ulong();
				for(size_t i=0;i<num_in+num_ker-1;i++){
					FieldT y = zeroElement;
					for(size_t k=0;k<num_in;k++){
						for(size_t l=0;l<num_ker;l++){
							if((k+l)==i){
								//std::cout<<"["<<k<<"]["<<l<<"]("<<(inValues[k+1]*inValues[1+num_in+1+l]).as_ulong()<<")";
								y += inValues[k+1] * inValues[1+num_in+l+1];
							}
						}
					}
					//std::cout<<y.as_ulong()<<"\t";
					wireValues[outWires[i]] = y;
				}
				/*
				std::cout<<std::endl;
				for(size_t i=0;i<num_in+num_ker-1;i++){
					std::cout<<outWires[i]<<"\t";
				}
				std::cout<<std::endl;

				*/
			}
			
			
			//end =  libff::get_nsec_time();
			//evalTime += (end - begin);
		} else {
			printf("Error: unrecognized line: %s\n", line.c_str());
			assert(0);
		}
		delete[] inputStr;
		delete[] outputStr;
		delete[] stateStr;
	}
	arithfs.close();

	end = clock();
	// printf("\t Evaluation Done in %lf seconds \n", (double) (evalTime) * 1e-9);
	 libff::leave_block("Parsing and Evaluating the circuit");
}

void CircuitReader::constructCircuit(char* arithFilepath) {

	struct proc_t usage1, usage2;
	cout << "Translating Constraints ... " << endl;
	look_up_our_self(&usage1);
	unsigned int i;

	currentVariableIdx = currentLinearCombinationIdx = 0;

	Wire wire;
	for (auto &wire : inputWireIds) {
		// std::cout<<"input :"<<wire<<std::endl;
		variables.push_back(make_shared<Variable>("input"));
		variableMap[wire] = currentVariableIdx;
		//std::cout<<"variableMap["<<wire<<"]="<<variableMap[wire]<<std::endl;
		currentVariableIdx++;
	}
	for (auto &wire : outputWireIds) {
		variables.push_back(make_shared<Variable>("output"));
		variableMap[wire] = currentVariableIdx;
		//std::cout<<"variableMap["<<wire<<"]="<<variableMap[wire]<<std::endl;

		currentVariableIdx++;
	}
	for (auto &wire : nizkWireIds) {
		variables.push_back(make_shared<Variable>("nizk input"));
		variableMap[wire] = currentVariableIdx;
		currentVariableIdx++;
	}
	for (auto &wire : cminputWireIds) {
		// std::cout<<"cminput :"<<wire<<std::endl;
		variables.push_back(make_shared<Variable>("cminput"));
		variableMap[wire] = currentVariableIdx;
		//std::cout<<"variableMap["<<wire<<"]="<<variableMap[wire]<<std::endl;
		currentVariableIdx++;
		inputWireIds.push_back(wire);
	}
	

	char type[200];
	char* inputStr;
	char* outputStr;
	char* stateStr;
	string line;
	unsigned int numGateInputs, numGateOutputs, numStates;

	ifstream ifs2(arithFilepath, ifstream::in);

	if (!ifs2.good()) {
		printf("Unable to open circuit file:\n");
		exit(5);
	}

	// Parse the circuit: few lines were imported from Pinocchio's code.

	//getline(ifs2, line);

	int lineCount = 0;
	while (getline(ifs2, line)) {
		lineCount++;
//		if (lineCount % 100000 == 0) {
//			printf("At Line:: %d\n", lineCount);
//		}

		if (line.length() == 0) {
			continue;
		}
		if (line[0] == '#')
			continue;
		if (line.compare(0,5,"total")==0)
			continue;

		inputStr = new char[line.size()];
		outputStr = new char[line.size()];
		stateStr = new char[line.size()];
		if (7 == sscanf(line.c_str(), "%s in %d <%[^>]> out %d <%[^>]> state %d <%[^>]>", type,
						&numGateInputs, inputStr, &numGateOutputs, outputStr, &numStates, stateStr)){
			/*
			std::cout<<"################"<<"\ntype : "<<type<<"\n#in : "<<numGateInputs<<"\ninStr : "<<inputStr
			<<"\n#Out : "<<numGateOutputs<<"\noutStr : "<<outputStr<<"\n#state : "<<numStates<<"\nstStr : "<<stateStr
			<<"\n###################"<<std::endl;
			*/
			if(strstr(type, "convol") && numStates == 4){
				//std::cout<<"ok convol2d\n"<<stateStr<<std::endl;
				addConvol2DConstraint(inputStr, outputStr, stateStr, numGateInputs, numGateOutputs, numStates);
			}
		}
		else if (5 == sscanf(line.c_str(), "%s in %d <%[^>]> out %d <%[^>]>", type,
						&numGateInputs, inputStr, &numGateOutputs, outputStr)) {
			/*
			std::cout<<"################"<<"\ntype : "<<type<<"\n#in : "<<numGateInputs<<"\ninStr : "<<inputStr
			<<"\n#Out : "<<numGateOutputs<<"\noutStr : "<<outputStr
			<<"\n###################"<<std::endl;
			*/
			if (strcmp(type, "add") == 0) {
				assert(numGateOutputs == 1);
				handleAddition(inputStr, outputStr);
			} else if (strcmp(type, "mul") == 0) {
				assert(numGateInputs == 2 && numGateOutputs == 1);
				addMulConstraint(inputStr, outputStr);
			} else if (strcmp(type, "xor") == 0) {
				assert(numGateInputs == 2 && numGateOutputs == 1);
				addXorConstraint(inputStr, outputStr);
			} else if (strcmp(type, "or") == 0) {
				assert(numGateInputs == 2 && numGateOutputs == 1);
				addOrConstraint(inputStr, outputStr);
			} else if (strcmp(type, "assert") == 0) {
				assert(numGateInputs == 2 && numGateOutputs == 1);
				addAssertionConstraint(inputStr, outputStr);
			} else if (strstr(type, "const-mul-neg-")) {
				assert(numGateInputs == 1 && numGateOutputs == 1);
				handleMulNegConst(type, inputStr, outputStr);
			} else if (strstr(type, "const-mul-")) {
				assert(numGateInputs == 1 && numGateOutputs == 1);
				handleMulConst(type, inputStr, outputStr);
			} else if (strcmp(type, "zerop") == 0) {
				assert(numGateInputs == 1 && numGateOutputs == 2);
				addNonzeroCheckConstraint(inputStr, outputStr);
			} else if (strstr(type, "split")) {
				assert(numGateInputs == 1);
				addSplitConstraint(inputStr, outputStr, numGateOutputs);
			} else if (strstr(type, "pack")) {
				assert(numGateOutputs == 1);
				addPackConstraint(inputStr, outputStr, numGateInputs);
			}
			else if(strstr(type, "convol1d")){
				addConvol1DConstraint(inputStr, outputStr, numGateInputs, numGateOutputs);
			}
		} else {
//			assert(0);
		}

		delete[] inputStr;
		delete[] outputStr;
		delete[] stateStr;
		clean();
	}

	ifs2.close();

	printf("\tConstraint translation done\n");
	look_up_our_self(&usage2);
	unsigned long diff = usage2.vsize - usage1.vsize;
	printf("\tMemory usage for constraint translation: %lu MB\n", diff >> 20);

}

void CircuitReader::mapValuesToProtoboard() {

	int zeropGateIndex = 0;
	for (WireMap::iterator iter = variableMap.begin();
			iter != variableMap.end(); ++iter) {
		Wire wireId = iter->first;
		pb->val(*variables[variableMap[wireId]]) = wireValues[wireId];
		// std::cout<<"pb : "<<variableMap[wireId] <<" wire : "<<wireId<<std::endl;
		// std::cout<<"wire val : ";
		// wireValues[wireId].as_bigint().print();
		if (zeropMap.find(wireId) != zeropMap.end()) {
			LinearCombination l = *zeroPwires[zeropGateIndex++];
			if (pb->val(l) == 0) {
				pb->val(*variables[zeropMap[wireId]]) = 0;
			} else {
				pb->val(*variables[zeropMap[wireId]]) = pb->val(l).inverse(
						pb->fieldType_);
			}
		}
		//std::cout<<"pb : "<<variableMap[wireId] <<" wire : "<<wireId<<"("<<pb->val(*variables[variableMap[wireId]]).asLong()<<")\n";

	}
	if (!pb->isSatisfied(PrintOptions::DBG_PRINT_IF_NOT_SATISFIED)) {
		printf("Note: Protoboard Not Satisfied .. \n");
		// assert(false);
	}
	printf("Assignment of values done .. \n");

}

int CircuitReader::find(const Wire& wireId, LinearCombinationPtr& lc,
		bool intentionToEdit) {

	LinearCombinationPtr p = wireLinearCombinations[wireId];
	if (p) {
		//std::cout<<"p exist\n";
		wireUseCounters[wireId]--;
		if (wireUseCounters[wireId] == 0) {
			toClean.push_back(wireId);
			lc = p;
		} else {
			if (intentionToEdit) {
				lc = make_shared<LinearCombination>(*p);
			} else {
				lc = p;
			}
		}
		return 1;
	} else {
		//std::cout<<"p not exist\n";
		//std::cout<<"wire : "<<wireId<<", variables : "<<variableMap[wireId]<<"\n";
		wireUseCounters[wireId]--;
		lc = make_shared<LinearCombination>(
				LinearCombination(*variables[variableMap[wireId]]));
		//std::cout<<(*lc).asString()<<"\n";
		if (wireUseCounters[wireId] == 0) {
			toClean.push_back(wireId);
		}
		return 2;
	}
}

void CircuitReader::clean() {

	for (Wire wireId : toClean) {
		wireLinearCombinations[wireId].reset();
	}
	toClean.clear();
}

void CircuitReader::addMulConstraint(char* inputStr, char* outputStr) {

	Wire outputWireId, inWireId1, inWireId2;
	//std::cout<<inputStr<<std::endl;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;
	//std::cout<<inWireId1<<std::endl;
	iss_i >> inWireId2;
	//std::cout<<inWireId2<<std::endl;

	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;
	//std::cout<<outputWireId<<std::endl;



	LinearCombinationPtr l1, l2;
	find(inWireId1, l1);
	find(inWireId2, l2);

	//std::cout<<"l1 : "<<(*l1).asString()<<std::endl;

	if (variableMap.find(outputWireId) == variableMap.end()) {
		variables.push_back(make_shared<Variable>("mul out"));
		variableMap[outputWireId] = currentVariableIdx;
		//std::cout<<currentVariableIdx<<std::endl;

		pb->addRank1Constraint(*l1, *l2, *variables[currentVariableIdx],
				"Mul ..");
		currentVariableIdx++;
	} else {
		pb->addRank1Constraint(*l1, *l2, *variables[variableMap[outputWireId]],
				"Mul ..");
	}
	//std::cout<<"ok mul\n";
}

void CircuitReader::addXorConstraint(char* inputStr, char* outputStr) {

	Wire outputWireId, inWireId1, inWireId2;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;
	iss_i >> inWireId2;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;

	LinearCombinationPtr lp1, lp2;
	find(inWireId1, lp1);
	find(inWireId2, lp2);
	LinearCombination l1, l2;
	l1 = *lp1;
	l2 = *lp2;
	if (variableMap.find(outputWireId) == variableMap.end()) {
		variables.push_back(make_shared<Variable>("xor out"));
		variableMap[outputWireId] = currentVariableIdx;
		pb->addRank1Constraint(2 * l1, l2,
				l1 + l2 - *variables[currentVariableIdx], "XOR ..");
		currentVariableIdx++;
	} else {
		pb->addRank1Constraint(2 * l1, l2,
				l1 + l2 - *variables[variableMap[outputWireId]], "XOR ..");
	}
}

void CircuitReader::addOrConstraint(char* inputStr, char* outputStr) {

	Wire outputWireId, inWireId1, inWireId2;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;
	iss_i >> inWireId2;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;

	LinearCombinationPtr lp1, lp2;
	find(inWireId1, lp1);
	find(inWireId2, lp2);
	LinearCombination l1, l2;
	l1 = *lp1;
	l2 = *lp2;
	if (variableMap.find(outputWireId) == variableMap.end()) {
		variables.push_back(make_shared<Variable>("or out"));
		variableMap[outputWireId] = currentVariableIdx;
		pb->addRank1Constraint(l1, l2, l1 + l2 - *variables[currentVariableIdx],
				"OR ..");
		currentVariableIdx++;
	} else {
		pb->addRank1Constraint(l1, l2,
				l1 + l2 - *variables[variableMap[outputWireId]], "OR ..");
	}
}

void CircuitReader::addAssertionConstraint(char* inputStr, char* outputStr) {

	Wire outputWireId, inWireId1, inWireId2;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;
	iss_i >> inWireId2;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;

	LinearCombinationPtr lp1, lp2, lp3;
	find(inWireId1, lp1);
	find(inWireId2, lp2);
	find(outputWireId, lp3);

	LinearCombination l1, l2, l3;
	l1 = *lp1;
	l2 = *lp2;
	l3 = *lp3;
	pb->addRank1Constraint(l1, l2, l3, "Assertion ..");

}

void CircuitReader::addSplitConstraint(char* inputStr, char* outputStr,
		unsigned short n) {

	Wire inWireId;
	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId;

	LinearCombinationPtr l;
	find(inWireId, l);

	istringstream iss_o(outputStr, istringstream::in);

	LinearCombination sum;
	FElem two_i = libff::Fr<libff::default_ec_pp> ("1");

	/*
	for (int i = 0; i < n; i++) {
		Wire bitWireId;
		iss_o >> bitWireId;
		variables.push_back(make_shared<Variable>("bit out"));
		variableMap[bitWireId] = currentVariableIdx;
		VariablePtr vptr = variables[currentVariableIdx];
		pb->enforceBooleanity(*vptr);
		sum += LinearTerm(*vptr, two_i);
		two_i += two_i;
		currentVariableIdx++;
	} */

	for (int i = 0; i < n; i++) {
		Wire bitWireId;
		iss_o >> bitWireId;
		VariablePtr vptr;
		if (variableMap.find(bitWireId) == variableMap.end()) {
			variables.push_back(make_shared<Variable>("bit out"));
			variableMap[bitWireId] = currentVariableIdx;
			vptr = variables[currentVariableIdx];
			currentVariableIdx++;
		} else {
			vptr = variables[variableMap[bitWireId]];
		}
		pb->enforceBooleanity(*vptr);
		sum += LinearTerm(*vptr, two_i);
		two_i += two_i;
	}


	pb->addRank1Constraint(*l, 1, sum, "Split Constraint");
}

void CircuitReader::addPackConstraint(char* inputStr, char* outputStr,
		unsigned short n) {

	Wire outputWireId;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;

	istringstream iss_i(inputStr, istringstream::in);
	LinearCombination sum;
	FElem two_i = libff::Fr<libff::default_ec_pp> ("1");
	for (int i = 0; i < n; i++) {
		Wire bitWireId;
		iss_i >> bitWireId;
		LinearCombinationPtr l;
		find(bitWireId, l);
		sum += two_i * (*l);
		two_i += two_i;
	}

	VariablePtr vptr;
	if (variableMap.find(outputWireId) == variableMap.end()) {
		variables.push_back(make_shared<Variable>("pack out"));
		variableMap[outputWireId] = currentVariableIdx;
		vptr = variables[currentVariableIdx];
		currentVariableIdx++;
	} else {

		vptr = variables[variableMap[outputWireId]];
	}

	pb->addRank1Constraint(*vptr, 1, sum, "Pack Constraint");

}

void CircuitReader::addNonzeroCheckConstraint(char* inputStr, char* outputStr) {

	Variable auxConditionInverse_;
	Wire outputWireId, inWireId;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;
	iss_o >> outputWireId;
	LinearCombinationPtr l;

	find(inWireId, l);
	VariablePtr vptr;
	if (variableMap.find(outputWireId) == variableMap.end()) {
		variables.push_back(make_shared<Variable>("zerop out"));
		variableMap[outputWireId] = currentVariableIdx;
		vptr = variables[currentVariableIdx];
		currentVariableIdx++;
	} else {
		vptr = variables[variableMap[outputWireId]];
	}
	variables.push_back(make_shared<Variable>("zerop aux"));
	pb->addRank1Constraint(*l, 1 - *vptr, 0, "condition * not(output) = 0");
	pb->addRank1Constraint(*l, *variables[currentVariableIdx], *vptr,
			"condition * auxConditionInverse = output");

	zeroPwires.push_back(l);
	zeropMap[outputWireId] = currentVariableIdx;
	currentVariableIdx++;

}

void CircuitReader::handleAddition(char* inputStr, char* outputStr) {

	Wire inWireId, outputWireId;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;

	istringstream iss_i(inputStr, istringstream::in);
	LinearCombinationPtr s, l;
	iss_i >> inWireId;
	find(inWireId, l, true);
	s = l;
	while (iss_i >> inWireId) {
		find(inWireId, l);
		*s += *l;
	}
	wireLinearCombinations[outputWireId] = s;
}

void CircuitReader::handleMulConst(char* type, char* inputStr,
		char* outputStr) {

	char* constStr = type + sizeof("const-mul-") - 1;
	Wire outputWireId, inWireId;

	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;
	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId;
	LinearCombinationPtr l;
	find(inWireId, l, true);
	wireLinearCombinations[outputWireId] = l;
	//std::cout<<"id : "<<outputWireId<<" value : "<< wireValues[outputWireId].as_ulong()<<std::endl;
	*(wireLinearCombinations[outputWireId]) *= readFieldElementFromHex(
			constStr);
	//std::cout<<"const : "<<readFieldElementFromHex(
	//		constStr).as_ulong()<<std::endl;
}

void CircuitReader::handleMulNegConst(char* type, char* inputStr,
		char* outputStr) {

	char* constStr = type + sizeof("const-mul-neg-") - 1;
	Wire outputWireId, inWireId;
	istringstream iss_o(outputStr, istringstream::in);
	iss_o >> outputWireId;
	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId;

	LinearCombinationPtr l;
	find(inWireId, l, true);

	wireLinearCombinations[outputWireId] = l;
	*(wireLinearCombinations[outputWireId]) *= readFieldElementFromHex(
			constStr);
	*(wireLinearCombinations[outputWireId]) *= FieldT(-1); //TODO: make shared FieldT constants

}


/*add function*/
void CircuitReader::addConvol1DConstraint(char* inputStr, char* outputStr, unsigned short num_in, unsigned short num_out) {

	Wire outputWireId, inWireId1, inWireId2;
	//std::cout<<inputStr<<std::endl;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;

	
	//std::vector<Wire> inputIds;
	//std::cout<<wireValues[inWireId1].as_ulong()<<"\n";
	size_t num = wireValues[inWireId1].as_ulong();
	//inputIds.reserve(num);

	LinearCombinationPtr input_lc, kernel_lc, temp_lc;
	Wire wireTemp;
	iss_i >> wireTemp;
	find(wireTemp, temp_lc, true);
	input_lc = temp_lc;
	for(size_t i=0; i<num-1;i++){
		iss_i >> wireTemp;
		//std::cout<<"input id : "<<wireTemp<<std::endl;

		find(wireTemp, temp_lc);
		*input_lc += (*temp_lc)*(i+2);
		//inputIds.push_back(wireTemp);
	}

	//std::cout<<"id1 : "<<inWireId1<<"\ninput size : "<<num<<std::endl;
	iss_i >> inWireId2;

	//std::vector<Wire> kernelIds; 
	num = wireValues[inWireId2].as_ulong();
	//kernelIds.reserve(num);
	
	LinearCombinationPtr temp_lc2;
	iss_i >> wireTemp;
	find(wireTemp, temp_lc2, true);
	kernel_lc = temp_lc2;
	for(size_t i=0; i<num-1;i++){
		iss_i >> wireTemp;
		//std::cout<<"kernel id : "<<wireTemp<<std::endl;

		find(wireTemp, temp_lc2);
		*kernel_lc += (*temp_lc2)*(i+2);
		//kernelIds.push_back(wireTemp);
	}
	//std::cout<<"id2 : "<<inWireId2<<"\nkernel size : "<<num<<std::endl;;
	
	istringstream iss_o(outputStr, istringstream::in);
	//iss_o >> outputWireId;
	//std::cout<<outputWireId<<std::endl;
	//num = wireValues[outputWireId].as_ulong();

	LinearCombinationPtr output_lc;//, temp_lc3;
	//VariablePtr output_v;
	iss_o >> wireTemp;
	//std::cout<<wireTemp<<std::endl;

	//std::cout<<"curIdx : "<<currentVariableIdx<<"\n";
	//find(wireTemp, temp_lc3, true);
	//output_lc = temp_lc3;
	if (variableMap.find(wireTemp) == variableMap.end()) {
		//std::cout<<"1"<<std::endl;

		variables.push_back(make_shared<Variable>("convol out"+std::to_string(0)));
		//std::cout<<"11"<<std::endl;

		variableMap[wireTemp] = currentVariableIdx;
		//std::cout<<"12"<<std::endl;
		//std::cout<<currentVariableIdx<<std::endl;
		//output_v = variables[currentVariableIdx];
		output_lc = make_shared<LinearCombination>(LinearCombination(*variables[currentVariableIdx]));
		//std::cout<<"out = cv["<<currentVariableIdx<<"]";
		currentVariableIdx++;
	} else {
		//std::cout<<"2"<<std::endl;

		//output_v = variables[variableMap[outputWireId]];
		output_lc = make_shared<LinearCombination>(LinearCombination(*variables[variableMap[wireTemp]]));
		//std::cout<<"out = v["<<variableMap[wireTemp]<<"]";
		//output_lc = *variables[variableMap[outputWireId]];
	}
	//std::cout<<"first ok"<<std::endl;

	for(size_t i=0; i<num_out-1;i++){
		iss_o >> wireTemp;
		//find(wireTemp, temp_lc3);
		//*output_lc += *temp_lc3;
		//kernelIds.push_back(wireTemp);
		if (variableMap.find(wireTemp) == variableMap.end()) {
			variables.push_back(make_shared<Variable>("convol out"+std::to_string(i+1)));
			variableMap[wireTemp] = currentVariableIdx;
			//std::cout<<currentVariableIdx<<std::endl;
			//output_v = variables[currentVariableIdx];
			*output_lc += (*variables[currentVariableIdx])*(i+2);
			//cout<<"+cv["<<currentVariableIdx<<"]";
			currentVariableIdx++;
		} else {
			*output_lc += (*variables[variableMap[wireTemp]])*(i+2);
			//cout<<"+v["<<variableMap[wireTemp]<<"]";
		}
	}
	//std::cout<<wireTemp<<std::endl;

	//std::cout<<"output size : "<<num_out<<std::endl;;

	pb->convol_outputs_size.insert(pb->convol_outputs_size.end(),num_out);
	pb->convol_size++;
	LinearCombination dum1, dum2, dum3;
	pb->addRank1Constraint(dum1, dum2, dum3, *input_lc, *kernel_lc, *output_lc,"Convol ..");
	
}

void CircuitReader::addConvol2DConstraint(char* inputStr, char* outputStr, char* stateStr, unsigned short num_in, unsigned short num_out, unsigned short num_state) {

	//Wire inputHeightW, inputWidthW, kernelHeightW, kernelWidthW;
	//std::cout<<"ok2d\n";
	//std::cout<<stateStr<<std::endl;
	size_t inputHeight, inputWidth, kernelHeight, kernelWidth;
	istringstream iss_st(stateStr, istringstream::in);
	iss_st >> inputHeight;
	iss_st >> inputWidth;
	iss_st >> kernelHeight;
	iss_st >> kernelWidth;

	//std::cout<<"in height : "<<inputHeight<<" width : "<<inputWidth<<" kernel : "<<kernelHeight<<" width :"<<kernelWidth<<std::endl;

	istringstream iss_i(inputStr, istringstream::in);

	LinearCombinationPtr input_lcH, input_lcW, temp_lcH, temp_lcW;
	Wire wireTemp;
	iss_i >> wireTemp;
	find(wireTemp, temp_lcH, true);
	input_lcH = temp_lcH;
	input_lcW = make_shared<LinearCombination>(*temp_lcH);
	for(size_t i=0;i<inputHeight;i++){
		for(size_t j=0;j<inputWidth;j++){
			if(i==0 && j==0){
				;
			}else{
			iss_i >> wireTemp;
			//std::cout<<"input id : "<<wireTemp<<" : "<<variableMap[wireTemp]<<std::endl;

			find(wireTemp, temp_lcH);
			temp_lcW = make_shared<LinearCombination>(*temp_lcH);
			*input_lcH += (*temp_lcH)*(i+1);
			*input_lcW += (*temp_lcW)*(j+1);
			}
		}	
	}

	LinearCombinationPtr kernel_lcH, kernel_lcW;
	iss_i >> wireTemp;
	find(wireTemp, temp_lcH, true);
	kernel_lcH = temp_lcH;
	kernel_lcW = make_shared<LinearCombination>(*temp_lcH);
	for(size_t i=0;i<kernelHeight;i++){
		for(size_t j=0;j<kernelWidth;j++){
			if(i==0 && j==0) {;}
			else{
			iss_i >> wireTemp;
			//std::cout<<"kernel id : "<<wireTemp<<std::endl;

			find(wireTemp, temp_lcH);
			temp_lcW = make_shared<LinearCombination>(*temp_lcH);
			*kernel_lcH += (*temp_lcH)*(i+1);
			*kernel_lcW += (*temp_lcW)*(j+1);
			}
		}
	}

	istringstream iss_o(outputStr, istringstream::in);

	LinearCombinationPtr output_lcW, output_lcH;
	iss_o >> wireTemp;

	//std::cout<<"curIdx : "<<currentVariableIdx<<"\n";

	if (variableMap.find(wireTemp) == variableMap.end()) {

		variables.push_back(make_shared<Variable>("convol out2d"+std::to_string(0)));

		variableMap[wireTemp] = currentVariableIdx;
		output_lcH = make_shared<LinearCombination>(LinearCombination(*variables[currentVariableIdx]));
		output_lcW = make_shared<LinearCombination>(LinearCombination(*variables[currentVariableIdx]));
		//std::cout<<"out = cv["<<currentVariableIdx<<"]";
		currentVariableIdx++;
	} else {
		//std::cout<<"2"<<std::endl;

		//output_v = variables[variableMap[outputWireId]];
		output_lcW = make_shared<LinearCombination>(LinearCombination(*variables[variableMap[wireTemp]]));
		output_lcH = make_shared<LinearCombination>(LinearCombination(*variables[variableMap[wireTemp]]));

		//std::cout<<"out = v["<<variableMap[wireTemp]<<"]";
		//output_lc = *variables[variableMap[outputWireId]];
	}

	size_t outputHeight = inputHeight + kernelHeight - 1;
	size_t outputWidth = inputWidth + kernelWidth - 1;
	for(size_t i=0;i<outputHeight;i++){
		for(size_t j=0;j<outputWidth;j++){
			if(i == 0 && j==0){;}
			else{
				iss_o >> wireTemp;
				if (variableMap.find(wireTemp) == variableMap.end()) {
					variables.push_back(make_shared<Variable>("convol2d out"+std::to_string(i+1)));
					variableMap[wireTemp] = currentVariableIdx;
					*output_lcH += (*variables[currentVariableIdx])*(i+1);
					*output_lcW += (*(make_shared<Variable>(*variables[currentVariableIdx])))*(j+1);
					//cout<<"+cv["<<currentVariableIdx<<"]";
					currentVariableIdx++;
				} else {
					*output_lcH += (*variables[variableMap[wireTemp]])*(i+1);
					*output_lcW += (*(make_shared<Variable>(*variables[variableMap[wireTemp]])))*(j+1);

					//cout<<"+v["<<variableMap[wireTemp]<<"]";
				}
			}
			
		}
	}

	//std::cout<<"output size : "<<outputHeight*outputWidth<<std::endl;;

	pb->convol_outputs_size.insert(pb->convol_outputs_size.end(), inputHeight + kernelHeight-1);
	pb->convol_outputs_size2.insert(pb->convol_outputs_size2.end(), inputWidth + kernelWidth -1);

	pb->convol_input_height.insert(pb->convol_input_height.end(),inputHeight);
	pb->convol_input_width.insert(pb->convol_input_width.end(),inputWidth);
	pb->convol_kernel_height.insert(pb->convol_kernel_height.end(),kernelHeight);
	pb->convol_kernel_width.insert(pb->convol_kernel_width.end(),kernelWidth);

	pb->convol_dimensions.insert(pb->convol_dimensions.end(),2);
	pb->convol_size++;
	LinearCombination dum1, dum2, dum3;
	pb->addRank1Constraint(dum1, dum2, dum3, *input_lcH, *kernel_lcH, *output_lcH, *input_lcW, *kernel_lcW, *output_lcW,"Convol2D ..");
	
}

void CircuitReader::addConvolConstraint(char* inputStr, char* outputStr, unsigned short num_in, unsigned short num_out) {

	Wire outputWireId, inWireId1, inWireId2;
	//std::cout<<inputStr<<std::endl;

	istringstream iss_i(inputStr, istringstream::in);
	iss_i >> inWireId1;

	
	//std::vector<Wire> inputIds;
	//std::cout<<wireValues[inWireId1].as_ulong()<<"\n";
	size_t num = wireValues[inWireId1].as_ulong();
	//inputIds.reserve(num);

	LinearCombinationPtr input_lc, kernel_lc, temp_lc;
	Wire wireTemp;
	iss_i >> wireTemp;
	find(wireTemp, temp_lc, true);
	input_lc = temp_lc;
	for(size_t i=0; i<num-1;i++){
		iss_i >> wireTemp;
		//std::cout<<"input id : "<<wireTemp<<std::endl;

		find(wireTemp, temp_lc);
		*input_lc += (*temp_lc)*(i+2);
		//inputIds.push_back(wireTemp);
	}

	//std::cout<<"id1 : "<<inWireId1<<"\ninput size : "<<num<<std::endl;
	iss_i >> inWireId2;

	//std::vector<Wire> kernelIds; 
	num = wireValues[inWireId2].as_ulong();
	//kernelIds.reserve(num);
	
	LinearCombinationPtr temp_lc2;
	iss_i >> wireTemp;
	find(wireTemp, temp_lc2, true);
	kernel_lc = temp_lc2;
	for(size_t i=0; i<num-1;i++){
		iss_i >> wireTemp;
		//std::cout<<"kernel id : "<<wireTemp<<std::endl;

		find(wireTemp, temp_lc2);
		*kernel_lc += (*temp_lc2)*(i+2);
		//kernelIds.push_back(wireTemp);
	}
	//std::cout<<"id2 : "<<inWireId2<<"\nkernel size : "<<num<<std::endl;;
	
	istringstream iss_o(outputStr, istringstream::in);
	//iss_o >> outputWireId;
	//std::cout<<outputWireId<<std::endl;
	//num = wireValues[outputWireId].as_ulong();

	LinearCombinationPtr output_lc;//, temp_lc3;
	//VariablePtr output_v;
	iss_o >> wireTemp;
	//std::cout<<wireTemp<<std::endl;

	//find(wireTemp, temp_lc3, true);
	//output_lc = temp_lc3;
	if (variableMap.find(wireTemp) == variableMap.end()) {
		//std::cout<<"1"<<std::endl;

		variables.push_back(make_shared<Variable>("convol out"+std::to_string(0)));
		//std::cout<<"11"<<std::endl;

		variableMap[wireTemp] = currentVariableIdx;
		//std::cout<<"12"<<std::endl;
		//std::cout<<currentVariableIdx<<std::endl;
		//output_v = variables[currentVariableIdx];
		output_lc = make_shared<LinearCombination>(LinearCombination(*variables[currentVariableIdx]));
		//std::cout<<"out = cv["<<currentVariableIdx<<"]";
		currentVariableIdx++;
	} else {
		//std::cout<<"2"<<std::endl;

		//output_v = variables[variableMap[outputWireId]];
		output_lc = make_shared<LinearCombination>(LinearCombination(*variables[variableMap[wireTemp]]));
		//std::cout<<"out = v["<<variableMap[wireTemp]<<"]";
		//output_lc = *variables[variableMap[outputWireId]];
	}
	//std::cout<<"first ok"<<std::endl;

	for(size_t i=0; i<num_out-1;i++){
		iss_o >> wireTemp;
		//find(wireTemp, temp_lc3);
		//*output_lc += *temp_lc3;
		//kernelIds.push_back(wireTemp);
		if (variableMap.find(wireTemp) == variableMap.end()) {
			variables.push_back(make_shared<Variable>("convol out"+std::to_string(i+1)));
			variableMap[wireTemp] = currentVariableIdx;
			//std::cout<<currentVariableIdx<<std::endl;
			//output_v = variables[currentVariableIdx];
			*output_lc += (*variables[currentVariableIdx])*(i+2);
			//cout<<"+cv["<<currentVariableIdx<<"]";
			currentVariableIdx++;
		} else { 
			*output_lc += (*variables[variableMap[wireTemp]])*(i+2);
			//cout<<"+v["<<variableMap[wireTemp]<<"]";
		}
	}
	std::cout<<wireTemp<<std::endl;

	//std::cout<<"output size : "<<num_out<<std::endl;;

	pb->convol_outputs_size.insert(pb->convol_outputs_size.end(),num_out);
	pb->convol_size++;
	LinearCombination dum1, dum2, dum3;
	pb->addRank1Constraint(dum1, dum2, dum3, *input_lc, *kernel_lc, *output_lc,"Convol ..");
	
}
