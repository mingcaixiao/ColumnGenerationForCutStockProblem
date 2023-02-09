#include <vector>
#include <iostream>
#include <fstream>

#include "absl/flags/flag.h"
#include "ortools/base/flags.h"
#include "ortools/base/init_google.h"
#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/model_exporter.h"

//#define WriteFile

using pattern = std::map<int,int>;
struct Data
{
	std::vector<pattern> pats;
	double rollLen;
	std::vector<double> len;
	std::vector<double> demand;
};

void writeStringToFile(const std::string& string,const std::string& filename)
{
	std::ofstream os(filename);
	if (os.is_open()){
		os << string;
		os.close();
	}else{
		std::cout << "filename is not open for write!" << std::endl;
	}
}
void cutstock(double rollLen,std::vector<double>& len,std::vector<double>&demand)
{
	using namespace operations_research;
	Data data;
	data.rollLen = rollLen;
	data.demand = demand;
	data.len = len;
	for (int i = 0; i < len.size(); i++){
		pattern pat;
		pat.insert({i,int(rollLen/len.at(i))});
		data.pats.push_back(pat);
	}
	auto solveMasterProblem = [](Data& data)->std::vector<double> {
		std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
		if (!solver) {
			LOG(WARNING) << "CLP solver unavailable.";
		}
		const double infinity = solver->infinity();
		//create variable
		std::vector<MPVariable*> vars(data.pats.size());
		for (int j=0;j<data.pats.size();j++){
			MPVariable* const x = solver->MakeNumVar(0.0,infinity,"");
			vars.at(j)=x;
		}
		//create cons
		std::vector<MPConstraint*> cons(data.demand.size());
		int rowIndex = 0;
		for (int i=0;i<data.demand.size();i++){
			MPConstraint* const c = solver->MakeRowConstraint(data.demand.at(i), infinity);
			int colIndex = 0;
			for (auto& pat : data.pats){
				if(pat.find(rowIndex)!=pat.end())
					c->SetCoefficient(vars.at(colIndex),pat.at(rowIndex));
				colIndex++;
			}
			cons.at(i)=c;
			rowIndex++;
		}
		//create objective
		MPObjective* const obj = solver->MutableObjective();
		for (int j=0;j<data.pats.size();j++){
			obj->SetCoefficient(vars.at(j),1);
		}
		obj->SetMinimization();
#ifdef WriteFile
		std::string lpString;
		if (solver->ExportModelAsLpFormat(true, &lpString)){
			writeStringToFile(lpString,"master_prob.lp");
			std::cout << "export lp file succeed" << std::endl;
		}else{
			std::cout << "export lp file failed" << std::endl;
		}
#endif
		MPSolverParameters param;
		param.SetIntegerParam(MPSolverParameters::LP_ALGORITHM, MPSolverParameters::LpAlgorithmValues::DUAL);
		MPSolver::ResultStatus resultStatus=solver->Solve(param);
		if (resultStatus != MPSolver::ResultStatus::OPTIMAL) {
			LOG(INFO) << "master problem don't get optimal!";
			std::string lpString;
			solver->ExportModelAsMpsFormat(true, false,&lpString);
			writeStringToFile(lpString,"wrong_master_prob.mps");
		}
		std::vector<double> dualPrices(cons.size(),0.0);
		for (int i=0;i<cons.size();i++){
			dualPrices.at(i) = cons.at(i)->dual_value();
		}
		std::cout << "linear releax value is "<< obj->Value()<< std::endl;
		return dualPrices;
	};
	auto solveSubProblem = [](std::vector<double>& dualPrices,Data& data,pattern& pat)->double {
		std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
		if (!solver){
			LOG(WARNING) << "Solver SCIP unavailable!";
		}
		//create variables
		std::vector<MPVariable*> vars(data.demand.size());
		for (int i = 0; i < data.demand.size(); i++){
			MPVariable* const x = solver->MakeIntVar(0.0,int(data.rollLen/data.len.at(i)),"");
			vars.at(i) = x;
		}
		//create constraint
		MPConstraint* const cons = solver->MakeRowConstraint(0.0,data.rollLen);
		for (int i = 0; i < vars.size(); i++){
			cons->SetCoefficient(vars.at(i),data.len.at(i));
		}

		//create obj
		MPObjective* const obj = solver->MutableObjective();
		for (int i = 0; i < vars.size(); i++){
			obj->SetCoefficient(vars.at(i), dualPrices.at(i));
		}
		obj->SetMaximization();
#ifdef WriteFile
		std::string lpString;
		if (solver->ExportModelAsLpFormat(true, &lpString)){
			writeStringToFile(lpString, "sub_prob.lp");
			std::cout << "export lp file succeed" << std::endl;
		}else{
			std::cout << "export lp file failed" << std::endl;
		}
#endif
		solver->Solve();
		
		double reduceCost = 1;
		pat.clear();
		for (int i = 0; i < vars.size(); i++){
			double value=vars.at(i)->solution_value();
			if (value != 0)
			{
				pat.insert({i,value});
			}
			reduceCost -= dualPrices.at(i)*value;
		}
		return reduceCost;
	};
	auto solveFinalProblem = [](Data& data) {
		std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
		if (!solver) {
			LOG(WARNING) << "GLOP solver unavailable.";
		}
		const double infinity = solver->infinity();
		//create variable
		std::vector<MPVariable*> vars(data.pats.size());
		for (int i=0;i<data.pats.size();i++){
			MPVariable* const x = solver->MakeIntVar(0.0, infinity, "");
			vars.at(i)=x;
		}
		//create cons
		std::vector<MPConstraint*> cons(data.demand.size());
		int rowIndex = 0;
		for (int i=0;i<data.demand.size();i++){
			MPConstraint* const c = solver->MakeRowConstraint(data.demand.at(i), infinity);
			int colIndex = 0;
			for (auto& pat : data.pats){
				if (pat.find(rowIndex) != pat.end())
					c->SetCoefficient(vars.at(colIndex), pat.at(rowIndex));
				colIndex++;
			}
			cons.at(i)=c;
			rowIndex++;
		}
		//create objective
		MPObjective* const obj = solver->MutableObjective();
		for (int j = 0; j < data.pats.size(); j++){
			obj->SetCoefficient(vars.at(j), 1);
		}
		obj->SetMinimization();
		solver->set_time_limit(1000*60);
		MPSolver::ResultStatus resultStatus = solver->Solve();
		double objValue = obj->Value();
		std::cout << "use roll num is " << objValue<<std::endl;
		//print solution
		std::cout << "solution is :" << std::endl;
		int validIndex = 0;
		int i = 0;
		std::cout << "pattern" << "  num" << std::endl;
		for (auto& var:vars){
			if (var->solution_value() != 0){
				//std::cout << "pattern_" << validIndex++ << ":";
				for (auto pat : data.pats.at(i))
				{
					std::cout << data.len.at(pat.first)<<"*"<<pat.second<<" ";
				}
				std::cout << var->solution_value();
				std::cout << std::endl;
			}
			i++;
		}
	};
	while (true){
		auto dualPrices=solveMasterProblem(data);
		pattern newPattern;
		double reduceCost=solveSubProblem(dualPrices,data,newPattern);
		std::cout << "reduce cost is " << reduceCost << std::endl;
		if (reduceCost >= -1e-3){
			break;
		}else{
			data.pats.push_back(newPattern);
		}
	}
	solveFinalProblem(data);
}

std::ostream& operator <<(std::ostream& os, std::vector<double>& a)
{
	os << "[";
	for (auto& i : a) {
		os << i<<",";
	}
	os << "]";
	return os;
}

std::istream& operator>>(std::istream& in, std::vector<double>& a)
{
	a.clear();
	char ch;
	in >> ch;
	while (true){
		double x;
		in >> x >> ch;
		a.push_back(x);
		if (ch == ']')break;
	}
	return in;
}
void readData(char* filename,double& rollLen,std::vector<double>&len,std::vector<double>&demand)
{
	std::ifstream in(filename);
	if (in.is_open()){
		in >> rollLen >> len >> demand;
		in.close();
	}else{
		std::cout << "file format is incorrect!" << std::endl;
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[])
{
	//std::cout << argv[1] << std::endl;
	if (argc != 2){
		std::cout << "usage:cutstock.exe inputFile" << std::endl;
		exit(EXIT_FAILURE);
	}else{
		double rollLen;
		std::vector<double>len;
		std::vector<double>demand;
		readData(argv[1],rollLen,len,demand);
		
		std::cout << "stockLen is " << rollLen << "\n"
			<<"len is" << len << "\n" 
			<< "demand is"<<demand << std::endl;

		cutstock(rollLen, len, demand);
	}
	return EXIT_SUCCESS;
}