#include "TCut.h"
#include "inc/ReadTree.hh"

using namespace std;

struct TrainingVars
{
  vector<string> variables;
  string weight_signal;
  string weight_background;
};

TrainingVars ReadVariables(string variable_list){
  ifstream in(variable_list);
  string str;
  TrainingVars GetVars;

  if (in.fail()) {cout << "Your variable list is missing!\n"; exit(0);}
  
  while (in >> str){
    if (str.compare("var:")==0){
      in >> str;
      GetVars.variables.push_back(str);
    }
    else if (str.compare("sgw:")==0){
      in >> str;
      GetVars.weight_signal=str;
    }
    else if (str.compare("bgw:")==0){
      in >> str;
      GetVars.weight_background=str;
    }
    else
      cout << "Your file is badly formatted.." << endl;
  }

  return GetVars;
}

void ntuple_to_csv(string ntuple = "/home/asonay/ATLAS/MVA_result/SM4t-212560_ntuples/offline/1LOS/mvatuple_7ji2bi/4top_lo_simon.root",
		   string tree_name = "nominal_Loose",
		   string label = "s",
		   string variable_list = "list/VariableList.txt",
		   string scut = "nBTags_MV2c10_70>=3&&nJets>=10",
		   string output_file_name = "data/training_signal_1l_10ji3bi.csv"
		   )
{
  

  TrainingVars GetVars = ReadVariables(variable_list);
  vector<string> variables = GetVars.variables;
  string weight_signal = GetVars.weight_signal;
  string weight_background = GetVars.weight_background;

  TCut cuts = scut.c_str();

  TFile *f = new TFile(ntuple.c_str());

  TTree *tr = (TTree*)f->Get(tree_name.c_str());

  TTree *tree = NULL;
  gROOT->cd();
  cout << "Copying trees patiance.." << endl;
  if (scut.size()>0)
    tree = (TTree*)tr->CopyTree(cuts);
  else
    tree = tr;

  
  unsigned n = tree->GetEntries();


  ReadTree read("Data",tree,variables);

  read.SetSingleVariable(weight_signal);

  vector<double> var;
  double w;

  ofstream myfile;
  myfile.open (output_file_name);

  myfile << "id" << ",";
  for (int i=0;i<variables.size();i++)
    myfile << variables[i] << ",";
   
  myfile << "weights" << ",";
  myfile << "label" << endl;

    
  int cous=0,coub=0;
  for (int i=0;i<n;i++){
      var = read.GetInput(i);
      w = read.GetInputSingle(i);

      myfile << i << ",";
      for (int j=0;j<variables.size();j++)
	myfile << var[j] << ",";
      myfile << w << ",";
      myfile << label << endl;
 
  }
  
  myfile.close();


}
