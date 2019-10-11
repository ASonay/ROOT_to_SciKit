#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TTreeFormula.h"

using namespace std;

class ReadTree
{

public:
  ReadTree(string filoc, string finame, string tname, vector<string> var);
  ReadTree(string reader_name, TTree *tr,  vector<string> var);
  ~ReadTree();

  vector<double> GetInput(int n);
  double GetInputSingle(int n);
  void SetSingleVariable(string variable);
  void Close();

private:

  TFile *fi;
  TTree *tree_raw;
  TTree *tree;
  
  vector<string> variables;
  vector<TTreeFormula*> tf;
  TTreeFormula *tf_s;
  TTreeFormula *formx;
  
  string name;

  vector<double> xv;
  double xd;
  
};


ReadTree::ReadTree(string filoc, string finame, string tname, vector<string> var):
  tree(NULL),
  tf_s(NULL),
  formx(NULL)
{
  variables = var;
  name = finame;

  fi = new TFile( (filoc+finame+".root").c_str() );
  
  tree = (TTree*)fi->Get(tname.c_str());

  tf.clear();
  for (int i=0;i<variables.size();i++){
    formx = new TTreeFormula((name+"_"+to_string(i)).c_str(),variables[i].c_str(),tree);
    tf.push_back(formx);
  }

}

ReadTree::ReadTree(string reader_name, TTree *tr, vector<string> var):
  tree(NULL),
  tf_s(NULL),
  formx(NULL)
{
  variables = var;
  name = reader_name;
  
  tree = tr;

  tf.clear();
  for (int i=0;i<variables.size();i++){
    formx = new TTreeFormula((name+"_"+to_string(i)).c_str(),variables[i].c_str(),tree);
    tf.push_back(formx);
  }
 
}

ReadTree::~ReadTree(){
}

vector<double> ReadTree::GetInput(int n=0){

  xv.clear();
  tree->GetEntry(n);
  for (int i=0;i<variables.size();i++)
    xv.push_back(tf[i]->EvalInstance());

  return xv;
}

void ReadTree::SetSingleVariable(string variable){
  delete tf_s;
  tf_s = new TTreeFormula(name.c_str(),variable.c_str(),tree);
}

double ReadTree::GetInputSingle(int n=0){
  tree->GetEntry(n);
  xd = tf_s->EvalInstance();
  return xd;
}


void ReadTree::Close(){
  delete tree;
  delete formx;
  delete tf_s;
}
