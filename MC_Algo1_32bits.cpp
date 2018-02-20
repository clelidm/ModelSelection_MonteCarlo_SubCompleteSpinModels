//g++ -O3 main.cpp
//time ./a.out
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <map>
#include <random>  // for Mersen Twister
#include <set> 

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
//info datafile: TO SPECIFY:
    const string data_file_name = "dir1/dataset.txt"; // directory/name of the datafile
    const unsigned int n = 9;                    // number of spins
    const unsigned int m = 3;                    // chosen number of independent operators

//MC output file: 
    const string directory = "dir2/";  // directory for the outputfile
    const string MC_fileOUT = directory + "outputfile_name";  // name of the output file

//for Monte Carlo: 
    const int N_MCsample = 1e5;       // number of MC steps
    const double beta = 0.1;          // 1./temperature
    const int N_shuffling = 5;        // see doc.
    const double eps = 1e-2;    //for post-treatment: function "post_treatment_bases()"

//do not modify:
    const uint32_t un = 1;
    const uint32_t Nop = (un << m) - 1;          // number of operators = 2^m - 1
    unsigned int N = 0;                          // data set size

//Structure containing the final information on the best model, at the end of the MC
struct FinalModel
{
  set<uint32_t> basis;   // one basis of the model
  set<uint32_t> all_op;  // all operators of the model
  double LogL;           // Max-logL + N(n-m)log(2) 
  bool operator < (const FinalModel &other) const { return all_op < other.all_op; }
};
/********************************************************************/
/***********************   RANDOM GENERATORS    *********************/
/********************************************************************/
std::mt19937 gen_32;   // Mersenne Twister
uniform_int_distribution<uint32_t> rand_index(1, n);
uniform_int_distribution<uint32_t> rand_indexm(1, m);
// INFO: for a random double in (0, 1)  --> use drand48();

void init_rand()
{
  int seed_val = (unsigned)time(NULL);
  srand48(seed_val);
  gen_32.seed(seed_val);
}

void rand_2index(uint32_t *i, uint32_t *j)     //sample 2 different indices (i,j) in [1:n]
{
  *i = rand_index(gen_32);
  *j = rand_index(gen_32); 
  while( (*i) == (*j) ) { *j = rand_index(gen_32); }
}

void rand_2indexm(uint32_t *i, uint32_t *j)     //sample 2 different indices (i,j) in [1:n]
{
  *i = rand_indexm(gen_32);
  *j = rand_index(gen_32); 
  while( (*i) == (*j) ) { *j = rand_index(gen_32); }
}

/********************************************************************/
/****************************    TOOLS    ***************************/
/********************************************************************/
/*****************    COMPARATOR for SET of BASES    ****************/
struct Cmp
{
    bool operator()(const pair< double, set<uint32_t> > &a, const pair< double, set<uint32_t> > &b)
    {
      return (a.first > b.first || (a.first == b.first && a.second < b.second) );
    }
};

struct Cmp_models
{
    bool operator()(const pair< double, set<uint32_t> > &a, const pair< double, set<uint32_t> > &b)
    {
      return (a.first > b.first || (a.first == b.first && a.second < b.second) );
    }
};
/**************    FILENAMES : Precision of double   ****************/
template <typename T>
string to_string_with_precision(const T a_value, const int n = 4)
{
    ostringstream out;
    out << setprecision(n) << a_value;
    return out.str();
}
/**************    FILENAMES : Scientific writing   *****************/
string to_string_scientific(int num)
{
  ostringstream streamObj;        // Create an output string stream
  streamObj << setprecision(1);   // Set precision to 2 digits
  streamObj << (double) num;      //Add double to stream
  return streamObj.str();  // Get string from output string stream
}

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
map<uint32_t, unsigned int> read_datafile()   // O(N)  //N = data set size
{
  string line, line2;     uint32_t nb = 0;
  N = 0;            // N = dataset size (global variable)
  cout << endl << "--->> Read \"" << data_file_name << "\",\t Build Nset...";

// ***** data are store in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set
  
  ifstream myfile (data_file_name.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line
      nb = bitset<n>(line2).to_ulong();   //convert string line2 into a binary integer
      Nset[nb] += 1;
      //cout << line << endl;   //cout << nb << " :  " << bitset<n>(nb) << endl;
      N++;
    }
    myfile.close();
  }
  else cout << "Unable to open file"; 

  cout << "\t\t data size N = " << N << endl;

  return Nset;
}

/******************************************************************************/
/************************** N_SET *********************************************/
/******************************************************************************/
void read_Nset (map<uint32_t, unsigned int> Nset)
// map.second = nb of time that the state map.first appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  int Ncontrol = 0;

  cout << endl << "Nset:" << endl;
  cout << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    cout << it->first << ":\t" << bitset<n>(it->first) << " => " << it->second; // << endl;
    cout << "  \t  P = " << it->second / (float) N << endl;
    Ncontrol += it->second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }
}

void distrib_Nset(map<uint32_t, unsigned int> Nset)  // Build the HISTOGRAM of Nset
{
  map<unsigned int, unsigned int> histo_set;  // histo_set[k] = #of states that appears k times in the data set
  map<uint32_t, unsigned int>::iterator it;
  int Ncontrol = 0;

// ***** Build Histogram:  ********************************
  cout << endl << "--->> Build distribution of Nset..." << endl;
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    histo_set[(it->second)] += 1;
    //cout << it->first << ":  " << bitset<n>(it->first) << " => " << it->second << '\n';
  }

// ***** Print Histogram into a file:  *********************
  fstream fichier( (data_file_name + "_histo_Nset.dat").c_str() , ios::out); //nom_fich.c_str()
  map<unsigned int, unsigned int>::iterator histo_it;

  fichier << "# 1: nb of occurence of a state \t 2: nb of states that occur [$1] times " << endl;
  for(histo_it = histo_set.begin(); histo_it!=histo_set.end(); ++histo_it)
    {
      fichier << histo_it->first << " " << histo_it->second << endl;
      Ncontrol += (histo_it->first)*(histo_it->second);
    }
  fichier.close();

  if (Ncontrol != N) { cout << "Error function \'distrib_Nset\': Ncontrol != N" << endl;  }
}

/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/
map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > 
        build_Kset(map<uint32_t, unsigned int> Nset, bool print_bool)
// mu_m = mu cut on the n first spins 
// Kset[mu_m].first = #of time state mu_m appears in the data set
// Kset[mu_m].second = list of pair<mu, N_mu>, 
//                          where mu are the states such that mu|m = mu_m
//                          and N_mu is the # of times state mu appears in the data
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset;

  uint32_t mu;
  unsigned int N_mu=0;
  uint32_t mu_m;    // mu_m = mu cut to the m first spins

  cout << endl << "--->> Build Kset..." << endl;

//Build Kset:
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    mu = it->first;       // state mu
    N_mu = it->second;    // # of times mu appears in the data set
    mu_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {  cout << mu << ": " << bitset<n>(mu) << "\t" << mu_m << ": " << bitset<m>(mu_m) << endl; }

    Kset[mu_m].first += N_mu;
    Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  cout << endl;

  return Kset;
}

void read_Kset(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset)
{
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it;
  unsigned int Ncontrol = 0;
  uint32_t mu_control = 0;

  cout << "Kset:" << endl;
  cout << "#1: m-state \t #2: nb of pts in m-state \t #3: Pba m-state \t #4 Pba state" << endl;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ncontrol += (it->second).first;
    while(mu_control < (it->first))  
      { 
        cout << mu_control << ":\txxxx" << bitset<m>(mu_control) << " :  NULL"; // << endl;
        cout << "\t  P_{" << bitset<m>(mu_control) <<  "} = " << 0; //<< endl;
        cout << "  \t\t    P_{xxxx" << bitset<m>(mu_control) <<  "} = " << 0 << endl;        
        mu_control++;
      }
    cout << it->first << ":\txxxx" << bitset<m>(it->first) << " : " << (it->second).first << " pts"; //<< endl;
    cout << "\t  P_{" << bitset<m>(it->first) <<  "} = " << (it->second).first / (float) N; //<< endl;
    cout << "  \t    P_{xxxx" << bitset<m>(it->first) <<  "} = " << (it->second).first / (float) N / (2*2*2*2) << endl;
    mu_control++;
  }
  cout << endl;

  if (Ncontrol != N) { cout << "Error read_Kset function: Ncontrol != N" << endl;  }
}

void read_Kset_full(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset)
{
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it;
  list<pair<uint32_t, unsigned int> >::iterator it_list;
  unsigned int Ncontrol = 0;
  uint32_t mu_control = 0;

  cout << "read2 Kset:" << endl;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ncontrol += (it->second).first;
    while(mu_control < (it->first))  
      { 
        cout << mu_control << ":\txxxx" << bitset<m>(mu_control) << " :  NULL"; // << endl;
        cout << "\t  P_{" << bitset<m>(mu_control) <<  "} = " << 0; //<< endl;
        cout << "  \t\t    P_{xxxx" << bitset<m>(mu_control) <<  "} = " << 0 << endl;        
        mu_control++;
      }
    cout << it->first << ":\txxxx" << bitset<m>(it->first) << " : " << (it->second).first << " pts"; //<< endl;
    cout << "\t  P_{" << bitset<m>(it->first) <<  "} = " << (it->second).first / (float) N; //<< endl;
    cout << "  \t    P_{xxxx" << bitset<m>(it->first) <<  "} = " << (it->second).first / (float) N / (2*2*2*2) << endl;

    for (it_list = (it->second).second.begin(); it_list != (it->second).second.end(); ++it_list)
    {
      cout << "\t " << (*it_list).first << ": " << bitset<n>((*it_list).first) << " : \t " << (*it_list).second << "pts" << endl;
    }
    mu_control++;
  }
  cout << endl;

  if (Ncontrol != N) { cout << "Error read_Kset function: Ncontrol != N" << endl;  }
}

void distrib_Kset(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset)
{
  map<unsigned int, unsigned int> histo_set;    // histo_set[k] = #of states that appears k times in the data set
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it;
  int Ncontrol = 0;

  cout << "--->> Build distribution of Kset..." << endl;

//Build Histogram:
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    histo_set[(it->second).first] += 1;
    Ncontrol += (it->second).first;     
    //cout << it->first << ":  " << bitset<n>(it->first) << " => " << (it->second).first << '\n';
  }
  if (Ncontrol != N) { cout << "Error distrib_Kset function 1: Ncontrol != N" << endl;  }

//Print Histogram into a file:
  string file_name = "distrib_states_m"+to_string(m)+".dat";
  fstream fichier(file_name.c_str(), ios::out);
  map<unsigned int, unsigned int>::iterator histo_it;
  Ncontrol = 0;

  for(histo_it = histo_set.begin(); histo_it!=histo_set.end(); ++histo_it)
    {
      fichier << histo_it->first << " " << histo_it->second << endl;
      Ncontrol += (histo_it->first)*(histo_it->second);
    }
  fichier.close();
  if (Ncontrol != N) { cout << "Error distrib_Kset function 2: Ncontrol != N" << endl;  }

  cout << endl;
}

/******************************************************************************/
/************************** LIKELIHOOD ****************************************/
/******************************************************************************/
// Given the value of Kset --> return the Loglikelihood
double Likelihood(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset)  //
{
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it;
  double L = 0.;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second).first;
    if (Ks == 0) {cout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    L += (Ks * log((double) Ks / N) );
    Ncontrol += Ks;
  }
  if (Ncontrol != N) { cout << "Error Likelihood function: Ncontrol != N" << endl;  }

  return L;
}

/******************************************************************************/
/*****************************  MODELS ****************************************/
/******************************************************************************/
void print_basis(set<uint32_t> Operators)  // print a set of Operators, with bit representation
{
  set<uint32_t>::iterator it_op; 
  int count = 0;

  for (it_op = Operators.begin(); it_op != Operators.end(); ++it_op)
  {  
    cout << "\t" << "phi^{" << (*it_op) << "} = " << bitset<n>(*it_op);
    count++;
  } 
  cout << "\t---> " << count << " operators" << endl;
}

void print_Op(set<uint32_t> Operators)  // print a set of Operators
{
  set<uint32_t>::iterator it_op; 
  int count = 0;

  for (it_op = Operators.begin(); it_op != Operators.end(); ++it_op)
  {
    cout << (*it_op) << " ";
    count++;
  } 
  cout << "\t---> " << count << " operators" << endl;
}

//Construct all the operators of a sub-complete model from a basis
set<uint32_t> all_op(set<uint32_t> basis)
{
  set<uint32_t>::iterator it_basis;
  set<uint32_t> op;  //set of all the operators, excluded the identity operator Id

  set<uint32_t> C_Op;  //C_Op.insert(0);   //Courant Operators 
  vector< set<uint32_t> > vect_C_Op;     // vector of Courant Operators
  vector< set<uint32_t> > vect_C_Op_buffer;     // BUFFER:   vector of Courant Operators

  set<uint32_t>::iterator it_op;

  int i = 0, j = 0;
  uint32_t phi = 0;  // new operators

//******* OPERATORS OF THE BASIS: **************************
  // Phi_1 -->  1rst element of the basis
  it_basis = basis.begin();
  op.insert(*it_basis);
  ++it_basis;

  // Phi_i  for i > 1
  while(it_basis != basis.end())
  {
    op.insert(*it_basis);
    C_Op.insert(*it_basis);   vect_C_Op.push_back(C_Op);      C_Op.clear();
    ++it_basis; 
  }

  //******* ALL COMBINATIONS BASIS OPERATORS: **************************
  while(!vect_C_Op.empty())
  { // Phi_1 x (...)
    i = 0;    it_basis = basis.begin();      //phi_1 = *it_basis;
    
    while(i<vect_C_Op.size())
    {
      for(j = i; j < vect_C_Op.size(); j++)  // LOOP for a fixed PHI_(i+1) 
      {
        for(it_op = vect_C_Op[j].begin(); it_op!=vect_C_Op[j].end(); ++it_op)
        { 
          //cout << "phi_" << (i + 1) << " = " << (*it_basis) << "\t Op = " << (*it_op);
          phi = (*it_basis)^(*it_op);
          //cout << " \t phi = " << bitset<m>(*it_basis) << " ^ " << bitset<m>(*it_op) << " = " << phi <<  " = " << bitset<m>(phi)  << endl; // " phi = phi_1 ^ Op = "
          op.insert(phi);
          if(i > 0) { C_Op.insert(phi); }  //  PHI_(i+1) != PHI_1
        }
      }
      //cout << endl;
      if(i > 0) { vect_C_Op_buffer.push_back(C_Op);  C_Op.clear();  } 

    //Go to the next PHI_(i+1)  -->  PHI_(i+1) = (*it_basis);
      ++it_basis;  ++i;
    }

    vect_C_Op.clear();
    vect_C_Op.swap(vect_C_Op_buffer);
  }
  return op;
}

/******************************************************************************/
/************************  Post-treatment of the MC  **************************/
/******************************************************************************/
void print_all_FinalModel(set<FinalModel> all_Models)
{
  set<FinalModel>::iterator it_model;
  int count = 1;

  for (it_model = all_Models.begin(); it_model!=all_Models.end(); ++it_model)
  {
    cout << "\t Model " << count << ", "; count++;
    cout << "L_max = " << (*it_model).LogL << " : \t";
    print_Op((*it_model).all_op);    cout << endl;
  }
}

set<FinalModel> post_treatment_bases(set< pair< double, set<uint32_t> >, Cmp > basis_best_set) // pair(LogL, basis)
{
  set< pair< double, set<uint32_t> >, Cmp >::iterator it_set;

  FinalModel Model;   // contains all the operators of a model
  set<FinalModel> all_Models;

  int count = 0;
  it_set = basis_best_set.begin(); 
  double L_max = ((*it_set).first);

  while ( abs( (*it_set).first - L_max ) < eps && (it_set!=basis_best_set.end()) )
  {
    count++;   
    Model.LogL = (*it_set).first;     Model.basis = (*it_set).second;
    Model.all_op = all_op((*it_set).second);   // return all the operators generated by this basis

    cout << "\t "<< count << ": L_max = " << Model.LogL << ": ";  print_basis(Model.basis);
      //cout << " --->  " << count << endl;
      //cout << count << ": L = " << (*it_set).first << endl; count++;
    all_Models.insert(Model);
    ++it_set;
  }
  cout << count << " different bases in total." << endl;

  return all_Models;
}

/******************************************************************************/
/*****************************  Tools for ALGO 1  *****************************/
/******************************************************************************/
bool transform(int i, int j, uint32_t *mu) //vector<agent>& vec)
// index i and j:  si = si + sj (mod2) 
{
  uint32_t UN_j = un << (j-1); // binary nb with a single 1 bit at the j-th position
  uint32_t S_j = (*mu) & UN_j; // return the value of the j-th bit of mu

  if (S_j == 0)   { return false; }  //mu is not modified
  
  else // i-th bit of mu is flipped
  {
    // uint32_t UN_i = un << (i-1);  // binary nb with a single 1 bit at the i-th position
    *mu = (*mu) ^ ( un << (i-1) ); 
    return true;
  }
}

void OUT_transformation(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset, map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > &Kset_buffer, uint32_t i, uint32_t j)
{
  //cout << "OUT:   ";

  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it_Kset;
  list<pair<uint32_t, unsigned int> >::iterator it_list;

  //unsigned int Ks_buffer = 0;  // check

  uint32_t mu_m = 0;
  uint32_t mu = 0;
  unsigned int N_mu = 0;

  for (it_Kset = Kset.begin(); it_Kset!=Kset.end(); ++it_Kset)
  {
    mu_m = it_Kset->first;
    Kset_buffer[mu_m].first = (it_Kset->second).first;        //Ks = (it_Kset->second).first;
    //Ks_buffer = 0; // check

    for (it_list = ((it_Kset->second).second).begin(); it_list != ((it_Kset->second).second).end(); ++it_list)
    {
      mu = (*it_list).first;
      N_mu = (*it_list).second;
      if( transform(i, j, &mu) ) 
        {   Kset_buffer[mu_m].second.push_back(make_pair(mu, N_mu));    }
      else { Kset_buffer[mu_m].second.push_back(*it_list); }
      //Ks_buffer += N_mu;  // check
    }
    //if (Ks_buffer != Kset_buffer[mu_m].first) { cout << " Error: Kset_buffer of " << mu_m << endl;  }  // check
  }
} 

void IN_transformation(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset, map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > &Kset_buffer, uint32_t i, uint32_t j)
{
  //cout << "IN:    ";

  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it_Kset;
  list<pair<uint32_t, unsigned int> >::iterator it_list;

  uint32_t mu_m = 0;
  uint32_t mu = 0;
  unsigned int N_mu = 0;
  uint32_t UN_i = un << (i-1);

  for (it_Kset = Kset.begin(); it_Kset!=Kset.end(); ++it_Kset)
  {
    mu_m = it_Kset->first;      //cout << "i = " << i << ":\t " << "mu_m = " << bitset<n>(mu_m) << " = " << mu_m;
    if (!transform(i, j, &mu_m)) //mu_m remains unchanged after the transformation si = si+sj
    {
      Kset_buffer[mu_m] = (it_Kset->second);   //Kset[mu_m] is unchanged      //cout << " : no change --> mu_m = " << bitset<n>(mu_m) << " = " << mu_m << endl;
    }
    else  // mu_m has changed: the i-th bit has been flipped  --> flip the i-th bit of all the mu in the list, Kset[mu_m].second
    { //cout << " : change --> mu_m = " << bitset<n>(mu_m) << " = " << mu_m << endl;
      for (it_list = ((it_Kset->second).second).begin(); it_list != ((it_Kset->second).second).end(); ++it_list)    //flip bit i in the list of la the mu
      {
        mu = (*it_list).first;
        N_mu = (*it_list).second;

        Kset_buffer[mu_m].first = (it_Kset->second).first;
        Kset_buffer[mu_m].second.push_back(make_pair((mu^UN_i), N_mu));    //cout << "mu = " << bitset<n>(mu) << " ,\t mu_prim = " << bitset<n>(mu^UN_i) << endl;           
      }
    }
  }
} 

void CROSS_transformation(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset, map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > &Kset_buffer, uint32_t i, uint32_t j)
{ // i IN [1:m]  and j OUT
  //cout << "CROSS: ";

  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > >::iterator it_Kset;
  list<pair<uint32_t, unsigned int> >::iterator it_list;

  uint32_t mu_m = 0;
  uint32_t mu_m_flip_i = 0;
  uint32_t mu = 0;
  unsigned int N_mu = 0;
  uint32_t UN_i = un << (i-1);

  for (it_Kset = Kset.begin(); it_Kset!=Kset.end(); ++it_Kset)
  {
    mu_m = it_Kset->first;
    mu_m_flip_i = mu_m^UN_i; //i-th bit of mu_m is flipped
    //cout << "i = " << i << ":\t " << "mu_m = " << bitset<n>(mu_m) << " = " << mu_m;

    for (it_list = ((it_Kset->second).second).begin(); it_list != ((it_Kset->second).second).end(); ++it_list)    //flip bit i in the list of la the mu
    {
      mu = (*it_list).first;
      N_mu = (*it_list).second;

      if (!transform(i, j, &mu)) // mu remains unchanged, i.e. mu_m remains unchanged
      {  
        Kset_buffer[mu_m].first += N_mu;
        Kset_buffer[mu_m].second.push_back(*it_list);  //unchanged
      }
      else // mu has changed: the i-th bit has been flipped  --> i-th bit of mu_m must be flipped too
      { //new mu_m = mu_m^UN_i
        Kset_buffer[mu_m_flip_i].first += N_mu;
        Kset_buffer[mu_m_flip_i].second.push_back(make_pair(mu, N_mu));  //unchanged        
      }         
    }
  }
} 

/******************************************************************************/
/*****************************  ALGO 1  ***************************************/
/******************************************************************************/
    // sample two indices i and j
    // perform si = si + sj (mod2) on the data set (Kset)
    // the basis transformation is performed later if the new model is accepted by the MC
set<FinalModel> Algo_MC1(map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset, double *L_max)
{
// VARIABLES and BUFFER:
  uint32_t i, j;  rand_2index(&i, &j);
  uint32_t UN_i_buffer = un;
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset_buffer; // contains the temporary new Kset
  bool shuffling = false;

// Initial BASIS: (phi_1 to phi_n)
  uint32_t *basis = (uint32_t *)malloc(n * sizeof(uint32_t));
  for (int k = 0; k<n; k++) {  basis[k] = UN_i_buffer;  UN_i_buffer = UN_i_buffer << 1;  }
  //for (int k = 0; k<m; k++) {  cout << "Op[" << k << "] = " << bitset<n>(basis[k]) << " = " << basis[k] << endl; }

//********** PRINT INTO FILE ****************
  string OUTPUT_file_name = MC_fileOUT + "_m" + to_string(m) + "_B" + to_string_with_precision(beta, 3);
         OUTPUT_file_name += "_Nit" + to_string_scientific(N_MCsample) + ".dat";
  fstream fichier(OUTPUT_file_name.c_str(), ios::out);
  fichier << "## 1:MC_time \t 2:L \t 3:basis" << endl;

//********** PRINT TO TERM *******************
  cout << endl << "--->> General info:";
  cout << "\t N_MCsteps = " << N_MCsample << " (including shuffling)\t beta = " << beta << endl;
  cout << "\t\t\t Output file = \"" << OUTPUT_file_name << "\"" << endl;

  cout << endl << "--->> Run MC Algo1..." << endl;
  cout << "\t Record basis of the current complete model each time L_max increases, ";
  cout << "\t where L_max = max(LogL) + 4 N log(2.)  :" << endl << endl;
  cout << "\t All n=" << n << " operators == basis of the current data transformation." << endl;
  cout << "\t First m=" << m << " operators  == Basis of the complete model" << endl << endl;

// INITIALISATION:
  double L = Likelihood(Kset);
  double L_new = L;
  Kset_buffer.clear();

// MODEL SELECTION, RECORD MODELS with MAX-LIKELIHOOD:
  //double 
  *L_max = L;
  int L_max_round = ( (int) (*L_max) ) - 1; 
  set<uint32_t> basis_best;
  set< pair< double, set<uint32_t> >, Cmp > basis_best_set;  //likelihood + set of model bases

// START RECORDING:
  cout << "  t = 0     \t L_max = " << L << ": \t";
  fichier << 0 << "\t" << L ;   //print basis in file
  for (int k=0; k<m; k++) 
  {  
    fichier << "\t" << basis[k];
    basis_best.insert(basis[k]);
    cout << "\t" << basis[k] << ": " << bitset<n>(basis[k]);
  }
  cout << "\t  ||";
  for (int k=m; k<n; k++) 
    {    cout << "\t" << basis[k] << ": " << bitset<n>(basis[k]);   }
  fichier << endl;  cout << endl;

  basis_best_set.insert( make_pair(*L_max, basis_best) );

//*************************************
//**********  MC LOOP: ****************
//*************************************
  int count_Cross_UP = 0;  // nb of cross transf UP (Lnew > L) --> always accepted
  int count_shuffling = 0;

  for (int i_MC = 1; i_MC <= N_MCsample; i_MC++)
  {
    //select 2 indices:
    if(count_shuffling > N_shuffling) {  rand_2indexm(&i, &j);  count_shuffling = 0; }
    else {  rand_2index(&i, &j);  }
    //cout << "i,j = " << i << "  " << j << "\t";

    shuffling = false;
    if (i > m)  // si = si+sj ;  Nset is modified, but not Kset (so is L unmodified)
      {  //cout << "OUT  (" << i << "," << j << ")"<< endl; 
      OUT_transformation(Kset, Kset_buffer, i, j);  shuffling = true;  count_shuffling++;  L_new = L;  }
    else // i in [1:m]
    {
      if (j > m)  // si = si+sj ;  Nset, Kset and L are modified 
        {  CROSS_transformation(Kset, Kset_buffer, i, j);      L_new = Likelihood(Kset_buffer);    }
      else // (i,j) in [1:m]   // si = si+sj ;  Nset is modified, Kset is shuffled, but L remains unchanged
        {  IN_transformation(Kset, Kset_buffer, i, j);  shuffling= true;  count_shuffling++; L_new = L;  }
    }
  //********** DATA FOR USER ****************  
    //cout << "L = " << L_new << "\t";  //Likelihood
    //read_Kset(Kset_buffer);

  // SWAP Kset and Kset_buffer:
    if (shuffling || L_new >= L || (drand48() <= exp(beta*(L_new - L)) ) )   //accept the move
    {
      //CHANGE BASIS:
      basis[(i-1)] = (basis[(i-1)])^(basis[(j-1)]);   // sig_i_prim = sig_i + sig_j 

      //SAVE NEW BASIS:
      if (L_new > L)
      { 
        count_Cross_UP++; //cout << "Accepted: ";
        if (L_new >= L_max_round)
        {
          *L_max = L_new;
          L_max_round = ( (int) (*L_max) ) - 1;
          cout << "  t = " << i_MC << "   \t L_max = " << (*L_max) << ": \t"; //<< "; L_max_round = " << L_max_round << endl;

          basis_best.clear();
          for (int k = 0; k<m; k++)   //take the m first basis elements
          {  
            basis_best.insert(basis[k]);
            cout << "\t" << basis[k] << ": " << bitset<n>(basis[k]);
          }
          cout << "\t  ||";
          for (int k = m; k<n; k++)
            { cout << "\t" << basis[k] << ": " << bitset<n>(basis[k]); }
          cout << endl;

          if ( basis_best.size() != m ) 
            {  cout << "Error basis:";  for (int k = 0; k<n; k++) {  cout << "\t" << basis[k] << ": " << bitset<n>(basis[k]);  } }
        }
        basis_best_set.insert( make_pair(*L_max, basis_best) );
      }
      Kset.clear();    Kset.swap(Kset_buffer); 
      L = L_new;
      //cout << "\t A: L = " << L; //accept   

      //PRINT BASIS IN FILE:
      fichier << i_MC << "\t" << L;
      for (int k = 0; k<n; k++)
        {   fichier << "\t" << basis[k];    }
      fichier << endl;
    }
    else //reject the move
      {   Kset_buffer.clear();    
          //cout << "\t R: L = " << L; //reject
      }
  }
  fichier.close();

  cout << endl << "\t # of Cross transf UP (always Accepted) = " << count_Cross_UP << endl;
  cout << "\t # of Cross transf UP that gives a new possible best complete model = " << basis_best_set.size() << endl;

//*************************************
//*********  POST-TREATMENT: **********
//*************************************
  cout << endl << "--->> Post-treatement:\t all bases found for L_max at eps = " << eps << " from L_max = " << (*L_max) << ":" << endl;
  set<FinalModel> all_Models = post_treatment_bases(basis_best_set);

  cout << endl << "--->> Post-treatement:\t MC Algo1 found " << all_Models.size() ;
  cout << " model(s) at eps = " << eps << " from L_max = " << (*L_max) << ":" << endl << endl;
  print_all_FinalModel(all_Models);
  return all_Models;
}

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/
//Rem : ---> 2^30 ~ 10^9

int main()
{
//*******************************************
//*********** INITIALISATION  ***************
//*******************************************
  init_rand();
  bool print_bool = false;  // true  --> print in the terminal the construction of Nset, Kset, etc...
                            // false --> don't print

//*******************************************
//********** READ DATA FILE:  ***************     -->  data are put in Nset:
//*******************************************    // Nset[mu] = # of time state mu appears in the data set
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "***********************************  Read the data:  **************************************";
  cout << endl << "*******************************************************************************************" << endl;

  map<uint32_t, unsigned int> Nset = read_datafile();

//*******************************************
//********** PRINT INFORMATIONS:  ***********
//*******************************************
  // number m of relevant spins   // Rem: m must be strictly smaller than 32, i.e. N<=10^18, which is likely the case!
  unsigned int m_advised = (int) (0.5*log(N)/log(2.)) + 1;   // advised m obtained by inverting: 2^m - 1 ~ sqrt(N)

  cout << endl << "--->> General info: \t";
  cout << "n=" << n << ", \t m = " << m ;                                                     //number of spins
  cout << "\t \t [ Rem: advised m = " << m_advised ;                                          //advised m                                                     //chosen m
  cout << ", \t Nop = 2^m-1 = " << Nop << " ~ sqrt(N) ~ " << sqrt(N) << " ]" << endl; //number of operators

//*******************************************
//********** STATISTICS OF Nset:  ***********
//*******************************************
  cout << endl << "--->> Statistics of Nset:";
  cout << "\t N = " << N << " = nb of points in the data set" << endl;
  cout << "\t\t\t\t Ndiff = " << Nset.size() << " = nb of differents states appearing in the data set" << endl;
  cout << "\t\t\t\t N/Ndiff = " << float (N)/(Nset.size()) << " = average nb of times a states appear" << endl;
  //distrib_Nset(Nset);
  if (print_bool) { 
  read_Nset(Nset);  }   // print Nset in the terminal

//*******************************************    // mu_m = mu cut on the n first spins 
//******* BUILD and STATS of Kset:  *********
//*******************************************
  map<uint32_t, pair<unsigned int, list<pair<uint32_t, unsigned int> > > > Kset; 
    // Kset[mu_m].first = # of time state mu_m appears in the data set
    // Kset[mu_m].second = list of pair<mu, N_mu>, 
    //                          where mu are the states such that mu|m = mu_m
    //                          and N_mu is the # of times state mu appears in the data set
  Kset = build_Kset(Nset, false);  // true  --> print the construction of Kset (while building it)
  //distrib_Kset(Kset);
  if (print_bool) { 
  cout << "Initial ";   read_Kset(Kset);  }  // print Kset in the terminal

//*******************************************
//******* initial LOG-LIKELIHOOD:  **********
//*******************************************
  cout << "--->> Initial logLikelihood: \t";
  double LogL = Likelihood(Kset);
  cout << "logL + (n - m) N log(2.) = " << LogL << "\t\t logL = " << LogL - (n - m)*N*log(2.) ;
  cout << "\t\t logL/(-n N log(2)) = " << (LogL - (n - m)*N*log(2.))/(-1.*n*N*log(2.)) << endl;

//*******************************************
//********** MONTE CARLO ALGO 1:  ***********
//*******************************************
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "***********************************  Monte Carlo:  ****************************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << "--->> Create folder: ";
  system( ("mkdir " + directory).c_str() );
  cout << endl;
 
  set<FinalModel> all_Models = Algo_MC1(Kset, &LogL);
  cout << "--->> Common MaxLogLikelihood for the best model(s) found: \t" << endl;
  cout << "\t\t logL_max + (n - m) N log(2.) = " << LogL << "\t\t logL_max = " << LogL - (n-m)*N*log(2.) ;
  cout << "\t\t logL_max/(-n N log(2)) = " << (LogL - (n-m)*N*log(2.))/(-n*N*log(2.)) << endl << endl;

  return 0;
}


