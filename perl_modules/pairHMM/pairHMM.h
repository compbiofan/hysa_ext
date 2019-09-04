#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string> 
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>

using namespace std;

string version ("pairHMM-0.0");
string options ("");

double smallNum1 = -10e10;
double smallNum2 = -10e12;
double modification_error = 0.1;
double log_mod = log(modification_error);
double log_mod_minus = log(1-modification_error);
double log_p25 = log(0.25);

class tools{
public:

    void split(string rReads, string sep, vector<string> &splitted);
    string chomp(string str);
};


void print_aln(const vector<char> &aln, const string seqAF, const string seqBF);

void pairAln(const string seqAF, const string seqBF, const string transF);

void this_max3(double tmp1, double tmp2, double tmp3, char c1, char c2, char c3, char &select, double &score);

void this_max2(double tmp1, double tmp2, char c1, char c2, char &select, double &score);

void this_max_array(const vector<double> &tmp, const vector<char> &c, char &select, double &score);

void max_array(const vector<double> &array, char &c, double &score);

double ratio(char x, char y);

double viterbi(const string seqAF, const string seqBF, const double *th, vector<char> &aln);

string read_seq(const string seqF);

void read_trans(const string transF, double *th);
