#include "pairHMM.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string> 
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>

using namespace std;

void tools::split(string rReads, string sep, vector<string> &splitted) {
    size_t pos_start = 0;
    size_t pos;
    if(sep.length() > 0) {
        do {
            pos = rReads.find(sep, pos_start);
            if(pos != string::npos) {
                if(pos - pos_start > 0) {
                    string str_tmp = rReads.substr(pos_start, pos - pos_start);
                    splitted.push_back(str_tmp);
                }
                rReads = rReads.substr(pos + 1);
            }
            else if(rReads.length() != 0)
                splitted.push_back(rReads);
        } while (pos != string::npos);
    }
    else {
        for(int i = 0; i < int(rReads.length()); i++) {
            splitted.push_back(rReads.substr(i,1));
        }
    }

    return;
}

// move the last endl sign, and return back the string
string tools::chomp(string str) {
    string ret;
    string::size_type pos = str.find_last_not_of("\n");
    if(pos != string::npos)
        ret = str.substr(0, pos+1);
    else
        ret = str;
    return ret;
}


int main(int argc, char *argv[])
{

    int c;
    cerr << "pairHMM version-" << version << endl;
    while((c = getopt(argc, argv, "m:")) >= 0) {
        switch(c) {
            case 'm':
                modification_error = atof(optarg);
                log_mod = log(modification_error);
                log_mod_minus = log(1-modification_error);
                break;
            default:
                fprintf(stderr, "Unrecognized option '-%c'.\n", c);
        }
    }

    if(optind == argc){
        fprintf(stderr, "\n./pairHMM [options] <sequence1.fa> <sequence2.fa> <transition_matrix_file>\n");
        fprintf(stderr, "\nOptions: \n");
        fprintf(stderr, "   -m DOUBLE  Modification error penalty [%.2f]\n", modification_error);
        return EXIT_SUCCESS;
    }
    else{
        cerr << "#command:"; 
        for (int idx=0; idx<argc; idx++) cerr << " " << argv[idx];
        cerr << endl;
    }

    // read sequences and transition matrix
    string seqAF = argv[optind];
    string seqBF = argv[optind+1];
    string transF = argv[optind+2];

    if(transF.compare("NA") == 0){
        transF = "/scratch/bcb/xfan3/p/pairHMM/hmm_trans.txt";
    }

    pairAln(seqAF, seqBF, transF);

}

void print_aln(const vector<char> &aln, const string seqAF, const string seqBF){
    // print the alignment composed of 0, 1, and 2
    for(int i = aln.size() - 1; i >= 0; i --){
        printf("%d", aln[i]);
    }
    cout << endl;
    string seqA = read_seq(seqAF);
    string seqB = read_seq(seqBF);
    int pA = 0;
    int pB = 0;
    // print the second line (X)
    for(int i = aln.size() - 1; i >= 0; i --){
        if(aln[i] == 0 || aln[i] == 1 || aln[i] == 3){
            cout << seqA.at(pA ++);
        }
        else{
            cout << "-";
        }
    }
    cout << endl;
    // print the second line (Y)
    for(int i = aln.size() - 1; i >= 0; i --){
        if(aln[i] == 0 || aln[i] == 2 || aln[i] == 4){
            cout << seqB.at(pB ++);
        }
        else{
            cout << "-";
        }
    }
    cout << endl;
    return;
}

void pairAln(const string seqAF, const string seqBF, const string transF){
    double *th = new double[69];
    for(int i = 0; i < 69; i ++){
        th[i] = 0;
    }
    read_trans(transF, th);
    vector<char> aln;
    double score = viterbi(seqAF, seqBF, th, aln);
    print_aln(aln, seqAF, seqBF);
    cout << "score: " << score << "\n";
    delete []th;
    return;
}

//void this_max_array(double *tmp, char *c, char num, char &select, double &score){
void this_max_array(const vector<double> &tmp, const vector<char> &c, char &select, double &score){
    int num = tmp.size();
    score = smallNum2;
    for(int i = 0; i < num; i ++){
        if(tmp[i] > score){
            score = tmp[i];
            select = c[i];
        }
    }
    return;
}

void this_max3(double tmp1, double tmp2, double tmp3, char c1, char c2, char c3, char &select, double &score){
    if(tmp1 >= tmp2 && tmp1 >= tmp3){
        select = c1;
        score = tmp1;
    }
    else if(tmp2 >= tmp1 && tmp2 >= tmp3){
        select = c2;
        score = tmp2;
    }
    else if(tmp3 >= tmp1 && tmp3 >= tmp2){
        select = c3;
        score = tmp3;
    }
    return;
}

void this_max2(double tmp1, double tmp2, char c1, char c2, char &select, double &score){
    if(tmp1 >= tmp2){
        select = c1;
        score = tmp1;
    }
    else{
        select = c2;
        score = tmp2;
    }
    return;
}

void max_array(const vector<double> &array, char &c, double &score){
    score = smallNum2;
    for(int i = 0; i < int(array.size()); i ++){
        if(score < array[i]){
            score = array[i];
            c = char(i);
        }
    }
    return;
}

double ratio(char x, char y){
    if(x == y){
        return log_mod_minus;
    }
    else{
        return log_mod;
    }
}

double viterbi(const string seqAF, const string seqBF, const double *th, vector<char> &aln){
    string seqA = read_seq(seqAF);
    string seqB = read_seq(seqBF);
    int lX = seqA.length();
    int lY = seqB.length();

    // intiialization
    int lZ = 5;
    double*** s;
    char*** p;
    s = new double**[lX];
    p = new char**[lX];
    for(int i = 0; i < lX; i ++){
        s[i] = new double*[lY];
        p[i] = new char*[lY];
        for(int j = 0; j < lY; j ++){
            s[i][j] = new double[lZ];
            p[i][j] = new char[lZ];
            for(int k = 0; k < lZ; k ++){
                s[i][j][k] = smallNum1;
                p[i][j][k] = -1;
            }
        }
    }
    s[0][0][0] = 0;

    vector<double> tmpv(lZ, 0);
    char c[] = {'M', 'X', 'Y', 'D', 'I'};
    vector<char> c1(c, c+sizeof(c));
    // recurrence
    for(int i = 1; i < lX; i ++){
        for(int j = 1; j < lY; j ++){
            // all others would have move from M, X and Y
            double log_odd = ratio(seqA.at(i), seqB.at(j));
            // level 1

            tmpv[0] = s[i-1][j-1][0] + th[0] + log_odd;
            tmpv[1] = s[i-1][j-1][1] + th[8] + log_odd;
            tmpv[2] = s[i-1][j-1][2] + th[16] + log_odd;
            tmpv[3] = s[i-1][j-1][3] + th[32] + log_odd;
            tmpv[4] = s[i-1][j-1][4] + th[64] + log_odd;
            //this_max3(tmp1, tmp2, tmp3, 'M', 'X', 'Y', p[i][j][0], s[i][j][0]);
            this_max_array(tmpv, c1, p[i][j][0], s[i][j][0]);

            // level 2
            tmpv[0] = s[i-1][j][0] + th[1] + log_p25;
            tmpv[1] = s[i-1][j][1] + th[9] + log_p25;
            tmpv[2] = s[i-1][j][2] + th[17] + log_p25;
            tmpv[3] = s[i-1][j][3] + th[33] + log_p25;
            tmpv[4] = s[i-1][j][4] + th[65] + log_p25;
            this_max_array(tmpv, c1, p[i][j][1], s[i][j][1]);
            //this_max2(tmp1, tmp2, 'M', 'X', p[i][j][1], s[i][j][1]);

            // level 3
            tmpv[0] = s[i][j-1][0] + th[2] + log_p25;
            tmpv[1] = s[i][j-1][1] + th[10] + log_p25;
            tmpv[2] = s[i][j-1][2] + th[18] + log_p25;
            tmpv[3] = s[i][j-1][3] + th[34] + log_p25;
            tmpv[4] = s[i][j-1][4] + th[66] + log_p25;
            this_max_array(tmpv, c1, p[i][j][2], s[i][j][2]);
            //this_max2(tmp1, tmp2, 'M', 'Y', p[i][j][2], s[i][j][2]);

            // level 4
            tmpv[0] = s[i-1][j][0] + th[3] + log_p25;
            tmpv[1] = s[i-1][j][1] + th[11] + log_p25;
            tmpv[2] = s[i-1][j][2] + th[19] + log_p25;
            tmpv[3] = s[i-1][j][3] + th[35] + log_p25;
            tmpv[4] = s[i-1][j][4] + th[67] + log_p25;
            this_max_array(tmpv, c1, p[i][j][3], s[i][j][3]);

            // level 5
            tmpv[0] = s[i][j-1][0] + th[4] + log_p25;
            tmpv[1] = s[i][j-1][1] + th[12] + log_p25;
            tmpv[2] = s[i][j-1][2] + th[20] + log_p25;
            tmpv[3] = s[i][j-1][3] + th[36] + log_p25;
            tmpv[4] = s[i][j-1][4] + th[68] + log_p25;
            this_max_array(tmpv, c1, p[i][j][4], s[i][j][4]);
        }
    }

    vector<double> array;
    for(int k = 0; k < lZ; k ++){
        array.push_back(s[lX-1][lY-1][k]);
    }

    double score;
    char s_last;
    max_array(array, s_last, score);
    int i = lX - 1;
    int j = lY - 1;
    aln.push_back(s_last);
    char pre;
    while(i > 0 || j > 0){
        pre = p[i][j][int(s_last)];
        if(pre == 'M'){
            pre = 0;
        }
        else if(pre == 'X'){
            pre = 1;
        }
        else if(pre == 'Y'){
            pre = 2;
        }
        else if(pre == 'D'){
            pre = 3;
        }
        else if(pre == 'I'){
            pre = 4;
        }
        if(s_last == 0){
            // match
            i --;
            j --;
        }
        else if(s_last == 1 || s_last == 3){
            // deletion
            i --;
        }
        else if(s_last == 2 || s_last == 4){
            // insertion
            j --;
        }
        s_last = pre;
        aln.push_back(s_last);
    }
    // release memory
    for(int i = 0; i < lX; i ++){
        for(int j = 0; j < lY; j ++){
            delete [] s[i][j];
            delete [] p[i][j];
        }
        delete [] s[i];
        delete [] p[i];
    }
    delete [] s;
    delete [] p;
    return score;
}

string read_seq(const string seqF){
    string ret = "";
    tools tl;
    ifstream STR;
    STR.open(seqF.c_str());
    char line_[1000000];
    if(STR.is_open()){
        while(STR.good()){
            STR.getline(line_, 1000000);
            if(line_[0] == '>'){
                continue;
            }
            string line(line_);
            line = tl.chomp(line);
            ret += line;
        }
        STR.close();
    }
    return ret;
}

// th keys: 
// MM: 0000 (0)
// MX: 0001 (1)
// MY: 0010 (2)
// XM: 0100 (4)
// XX: 0101 (5)
// XY: 0110 (6)
// YM: 1000 (8)
// YX: 1001 (9)
// YY: 1010 (10)
void read_trans(const string transF, double *th){
    int m_num = 5;
    char i = 0;
    ifstream TRAN;
    TRAN.open(transF.c_str());

    tools tl;
    char line_[1000];
    if(TRAN.is_open()){
        while(TRAN.good()){
            TRAN.getline(line_, 1000);
            if(line_[0] == '#'){
                continue;
            }
            string line(line_);
            if(line.length() == 0)
                continue;
            i ++;
            line = tl.chomp(line);
            vector<string> fields;
            tl.split(line, "\t", fields);
            if(i == 1){
                for(int i = 0; i < m_num; i ++){
                    th[i] = log(atof(fields[i].c_str()));
                }
            }
            else if(i == 2){
                for(int i = 8; i < 8 + m_num; i ++){
                    th[i] = log(atof(fields[i-8].c_str()));
                }
            }
            else if(i == 3){
                cout << fields[0] << endl;
                for(int i = 16; i < 16 + m_num; i ++){
                    th[i] = log(atof(fields[i-16].c_str()));
                }
            }
            else if(i == 4){
                for(int i = 32; i < 32 + m_num; i ++){
                    th[i] = log(atof(fields[i-32].c_str()));
                }
            }
            else if(i == 5){
                for(int i = 64; i < 64 + m_num; i ++){
                    th[i] = log(atof(fields[i-64].c_str()));
                }
            }


        }
        TRAN.close();
    }
    return;
}







