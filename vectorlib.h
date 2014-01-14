#ifndef VECTORLIB
#define VECTORLIB

#include <vector>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include <algorithm>


const int min_int = std::numeric_limits<int>::min();
const int max_int = std::numeric_limits<int>::max();
const double min_double = std::numeric_limits<double>::min();
const double max_double = std::numeric_limits<double>::max();


//
//
//
/// Utility functions
//
//
//

void print_vector(std::vector<double> vec)
{
    for(int i=0;i<vec.size();i++)
        std::cout << vec[i] << " ";
    std::cout << std::endl;
}

void print_matrix(std::vector<std::vector<double > > vec)
{
    for(int i=0; i< vec.size(); i++)
    {
        for(int j =0; j < vec[i].size(); j++)
        {
            std::cout << vec[i][j] << " ";
        }
        std::cout << "." << std::endl;
    }
}


double norm(const std::vector<double>& op)
{
    double result;

    for(int i =0; i< op.size(); i++)
    {
        result += pow(op[i],2);
    }

    result = sqrt(result);

    return result;
}

const std::vector<double>& operator << (const std::vector<double>& op, const double s) 
{
    ((std::vector<double>&) op).push_back(s);
    return op;
}

const std::vector<double> operator * (const std::vector<double>& op, const double s) 
{
    std::vector<double> newvec;
    
    for(int i=0;i<op.size(); i++)
    {
        newvec.push_back(op[i]*s);
    }

    return newvec;
}

const std::vector<double> operator + (const std::vector<double>& op, const std::vector<double>& v) 
{
    std::vector<double> newvec;
    double res;
    for(int i=0;i<op.size(); i++)
    {
        newvec.push_back(op[i]+v[i]);
    }
    return newvec;
}

const std::vector<double> operator - (const std::vector<double>& op, const std::vector<double>& v) 
{
    std::vector<double> newvec;
    double res;
    for(int i=0;i<op.size(); i++)
    {
        newvec.push_back(op[i]-v[i]);
    }
    return newvec;
}

const double operator * (const std::vector<double>& op, const std::vector<double>& v2)
{
    double res = 0;

    if(op.size() != v2.size())
        return min_double;

    for(int i =0; i< op.size(); i++)
    {
        res += op[i]*v2[i];
    }

    return res;
}


#endif
