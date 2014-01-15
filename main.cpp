#include <pegasos-svm.h>
#include <iostream>

int main()
{
    //std::vector<double> y;
    //std::vector<std::vector<double>> S;

    /*
    if(argc < 2)
    {
        std::cout << "Usage: pegasos iterations learning_rate type [gamma]" << std::endl;
        std::cout << " type : 0 - linear, 1 - RBF kernel " << std::endl;
        return 0;
    }
    */

    auto data=load_data("nor.data");   
 
    print_vector(data.first);
    //print_matrix(data.second);

    std::vector<double> v1 = {1,1,1,1}; 
    std::vector<double> v2 = {1,2,3,4};

    print_vector(v1*3);
    std::cout <<  v1*v2 << std::endl;
    std::cout <<  (v1*3)*v2 << std::endl; 

    print_vector(v1 + v2);
    print_vector(v1 - v2);

    std::cout << norm(v1) << std::endl;   

    return 0;
}
