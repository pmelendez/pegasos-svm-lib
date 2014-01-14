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
    print_matrix(data.second);

    return 0;
}
