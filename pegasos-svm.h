#ifndef _PEGASOS_SVM
#define _PEGASOS_SVM

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
#include <vectorlib.h>
#include <utility>

std::pair<std::vector<double>, std::vector<std::vector<double>>> load_data (const char* s) 
{
    std::ifstream file;
    file.open(s);
    std::string line;
    std::istringstream instream;
    double val;

    std::vector<std::vector<double> > S;
    std::vector<double> labels;     

    // Only raw format and svmlight format are supported
    bool is_raw_format = true;    
    bool is_first_line = true;
    int line_count = 1;    

    while ( file.good() )
    {
        getline (file,line);
        std::vector<double> row;

        instream.clear();
        instream.str(line);

        if(is_first_line)
        {
            is_first_line = false;

            std::string::size_type loc = line.find( ":", 0 );

            if( loc != std::string::npos )
            {
                is_raw_format = false;
            }
        }


        if(line.length() > 0)
        {
            instream >> val;
            labels << val;

            if(is_raw_format)
            {

                while(!instream.eof())
                {
                    instream >> val;
                    row << val;
                }

                if(row.size()>0)
                    S.push_back(row);
            }
            else
            {
                // The file is parsed assuming a svm-light file format

                std::string item;
                int i=0;

                while (instream >> item) 
                {

                    if (item != "" )
                    {
                        int index;
                        double feature_val;
                        std::string strindex,strval;

                        std::stringstream ss(item);
                        std::getline(ss,strindex,':');
                        std::getline(ss,strval,':');

                        std::stringstream ssindex(strindex);
                        std::stringstream ssval(strval);

                        ssindex >> index;
                        ssval >> feature_val;

                        while ( i < (index-1) )
                        {
                            row << 0;
                            i++;
                        }

                        
                        row << feature_val;
                        i++;
                    }
                    //std::cout << "!" << item << "!";
                }

                if(row.size()>0)
                {
                    S.push_back(row);
                }

                //std::cout << "line:" << line_count << std::endl;
            }
        }
        line_count++;
    }

    file.close();
    return make_pair(labels,S);
}


//
//
// MODEL
//
//
//


class RBF
{
    private:
        double _gamma;

    public:
        double operator()(const std::vector<double>& v1, const std::vector<double>& v2)
        {
            double result;
            double _gamma2 = pow(_gamma, 2.0);

            result = exp(-1*(pow(norm(v1-v2), 2)/(2*_gamma2)));

            return result;
        }

        RBF(double gamma)
        {
            _gamma=gamma;
        }


};

class Model
{
	private:
		bool useKernel;
	public:
		Model(std::vector<double> w):useKernel(false)
		{
			_w = w;
		}
		
		Model(std::vector<double> alpha, double lambda ,int T, std::vector<double> _y, std::vector<std::vector<double> > S, double gamma):useKernel(true)
		{
			_alpha = alpha;
			coef = (1/(lambda*T));
			_t = T;
			y = _y;
			_x = S;
			_gamma = gamma;
			m = S.size();
		}

		double operator()(std::vector<double> x)
		{
			if (useKernel)
				return kpredict(x);
			else
				return predict(x);
		}

		void printWeights()
		{
			print_vector(_w);
		}

	protected:
		std::vector<double> _w;
				
		double predict(std::vector<double> x)
		{
			double val = _w*x;
			double result = (val < 0) ? -1 : 1;
			return result;
		}
		
		double kpredict(std::vector<double> x)
		{
			double val = coef;
			double sum =0;
			RBF K=RBF(_gamma);

			for(int j=0; j< m; j++)
			{
				sum += _alpha[j]*y[j]*K(x,_x[j]);
				//sum += _alpha[j]*K(x,_x[j]);
			}

			val = val * sum;

			double result = val; 
			return (result>0)?1:-1;
		}

		std::vector<double> _alpha;
		std::vector<double> y;
		std::vector<std::vector<double> > _x;
		double coef;
		int _t;
		double _gamma;
		int m;

};



//
//
//
// LINEAR SOLVER
//
//
//
//



class Pegasos
{
    protected:
        std::vector<double> _error;
        const static int error_samples = 4;

    public:
				static Model train(std::vector< std::vector < double > > S, std::vector<double> y, double lambda, int T  )
				{
					Pegasos p;
					return p._train(S,y,lambda,T);
				}

		private:
				Pegasos() {}
				~Pegasos() {}
	
				Model _train(std::vector< std::vector < double > > S, std::vector<double> y, double lambda, int T  )
        {

            int m = S.size();
            int n = S[0].size();
            std::vector<double> w(n);

            for ( int t = 1; t< T+1 ; t++)
            {
                double nt = 1 / (lambda * (double) t);
                double f1 = (1.0-nt*lambda);
                int i = rand() % m;
                std::vector<double>& x = S[i];
                //PM@TODO Are values in rate initialized in zero?
                //not working.. argh!
                std::vector<double> rate(n);
                print_vector(rate);
                double predict =  y[i] * ( w * x );
                
                if ( predict < 1)
                {
                    rate = (x*y[i])*nt;
                }
                
                w = w*f1; 
                w = w + rate;
                dump_error(S, y, lambda, t, w);
            }
            //print_vector(_error);
            return Model(w);
        }

    void dump_error(std::vector<std::vector<double> > &S, std::vector<double>& y, double lambda, int t, std::vector<double> w )
    {
        double error =0;
        double value = 0;
        int n = S.size();

        // Stochastic average error
        for ( int i=0; i < error_samples; i++ )
        {
            int r = rand() % n;
            value = w*S[r];

            error = error + pow(y[r] - value, 2)/2;
        }

        error = error / error_samples;

        _error << error;
    }

};

///
//
/// Kernelized version
//
///


class PegasosKernelized
{
    private:
        std::vector<double> _error;
        const static int error_samples = 4;
				PegasosKernelized() {}
        ~PegasosKernelized() {}
				
				Model _train(std::vector< std::vector < double > > S, std::vector<double> y, double lambda, int T, double gamma  )
        {

            int m = S.size();
            int n = S[0].size();
            std::vector<double> _alpha(m);
            RBF K=RBF(gamma);

            for ( int t = 1; t< T+1 ; t++)
            {
                int it = rand() % m;
                double nt = 1 / (lambda * (double) t);
                double f1 = (1.0-nt*lambda);

                double sum =0;

                for (int j=0; j < m; j++ )
                {
                    sum += _alpha[j]*y[it]*K(S[it], S[j]); 
                }


                if ( y[it]*nt*sum < 1 )
                {
                    _alpha[it] = _alpha[it] + 1;
                }

                dump_error(S, y, lambda, t, _alpha,gamma);
            }
            std::cout << "alpha: ";
            print_vector(_alpha);
            return Model(_alpha,lambda,T,y,S,gamma);
        }


        void dump_error(std::vector<std::vector<double> > &S, std::vector<double>& y, double lambda, int t, std::vector<double> w, double gamma )
        {
            double error =0;
            double value = 0;
            int n = S.size();

            // Stochastic average error
            for ( int i=0; i < error_samples; i++ )
            {
                int r = rand() % n;

                Model m = Model(w,lambda,t,y,S,gamma);
                value = m(S[r]);

                error = error + pow(y[r] - value, 2)/2;
            }

            error = error / error_samples;
            _error << error;
        }

    public:
			static Model train(std::vector< std::vector < double > > S, std::vector<double> y, double lambda, int T, double gamma  )
			{
				return PegasosKernelized()._train(S, y, lambda, T, gamma);
			}
};


class PegasosTrainer 
{
    public:
 
    Model train(const std::vector<std::vector<double>> &TS, std::vector<double> labels, double learning_rate,int iterations, bool kernel,  double gamma=-1.0)
    {
        
        if(!kernel)
        {
            return Pegasos::train(TS, labels, learning_rate, iterations );
        }
        else
        {
            return PegasosKernelized::train(TS, labels, learning_rate, iterations, gamma );
        }   
    }
};


#endif
