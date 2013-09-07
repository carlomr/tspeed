#include"Dunavant.hpp"
#include<iostream>
int main()
{
    
    int rule =4;
	int order_num = dunavant_order_num ( rule );

	double *xytab = new double[2*order_num];
	double *wtab = new double[order_num];

	dunavant_rule ( rule, order_num, xytab, wtab );
	for(int i = 0; i<2*order_num;++i)
	    std::cout << xytab[i] << std::endl;

}

