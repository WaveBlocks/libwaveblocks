#include <iostream>

int number(int dim, int sum)
{
    if (dim == 0 || sum = 0) {
        return 1;
    }
    else {
        int n = 0;
        for (int i = 0; i < sum; i++) {
            n += number(dim-1, i);
        }
        return n;
    }
}

int main()
{
    std::cout << number(2,2) << std::endl;
    
    return 0;
}