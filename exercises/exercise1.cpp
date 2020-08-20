#include <typeinfo>
#include<iostream>
#include<cmath>

using namespace std;

int array_size(int before_decimal, int number_of_digits){
    int biggest_number = pow(10,(before_decimal));
    int binary_digits = 0;
    int binary_number = 1;
    while (binary_number<biggest_number){
        binary_number  = binary_number*2;
        binary_digits++;
    }
    return binary_digits + number_of_digits;
};

int* allocate_array(int size){
    int* digits = new int[size];
    for (int i = 0; i<size; i++){
        digits[i] = 0;
    }
    return digits;        
};

void find_binary(int* digits, int whole_part, float float_part, int size, int binary_digits){

    
    int a = whole_part;
    for (int i = binary_digits; i<=size; i++){
        int digit_i = a%2;
        cout << digit_i << endl;
        a = a/2;
        digits[i] = digit_i;
    }

    float b = float_part;
    for (int i = binary_digits-1; i>=0; i--){
        int digit_i = int(floor(b*2));
        b = b*2 - digit_i;
        digits[i] = digit_i;
    }
};





//lese antall siffer før desimaltegnet og finne største binær-tallet.


int main(int argc, char** argv){

    int binary_digits = atoi(argv[2]);
    float number_to_convert = atof(argv[1]);
    int whole_part = int(number_to_convert);
    float float_part = number_to_convert - whole_part;
    char* s = argv[1];
    float number_decimal = atof(argv[1]);
    bool decimal_found = 0;
    int length_of_argument = 0;
    int before_decimal = 0;
    int after_decimal = 0;

    for (char* t = s; *t; t++) {
        char c = t[0];
        if (c == '.'){
            decimal_found = true;
            continue;
        }
        if (decimal_found == false){
            before_decimal++;
        }else{
            after_decimal++;
        }
    }
    int size = array_size(before_decimal, binary_digits);
    cout << size << endl;
    int* digits = allocate_array(size);
    find_binary(digits, whole_part, float_part, size, binary_digits);


    for (int i = size-1; i>=0; i--){
        cout << digits[i];
        if (binary_digits == i){
            cout << ".";
        }
    }
    cout << endl;

    delete digits;
}