#include <omp.h>
#ifndef A0_HPP
#define A0_HPP

#include <vector>


/* SOLUTION 1: Critical Section */

// template <typename T, typename Pred>
// void omp_extract_if(const std::vector<T>& in, std::vector<T>& out, Pred pred) {
//     #pragma omp parallel shared(out)
//     {
//         #pragma omp for schedule(auto)
//         for (int i = 0; i < in.size(); i++)
//         {
//             if(pred(in[i]))
//             {
//                 #pragma omp critical
//                 out.push_back(in[i]);
//             }
//         }
//     }
// } // omp_extract_if

/* SOLUTION 2: Custom Reduction */

// template <typename T>
// inline std::vector<T>&& append(std::vector<T> L1, const std::vector<T>& L2) {
// L1.insert(std::end(L1), std::begin(L2), std::end(L2));
// return std::move(L1);
// } // append

// template <typename T, typename Pred>
// void omp_extract_if(const std::vector<T>& in, std::vector<T>& out, Pred pred) {
    
//     using CustomVectorDT = std::vector<T>;

//     #pragma omp declare reduction \
//         (vector_reduce: CustomVectorDT :omp_out=append(omp_out,omp_in)) \
//         initializer(omp_priv = omp_orig)

//     #pragma omp parallel shared(out)
//     {

//         #pragma omp for schedule(auto) reduction(vector_reduce:out)
//         for (int i=0; i < in.size(); i++) 
//         {
//             if (pred(in[i]))
//             {
//                 out.push_back(in[i]);
//             }
//         }        
//     }
// } // omp_extract_if


/* SOLUTION 3: Parallel Prefix */

// template <typename T, typename Pred>
// void omp_extract_if(const std::vector<T>& in, std::vector<T>& out, Pred pred) {
    
//     int inputSize = in.size();

//     // std::cout << "Input Size : " << inputSize << std::endl;

//     int inputLogSize = std::ceil(std::log2(inputSize));
//     int bitVectorSize = 1 << inputLogSize;
//     int paddedSize = bitVectorSize - inputSize;

//     // if(paddedSize == 0)
//     // {
//     //     //Can't push into const!
//     //     in.push_back(0);
//     //     inputSize = in.size();
//     //     inputLogSize = std::ceil(std::log2(inputSize));
//     //     bitVectorSize = std::pow(2,inputLogSize);
//     // }
//     // std::cout << "Padded Size : " << paddedSize << std::endl;

//     std::vector<int> bitVector(bitVectorSize);

//     // std::cout << "Bit Vector Size : " << bitVectorSize << std::endl;

//     // for(auto it=input.begin(); it != input.end(); it++)
//     // {
//     //     std::cout <<  *it << " ";
//     // }
//     // std::cout <<  std::endl;

//     #pragma omp parallel for schedule(auto) shared(bitVector)
//     for(int i = 0; i < inputSize; i++)
//     {
//         if(pred(in[i]))
//             bitVector[i] = 1;
//         else
//             bitVector[i] = 0;
//     }
    
//     std::vector<int> bitSumVector = bitVector;

//     // for(auto it=bitVector.begin(); it != bitVector.end(); it++)
//     // {
//     //     std::cout << *it << " "; 
//     // }
//     // std::cout << std::endl;

//     int powI = 0;
//     int powI1 = 0; 
//     int temp = 0;
    
//     #pragma omp parallel shared(bitSumVector)
//     {
//         for(int i = 0; i < inputLogSize-1; i++)
//         {
//             powI = 1 << i;
//             powI1 = 1 << (i+1);
//             #pragma omp for schedule(auto) nowait
//             for(int j = 0; j < bitVectorSize-1; j += powI1)
//             {
//                 bitSumVector[j + powI1 - 1] = 
//                     bitSumVector[j +  powI - 1 ] + bitSumVector[j + powI1 - 1];
//             }
//         }   
//     }

//     // for(auto it=bitVector.begin(); it != bitVector.end(); it++)
//     // {
//     //     std::cout << *it << " "; 
//     // }
//     // std::cout << std::endl;

//     bitSumVector[bitVectorSize - 1] = 0;

//     #pragma omp parallel private(temp) shared(bitSumVector)
//     {
//         for(int i = inputLogSize-1; i >= 0 ; i--)
//         {
//             powI = 1 << i;
//             powI1 = 1 << (i+1);
//             #pragma omp for schedule(auto) nowait
//             for(int j = 0; j < bitVectorSize-1; j+=powI1)
//             {
//                 temp = bitSumVector[j + powI - 1];
//                 bitSumVector[j + powI - 1] = bitSumVector[j + powI1 - 1];
//                 bitSumVector[j + powI1 - 1] = temp + bitSumVector[j + powI1 - 1];
//             }
//         }
//     }

//     // for(auto it=bitSumVector.begin(); it != bitSumVector.end(); it++)
//     // {
//     //     std::cout <<  *it << " ";
//     // }
//     // std::cout << std::endl << "No of Primes : " << bitSumVector.back() << std::endl;

//     out.resize(bitSumVector.back());

//     // std::cout << "Out Size : " << out.size() << std::endl << std::endl;

//     #pragma omp parallel shared(in,out)
//     {
//         #pragma omp for schedule(auto) nowait
//         for(int i = 0; i <= bitVectorSize-1; i++)
//         {
//             if(bitVector[i] == 1)
//             {
//                 out[bitSumVector[i]] = in[i];
//             }
//         }
        
//         #pragma omp single 
//         if(paddedSize == 0 && bitVector[bitVectorSize - 1] == 1)
//         {
//             out.push_back(in[inputSize - 1]);
//         }
//     }

//     // for(auto it=out.begin(); it != out.end(); it++)
//     // {
//     //     std::cout <<  *it << " ";
//     // }
//     // std::cout <<  std::endl;
// }

/* SOLUTION 4: Mask in Parallel : Prefix In Serial : Compute out in Parallel*/

// template <typename T, typename Pred>
// void omp_extract_if(const std::vector<T>& in, std::vector<T>& out, Pred pred) {
//     std::vector<int> input;
//     for(int i=0; i < 10000; i++)
//     {
//         input.push_back(i);
//     }

//     int inputSize = input.size();
//     int inputLogSize = std::ceil(std::log2(inputSize));
//     int bitVectorSize = 1 << inputLogSize;
//     int paddedSize = bitVectorSize - inputSize;
//     std::vector<int> bitVector(bitVectorSize);    
//     int numberOfThreads = 0;
//     int chunkSize = 0;
//     #pragma omp parallel
//     {
//         numberOfThreads = omp_get_num_threads();
//         chunkSize = inputSize/numberOfThreads;
//     }

//     #pragma omp parallel for shared(bitVector)
//     for(int i = 0; i < inputSize/chunkSize; i++)
//     {
//         for(int j = chunkSize*i; j < (chunkSize*i)+chunkSize; j++ )
//         {
//             if(pred(input[j]))
//                 bitVector[j] = 1;
//             else
//                 bitVector[j] = 0;

//         }
//     }

//     std::vector<int> bitSumVector = bitVector;
    
//     int powI = 0;
//     int powI1 = 0;
//     int temp = 0;

//     for(int i = 0; i < inputLogSize-1; i++)
//     {
//         powI = 1 << i;
//         powI1 = 1 << (i+1);
//         for(int j = 0; j < bitVectorSize-1; j += powI1)
//         {
//             bitSumVector[j + powI1 - 1] = 
//                 bitSumVector[j +  powI - 1 ] + bitSumVector[j + powI1 - 1];
//         }
//     }

//     bitSumVector[bitVectorSize - 1] = 0;

//     for(int i = inputLogSize-1; i >= 0 ; i--)
//     {
//         powI = 1 << i;
//         powI1 = 1 << (i+1);
//         for(int j = 0; j < bitVectorSize-1; j+=powI1)
//         {
//             temp = bitSumVector[j + powI - 1];
//             bitSumVector[j + powI - 1] = bitSumVector[j + powI1 - 1];
//             bitSumVector[j + powI1 - 1] = temp + bitSumVector[j + powI1 - 1];
//         }
//     }

//     out.resize(bitSumVector.back(),-1);

//     chunkSize = bitVectorSize/numberOfThreads;

//     #pragma omp parallel for shared(bitVector,out)
//     for(int i = 0; i < bitVectorSize/chunkSize; i++)
//     {
//         for(int j = chunkSize*i; j < (chunkSize*i)+chunkSize; j++ )
//         {
//             if(bitVector[j] == 1)
//             {
//                 out[bitSumVector[j]] = input[j];
//             }
//         }
//     }

//     if(paddedSize == 0 && bitVector[bitVectorSize - 1] == 1)
//     {
//         out.push_back(input[inputSize - 1]);
//     }

//     std::cout << std::endl << "No of Primes : " << bitSumVector.back() << std::endl;

//     for(auto it=out.begin(); it != out.end(); it++)
//     {
//         std::cout <<  *it << " ";
//     }
//     std::cout <<  std::endl;
// }

template <typename T, typename Pred>
void omp_extract_if(const std::vector<T>& in, std::vector<T>& out, Pred pred) {
    int inputSize = in.size();
    std::vector<int> bitVector(inputSize);    
    int numberOfThreads = 0;
    int chunkSize = 0;
    #pragma omp parallel
    {
        numberOfThreads = omp_get_num_threads();
        chunkSize = inputSize/numberOfThreads;
    }

    #pragma omp parallel for shared(bitVector)
    for(int i = 0; i < inputSize/chunkSize; i++)
    {
        for(int j = chunkSize*i; j < (chunkSize*i)+chunkSize; j++ )
        {
            if(pred(in[j]))
                bitVector[j] = 1;
            else
                bitVector[j] = 0;

        }
    }

    std::vector<int> bitSumVector = bitVector;
    
    for(int i=1; i < inputSize; i++)
    {
        bitSumVector[i] = bitSumVector[i] + bitSumVector[i-1];
    }

    out.resize(bitSumVector.back(),-1);

    #pragma omp parallel for shared(bitVector,out)
    for(int i = 0; i < inputSize/chunkSize; i++)
    {
        for(int j = chunkSize*i; j < (chunkSize*i)+chunkSize; j++ )
        {
            if(bitVector[j] == 1)
            {
                out[bitSumVector[j-1]] = in[j];
            }
        }
    }
}

#endif // A0_HPP
