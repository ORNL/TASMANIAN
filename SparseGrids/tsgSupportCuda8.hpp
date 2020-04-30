#ifndef __TASMANIAN_SUPPORT_CUDA8_HPP
#define __TASMANIAN_SUPPORT_CUDA8_HPP

#define Tasmanian_CUDA8_COMPAT

/*
 * The nvcc compiler for CUDA 8 does not support C++ 2014 even if the host compiler supports the standard.
 * The cuda-specific files are deliberately kept clean from 2014 features, but the common-glue file is
 * tsgAcceleratedDataStructures.hpp which must be visible and compatible with both modes.
 * This adds some functionality that is otherwise missing from nvcc.
 */
namespace std{ // adding missing c++11 functionality
    template< bool B, class T = void >
    using enable_if_t = typename enable_if<B,T>::type;

    #ifdef __NVCC__
    template<class T, class U = T> T exchange(T& x, U&& new_x){
        T old_x = std::move(x);
        x = std::forward<U>(new_x);
        return old_x;
    }

    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args){
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
    #endif
}


#endif
