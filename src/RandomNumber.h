// Random number generator adapted from Numerical Recipes
#include <iostream>

typedef unsigned long long int Ullong;
typedef unsigned int Uint;
typedef double Doub;

class RandomNumber {
public:
    // Constructor. Call with any integer seed (except value of v above).
    RandomNumber( ) :
    v(4101842887655102017LL), 
    w(1),
    j( time(0) ) {      // Seed with current time (random seed)     
        u = j ^ v; 
        getLongInt(); 
        v = u; 
        getLongInt(); 
        w = v; 
        getLongInt(); 
    }  

    // Return 64-bit random integer. 
    inline Ullong getLongInt() {     
        u = u * 2862933555777941757LL + 7046029254386353087LL; 
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8; 
        w = 4294957665U*(w & 0xffffffff) + (w >> 32); 
        Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4; 
        return (x + v) ^ w; 
    } 

    // Return random double-precision ﬂoating value in the range 0. to 1.
    inline Doub getDouble() { return 5.42101086242752217E-20 * getLongInt(); } 

    // Return 32-bit random integer.  
    inline Uint getInt() { return (Uint)getLongInt(); } 

private:
    Ullong u, v, w;
    Ullong j;           // initial seed
};
