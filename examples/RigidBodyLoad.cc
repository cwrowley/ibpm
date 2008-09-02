#include <iostream>
#include "ibpm.h"

using namespace ibpm;
using namespace std;

int main() {
    cout << "Type commands to define a RigidBody:" << endl;
    RigidBody body;
    body.load( cin );

    int nPoints = body.getNumPoints();
    BoundaryVector pts = body.getPoints();
    cout << "The body [" << body.getName() << "] has " << nPoints << " points:" << endl;
    cout << pts;
    return 0;
}