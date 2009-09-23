int nx = 100;
int ny = 200;
int ngrid = 3;
double length = 10;
double xOffset = -5;
double yOffset = -5;
Grid grid( nx, ny, ngrid, length, xOffset, yOffset );
cout << "Grid spacing is " << grid.Dx() << endl;
